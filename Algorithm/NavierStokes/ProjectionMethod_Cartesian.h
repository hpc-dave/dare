/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_CARTESIAN_H_
#define ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_CARTESIAN_H_

#include <concepts>
#include <string>
#include <type_traits>
#include <utility>

#include "Grid/Cartesian.h"
#include "ProjectionMethod_freefunc.h"

namespace dare {

template <typename PM, std::size_t Dim>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<Dim>>)
void free_compile_time_check(PM*) {
    static_assert(dare::uses_newton_iterations_v<PM::ContinuityIterationsType>,
                  "Cartesian grid right now only uses newton iterations for enforcing continuity");
}
template <typename PM, std::size_t Dim, typename... Args>
void free_pm_initialize(PM* pm, dare::Cartesian<Dim>* grid, Args&&... bc_args) {
    static_assert(PM::dimension == Dim, "The projection method and grid do not have the same dimension!");  // NOLINT
    static_assert(PM::dimension < 3, "Not equipped for higher dimensions");
    static const std::size_t num_tsteps_momentum = PM::num_tsteps_momentum;
    std::string m_names[] = {"u", "v", "w"};
    typename PM::GridType::Options opt;
    for (auto& o : opt)
        o = 0.;

    for (std::size_t d{0}; d < Dim; d++) {
        auto opt_loc = opt;
        opt_loc[d] = 1;
        pm->InitializeMomentum(d,
                               m_names[d],
                               grid->GetRepresentation(opt_loc),
                               grid->GetExecutionManager(),
                               num_tsteps_momentum,
                               bc_args...);
    }
    pm->InitializeContinuity("pressure",
                             grid->GetRepresentation(opt),
                             grid->GetExecutionManager(),
                             2,
                             bc_args...);
}

template <typename PM, dare::NaturalNumber Direction>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
void free_pm_build_momentum(PM* pm, Direction direction) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method");  // NOLINT
    static const std::size_t dir = Direction::value;
    using GridType = typename PM::GridType;
    using CNB = dare::CartesianNeighbor;
    using LO = typename GridType::LocalOrdinalType;
    using SC = typename GridType::ScalarType;
    using DensityType = typename PM::DensityVariableType;
    using ViscosityType = typename PM::ViscosityVariableType;
    using PorosityType = typename PM::ViscosityVariableType;
    using FluxLimiter = typename PM::TVDScheme;
    using IndexLocal = typename GridType::Index;
    using GridVectorType = dare::GridVector<GridType, SC, 1>;
    using DDT = dare::DDT<GridType>;
    using TVD = dare::TVD<GridType, SC, FluxLimiter>;
    using DivergenceVStress = dare::Divergence<GridType, dare::EULER_BACKWARD>;
    using DivergenceAdvection = dare::Divergence<GridType, typename PM::ConvectiveTimeSchemeType>;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    static const CNB f_low[] = {CNB::WEST, CNB::SOUTH, CNB::BOTTOM};
    dare::Vector<PM::dimension, const GridVectorType*> velocities;
    for (std::size_t d{0}; d < PM::dimension; d++) {
        velocities[d] = &pm->GetMomentum(d)->GetField().GetDataVector(1);
    }

    auto BuildStrategy = [=](auto mblock) {
        auto g_r{mblock->GetRepresentation()};
        LO o_loc{mblock->GetLocalOrdinal()};  // this refers to the internal one without ghost/halo cells
        IndexLocal ind{mblock->GetIndex()};

        DDT ddt(*g_r, o_loc, pm->GetTimeStepSize());
        DivergenceAdvection div_a(*g_r, o_loc);
        DivergenceVStress div_v(*g_r, o_loc);
        // check here with is_pointer_v
        TVD tvd(*g_r, o_loc, velocities);

        const DensityType rho{pm->GetDensity()};
        const ViscosityType mu{pm->GetViscosity()};
        const PorosityType epsilon{pm->GetPorosity()};

        // FVStencil rho_f = dare::InterpolateToFaceStencil(*g_r, ind, rho);
        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_r, ind, mu);
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_r, ind, epsilon);
        // accumulation
        (*mblock) = ddt(epsilon, rho, *pm->GetMomentum(dir));

        // advection
        (*mblock) += div_a(tvd.Interpolate(epsilon),
                           tvd.Interpolate(rho),
                           tvd.interpolate(velocities[dir]->GetDataVector()),
                           *velocities[dir]);

        // pressure force
        mblock->GetRhs(0) += pm_pressure_force_Cartesian(pm, direction, ind, epsilon_f);

        // viscous stress
        auto [tau_im, tau_ex] = pm_viscious_stress_Cartesian(pm, direction, *g_r, ind, epsilon_f * mu_f, velocities);
        (*mblock) += div_v(tau_im + tau_ex);

        // explicit forcing
        mblock->GetRhs(0) += pm_explicit_force_Cartesian(pm, *pm->GetMomentum(dir), ind);
    };

    pm->GetMomentum(dir)->Build(BuildStrategy);
}

template<typename PM, dare::NaturalNumber Direction>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
typename PM::SC pm_pressure_force_Cartesian(
    PM* pm,
    Direction,
    const typename PM::IndexLocal& ind,
    const dare::FaceValueStencil<typename PM::GridType, typename PM::SC, 1>& epsilon) {
    static const std::size_t dir = Direction::value;
    using SC = typename PM::SC;
    using Index = typename PM::IndexLocal;
    using CNB = typename dare::CartesianNeighbor;
    Index ind_nb(ind);
    ind_nb[dir] -= 1;
    SC eps{0.};
    if constexpr(dir == 0)
        eps = epsilon.GetValue(CNB::WEST, 0);
    else if constexpr(dir == 1)
        eps = epsilon.GetValue(CNB::SOUTH, 0);
    else if constexpr(dir == 2)
        eps = epsilon.GetValue(CNB::BOTTOM, 0);
    else
        static_assert(dare::always_false<PM>, "ONLY UP TO 3D, STUPID!");

    SC delta_p = pm->GetContinuity()->GetPressure().At(ind) - pm->GetContinuity()->GetPressure().At(ind_nb);
    SC dV = pm->GetContinuity()->GetRepresentation()->GetCellVolume();
    SC dx = pm->GetContinuity()->GetRepresentation()->GetDistances()[dir];
    return -eps * delta_p / dx * dV;
}

template <typename PM, dare::NaturalNumber Direction>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
std::pair<dare::FaceMatrixStencil<typename PM::GridType, typename PM::SC, 1>,
          dare::FaceValueStencil<typename PM::GridType, typename PM::SC, 1>>
pm_viscious_stress_Cartesian(
    PM* pm,
    Direction,
    const typename PM::GridType::Representation& grep,
    const typename PM::IndexLocal& ind,
    const dare::FaceValueStencil<typename PM::GridType, typename PM::SC, 1>& eps_mu_f,
    const dare::Vector<PM::dimension, const dare::GridVector<typename PM::GridType, typename PM::SC, 1>>& v) {
    static const std::size_t dim = PM::dimension;
    static const std::size_t dir = Direction::value;
    using Treatment = typename PM::ViscousStressTreatment;
    using FMStencil = dare::FaceMatrixStencil<typename PM::GridType, typename PM::SC, 1>;
    using FVStencil = dare::FaceValueStencil<typename PM::GridType, typename PM::SC, 1>;
    using CNB = dare::CartesianNeighbor;
    using Index = typename PM::IndexLocal;

    static_assert(dim < 4, "limited to 3 dimensions");

    if constexpr(dim < 2) {
        FMStencil m_empty;
        FVStencil f_empty;
        return std::make_pair(m_empty, f_empty);
    } else {
        auto dn_r = 1. / grep->GetDistances();

        // implicit component
        FVStencil coef_faces;

        // here we do it really verbose, all the other versions don't improve readability
        CNB f_low = ToFace<dir * 2>();     // lower face in momentum direction
        CNB f_up = ToFace<dir * 2 + 1>();  // upper face in momentum direction
        // if constexpr (dare::is_pm_dijkhuizen_stress_tensor_v<Treatment>) {
        if constexpr (dare::PMDijkhuizenStressTreatment<Treatment>) {
            for (auto face : grep.GetFaces())
                coef_faces(face, 0) = eps_mu_f(face, 0) * dn_r[dare::ToFace(face) / 2];
        } else if constexpr (dare::PMDefaultStressTreatment<Treatment>) {
            coef_faces(f_low, 0) = eps_mu_f(f_low, 0) * dn_r[dir];
            coef_faces(f_up, 0) = eps_mu_f(f_up, 0) * dn_r[dir];
        }
        // and again for the main diagonal
        coef_faces(f_low, 0) *= 2.;
        coef_faces(f_up, 0) *= 2.;

        FMStencil tau_im;
        for (auto face : grep.GetFaces()) {
            tau_im.SetValues(face, 0, -coef_faces(face, 0), coef_faces(face, 0));
        }

        // explicit components
        FVStencil tau_ex;
        Index ind_low(ind), ind_up(ind);
        if constexpr (dir == 0) {
            // x-direction
            if constexpr (dare::PMDefaultStressTreatment<Treatment>) {
                // du/dy
                ind_low.j() -= 1;
                tau_ex(CNB::SOUTH) = v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0);
                ind_low.j() += 1;
                ind_up.j() += 1;
                tau_ex(CNB::NORTH) = v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0);
                if constexpr(dim > 2) {
                // du/dz
                    ind_low = ind_up = ind;
                    ind_low.k() -= 1;
                    tau_ex(CNB::BOTTOM) = v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0);
                    ind_low.k() += 1;
                    ind_up.k() += 1;
                    tau_ex(CNB::TOP) = v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0);
                }
                ind_low = ind_up = ind;
            }
            ind_low.i() -= 1;
            // dv/dx
            tau_ex(CNB::SOUTH) += v[1]->At(ind_up, 0) = v[1]->At(ind_low, 0);  // check sign
            // dw/dx
            if constexpr (dim > 2)
                tau_ex(CNB::BOTTOM) += v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);  // check sign
            ind_low.i() += 1;
            ind_up.i() += 1;
            tau_ex(CNB::NORTH) += v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0);  // check sign
            if constexpr(dim > 2)
                tau_ex(CNB::TOP) += v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);  // check sign
        }
        if constexpr (dir == 1) {
            // y-direction
            if constexpr (dare::PMDefaultStressTreatment<Treatment>) {
                // dv/dx
                ind_low.i() -= 1;
                tau_ex(CNB::WEST) = v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0);
                ind_low.i() += 1;
                ind_up.i() += 1;
                tau_ex(CNB::EAST) = v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0);
                if constexpr (dim > 2) {
                    // dv/dz
                    ind_low = ind_up = ind;
                    ind_low.k() -= 1;
                    tau_ex(CNB::BOTTOM) = v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0);
                    ind_low.k() += 1;
                    ind_up.k() += 1;
                    tau_ex(CNB::TOP) = v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0);
                }
                ind_low = ind_up = ind;
            }
            ind_low.j() -= 1;
            // du/dy
            tau_ex(CNB::WEST) += v[0]->At(ind_up, 0) = v[0]->At(ind_low, 0);  // check sign
            // dw/dy
            if constexpr (dim > 2)
                tau_ex(CNB::BOTTOM) += v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);  // check sign
            ind_low.j() += 1;
            ind_up.j() += 1;
            tau_ex(CNB::EAST) += v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0);  // check sign
            if constexpr (dim > 2)
                tau_ex(CNB::TOP) += v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);  // check sign
        } if constexpr(dir > 2) {
            // z-direction
            if constexpr (dare::PMDefaultStressTreatment<Treatment>) {
                // dw/dx
                ind_low.i() -= 1;
                tau_ex(CNB::WEST) = v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);
                ind_low.i() += 1;
                ind_up.i() += 1;
                tau_ex(CNB::EAST) = v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);
                // dw/dy
                ind_low = ind_up = ind;
                ind_low.j() -= 1;
                tau_ex(CNB::SOUTH) = v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);
                ind_low.j() += 1;
                ind_up.j() += 1;
                tau_ex(CNB::NORTH) = v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0);
                ind_low = ind_up = ind;
            }
            ind_low.k() -= 1;
            // du/dz
            tau_ex(CNB::WEST) += v[0]->At(ind_up, 0) = v[0]->At(ind_low, 0);  // check sign
            // dv/dz
            tau_ex(CNB::SOUTH) += v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0);  // check sign
            ind_low.k() += 1;
            ind_up.k() += 1;
            tau_ex(CNB::EAST) += v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0);  // check sign
            tau_ex(CNB::NORTH) += v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0);  // check sign
        }
        for (auto face : grep.GetFaces()) {
            tau_ex *= dn_r[dare::ToFace(face) / 2];
        }
        tau_ex *= eps_mu_f;
        return std::make_pair(tau_im, tau_ex);
    }
}

template<typename PM, NaturalNumber Direction>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
typename PM::SC pm_explicit_force_Cartesian(PM* pm, Direction,
                                            const typename PM::MomentumType& m,
                                            typename PM::IndexLocal ind) {
    // static const std::size_t dir = Direction::value;
    using VType = std::remove_cv_t<std::remove_pointer_t<typename PM::ExplicitForceVariableType>>;
    using SC = typename PM::SC;
    if constexpr (is_none_v<VType>) {
        return 0.;
    } else {
        static_assert(FieldType<VType>, "Can only work with fields!");
        SC v{0.};
        for (auto& it : m.GetCustomMember().beta_ex) {
            v += it->GetDataVector().At(ind, 0);
        }
        return v;
    }
}

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_CARTESIAN_H_
