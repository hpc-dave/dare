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
    static_assert(PM::dimension < 4, "Not equipped for higher dimensions");
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
                               num_tsteps_momentum,
                               bc_args...);
    }
    pm->InitializeContinuity("pressure",
                             grid->GetRepresentation(opt),
                             2,
                             bc_args...);
}

template <typename PM, dare::NaturalNumber Direction, std::integral LO>
dare::CenterMatrixStencil<typename PM::GridType, typename PM::SC, 1>
free_pm_ddt_Cartesian(PM* pm,
            Direction dir,
            const typename PM::GridType::Representation& g_r,
            LO o_loc) {
    static_assert(std::is_same_v<LO, typename PM::GridType::LocalOrdinalType>);
    using GridType = typename PM::GridType;
    using DDT = dare::DDT<GridType>;

    DDT ddt(g_r, o_loc, pm->GetTimeStepSize());

    return ddt(pm->GetPorosity(), pm->GetDensity(), pm->GetMomentum(dir)->GetField());
}

template <typename PM, dare::NaturalNumber Direction, std::integral LO>
dare::CenterMatrixStencil<typename PM::GridType, typename PM::SC, 1>
free_pm_convection_Cartesian(
    PM* pm,
    Direction dir,
    const typename PM::GridType::Representation& g_r,
    LO o_loc,
    const typename PM::PorosityVariableType epsilon,
    const typename PM::DensityVariableType rho,
    const dare::Vector<PM::dimension, const dare::GridVector<typename PM::GridType, typename PM::SC, 1>*>& velocities) {
    static_assert(std::is_same_v<LO, typename PM::GridType::LocalOrdinalType>);
    static_assert(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>);
    using GridType = typename PM::GridType;
    using SC = typename PM::SC;
    using Divergence = dare::Divergence<GridType, typename PM::ConvectiveTimeSchemeType>;
    using FluxLimiter = typename PM::TVDScheme;
    using TVD = dare::TVD<GridType, SC, FluxLimiter>;

    Divergence div(g_r, o_loc);
    TVD tvd(g_r, o_loc, velocities);
    return div(tvd(epsilon, rho, pm->GetMomentum(dir)->GetField()->GetDataVector(1)));
}

template <typename PM, dare::NaturalNumber Direction>
dare::CenterMatrixStencil<typename PM::GridType, typename PM::SC, 1> free_pm_pressure_force_Cartesian(
    PM* pm,
    Direction direction,
    const typename PM::Index& ind,
    const typename PM::PorosityVariableType epsilon) {
    static_assert(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>, "Inconsistent dimensions");
    using PorosityType = std::remove_cv_t<std::remove_pointer_t<decltype(epsilon)>>;
    using SC = typename PM::SC;
    using Index = typename PM::Index;
    Index ind_nb(ind);
    ind_nb[direction] -= 1;
    SC eps{1.};
    if constexpr (dare::is_field_v<PorosityType>)
        eps = 0.5 * (epsilon->GetDataVector().At(ind, 0) + epsilon->GetDataVector().At(ind_nb, 0));
    else if constexpr(!dare::is_none_v<PorosityType>)
        eps = epsilon;

    SC delta_p = pm->GetContinuity()->GetPressure()->GetDataVector().At(ind, 0)
                 - pm->GetContinuity()->GetPressure()->GetDataVector().At(ind_nb, 0);
    SC dV = pm->GetContinuity()->GetGridRepresentation()->GetCellVolume();
    SC dx = pm->GetContinuity()->GetGridRepresentation()->GetDistances()[direction];
    dare::CenterMatrixStencil<typename PM::GridType, SC, 1> stencil;
    stencil.GetRhs(0) = -eps * delta_p / dx * dV;
    return stencil;
}

template <typename PM, dare::NaturalNumber Direction>
dare::CenterMatrixStencil<typename PM::GridType, typename PM::SC, 1>
free_pm_viscious_stress_Cartesian(
    PM* pm,
    Direction,
    const typename PM::GridType::Representation& grep,
    const typename PM::Index& ind,
    const dare::FaceValueStencil<typename PM::GridType, typename PM::SC, 1>& eps_mu_f,
    const dare::Vector<PM::dimension, const dare::GridVector<typename PM::GridType, typename PM::SC, 1>*>& v) {
    static const std::size_t dim = PM::dimension;
    static const std::size_t dir = Direction::value;
    using SC = typename PM::SC;
    using LO = typename PM::LO;
    using GridType = typename PM::GridType;
    using Treatment = typename PM::ViscousStressTreatment;
    using FMStencil = dare::FaceMatrixStencil<GridType, SC, 1>;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CMStencil = dare::CenterMatrixStencil<GridType, typename PM::SC, 1>;
    using CNB = dare::CartesianNeighbor;
    using Index = typename PM::Index;
    using Divergence = dare::Divergence<GridType, typename PM::ConvectiveTimeSchemeType>;
    static_assert(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>, "Inconsistent dimensions");
    static_assert(dim < 4, "limited to 3 dimensions");

    if constexpr(dim < 2) {
        return CMStencil{};
    } else {
        auto dn_r = 1. / grep.GetDistances();

        // implicit component
        FVStencil coef_faces;

        // here we do it really verbose, all the other versions don't improve readability
        CNB f_low = ToFace<dir * 2>();     // lower face in momentum direction
        CNB f_up = ToFace<dir * 2 + 1>();  // upper face in momentum direction
        // if constexpr (dare::is_pm_dijkhuizen_stress_tensor_v<Treatment>) {
        if constexpr (dare::PMDijkhuizenStressTreatment<Treatment>) {
            for (auto face : grep.GetFaces())
                coef_faces(face, 0) = dare::ToNormal(face) * eps_mu_f(face, 0) * dn_r[dare::ToFace(face) / 2];
        } else if constexpr (dare::PMDefaultStressTreatment<Treatment>) {
            coef_faces(f_low, 0) = -eps_mu_f(f_low, 0) * dn_r[dir];
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
                tau_ex(CNB::SOUTH, 0) = (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[1];
                ind_low.j() += 1;
                ind_up.j() += 1;
                tau_ex(CNB::NORTH, 0) = (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[1];
                if constexpr(dim > 2) {
                // du/dz
                    ind_low = ind_up = ind;
                    ind_low.k() -= 1;
                    tau_ex(CNB::BOTTOM, 0) = (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[2];
                    ind_low.k() += 1;
                    ind_up.k() += 1;
                    tau_ex(CNB::TOP, 0) = (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[2];
                }
                ind_low = ind_up = ind;
            }
            ind_low.i() -= 1;
            // dv/dx
            tau_ex(CNB::SOUTH, 0) += (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[0];  // check sign
            // dw/dx
            if constexpr (dim > 2)
                tau_ex(CNB::BOTTOM, 0) += (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[0];  // check sign
            ind_low.j() += 1;
            ind_up.j() += 1;
            tau_ex(CNB::NORTH, 0) += (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[0];  // check sign
            if constexpr (dim > 2) {
                ind_low.j() -= 1;
                ind_up.j() -= 1;
                ind_low.k() += 1;
                ind_up.k() += 1;
                tau_ex(CNB::TOP, 0) += (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[0];  // check sign
            }
        }
        if constexpr (dir == 1) {
            // y-direction
            if constexpr (dare::PMDefaultStressTreatment<Treatment>) {
                // dv/dx
                ind_low.i() -= 1;
                tau_ex(CNB::WEST, 0) = (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[0];
                ind_low.i() += 1;
                ind_up.i() += 1;
                tau_ex(CNB::EAST, 0) = (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[0];
                if constexpr (dim > 2) {
                    // dv/dz
                    ind_low = ind_up = ind;
                    ind_low.k() -= 1;
                    tau_ex(CNB::BOTTOM, 0) = (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[2];
                    ind_low.k() += 1;
                    ind_up.k() += 1;
                    tau_ex(CNB::TOP, 0) = (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[2];
                }
                ind_low = ind_up = ind;
            }
            ind_low.j() -= 1;
            // du/dy
            tau_ex(CNB::WEST, 0) += (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[1];
            // dw/dy
            if constexpr (dim > 2)
                tau_ex(CNB::BOTTOM, 0) += (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[1];
            ind_low.i() += 1;
            ind_up.i() += 1;
            tau_ex(CNB::EAST, 0) += (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[1];
            if constexpr (dim > 2) {
                ind_low.i() -= 1;
                ind_up.i() -= 1;
                ind_low.k() += 1;
                ind_up.k() += 1;
                tau_ex(CNB::TOP, 0) += (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[1];
            }
        } if constexpr(dir == 2) {
            // z-direction
            if constexpr (dare::PMDefaultStressTreatment<Treatment>) {
                // dw/dx
                ind_low.i() -= 1;
                tau_ex(CNB::WEST, 0) = (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[0];
                ind_low.i() += 1;
                ind_up.i() += 1;
                tau_ex(CNB::EAST, 0) = (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[0];
                // dw/dy
                ind_low = ind_up = ind;
                ind_low.j() -= 1;
                tau_ex(CNB::SOUTH, 0) = (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[1];
                ind_low.j() += 1;
                ind_up.j() += 1;
                tau_ex(CNB::NORTH, 0) = (v[2]->At(ind_up, 0) - v[2]->At(ind_low, 0)) * dn_r[1];
                ind_low = ind_up = ind;
            }
            ind_low.k() -= 1;
            // du/dz
            tau_ex(CNB::WEST, 0) += (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[2];  // check sign
            // dv/dz
            tau_ex(CNB::SOUTH, 0) += (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[2];  // check sign
            ind_low.i() += 1;
            ind_up.i() += 1;
            tau_ex(CNB::EAST, 0) += (v[0]->At(ind_up, 0) - v[0]->At(ind_low, 0)) * dn_r[2];  // check sign
            ind_low.i() -= 1;
            ind_up.i() -= 1;
            ind_low.j() += 1;
            ind_up.j() += 1;
            tau_ex(CNB::NORTH, 0) += (v[1]->At(ind_up, 0) - v[1]->At(ind_low, 0)) * dn_r[2];  // check sign
        }
        tau_ex *= eps_mu_f;

        LO o_loc = grep.MapIndexToOrdinalLocalInternal(grep.MapLocalToInternal(ind));
        Divergence div(grep, o_loc);

        dare::CenterMatrixStencil<GridType, SC, 1> s = div(tau_im);
        s.GetRhs() = div(tau_ex);
        return s;
    }
}

template <typename PM, NaturalNumber Direction>
dare::CenterMatrixStencil<typename PM::GridType, typename PM::SC, 1>
free_pm_explicit_force_Cartesian(PM* pm,
                            Direction dir,
                            const typename PM::MomentumType& m,
                            const typename PM::Index& ind) {
    static_assert(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>, "Inconsistent dimensions");
    using VType = std::remove_cv_t<std::remove_pointer_t<typename PM::ExplicitForceVariableType>>;
    using SC = typename PM::SC;
    dare::CenterMatrixStencil<typename PM::GridType, SC, 1> s;
    if constexpr (is_none_v<VType>) {
        // Do nothing
    } else {
        static_assert(FieldType<VType>, "Can only work with fields!");
        SC v{0.};
        for (auto& it : m.GetCustomMember().beta_ex) {
            // Interpolate! Or test at least! How about a field wrapper?
#ifndef DARE_NDEBUG
            auto opt_f = it->GetGridRepresentation().GetOptions();
            auto opt_m = pm->GetMomentum(dir)->GetField()->GetGridRepresentation().GetOptions();
            if (opt_f != opt_m) {
                ERROR << "At the moment, the explicit forces provided to the momentum equations "
                      << "need to be located on the cell center of the staggered grid" << ERROR_CLOSE;
                pm->GetExecutionManager()->Terminate(__func__, "Invalid option combination");
            }
#endif
            v += it->GetDataVector().At(ind, 0);
        }
        v *= pm->GetMomentum(dir)->GetField()->GetGridRepresentation().GetCellVolume();
        s.SetRHS(0, v);
    }
    return s;
}

template <typename PM, dare::NaturalNumber Direction>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
void free_pm_build_momentum(PM* pm, Direction direction) {
    static_assert(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>, "Inconsistent dimensions");
    static const std::size_t dir = Direction::value;
    using GridType = typename PM::GridType;
    using LO = typename GridType::LocalOrdinalType;
    using SC = typename GridType::ScalarType;
    using DensityType = typename PM::DensityVariableType;
    using ViscosityType = typename PM::ViscosityVariableType;
    using PorosityType = typename PM::PorosityVariableType;
    using IndexLocal = typename GridType::Index;
    using GridVectorType = dare::GridVector<GridType, SC, 1>;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    dare::Vector<PM::dimension, const GridVectorType*> velocities;
    for (std::size_t d{0}; d < PM::dimension; d++) {
        velocities[d] = &pm->GetMomentum(d)->GetField()->GetDataVector(1);
    }

    auto BuildStrategy = [=](auto mblock) {
        const typename GridType::Representation* g_r{mblock->GetRepresentation()};
        LO o_loc{mblock->GetLocalOrdinal()};  // this refers to the internal one without ghost/halo cells
        IndexLocal ind{mblock->GetIndex()};

        const DensityType rho{pm->GetDensity()};
        const ViscosityType mu{pm->GetViscosity()};
        const PorosityType epsilon{pm->GetPorosity()};

        // FVStencil rho_f = dare::InterpolateToFaceStencil(*g_r, ind, rho);
        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_r, ind, mu);
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_r, ind, epsilon);
        // accumulation
        (*mblock) += free_pm_ddt_Cartesian(pm, direction, *g_r, o_loc);

        // advection
        (*mblock) += free_pm_convection_Cartesian(pm, *g_r, o_loc, epsilon, rho, velocities);

        // pressure force
        (*mblock) += free_pm_pressure_force_Cartesian(pm, direction, ind, epsilon);

        // viscous stress
        (*mblock) += free_pm_viscious_stress_Cartesian(pm, direction, *g_r, ind, epsilon_f * mu_f, velocities);

        // explicit forcing
        (*mblock) += free_pm_explicit_force_Cartesian(pm, direction, *pm->GetMomentum(dir), ind);

        // Boundary conditions and normalization
    };

    pm->GetMomentum(dir)->Build(BuildStrategy);
}

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_CARTESIAN_H_
