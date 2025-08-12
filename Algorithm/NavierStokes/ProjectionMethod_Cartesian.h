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

template <typename PM, std::size_t Dim, typename BCcont, typename... Args>
void free_pm_initialize(PM* pm, dare::Cartesian<Dim>* grid, BCcont bc_cont, Args... bc_mom) {
    static_assert(PM::dimension == Dim, "The projection method and grid do not have the same dimension!");  // NOLINT
    static_assert(PM::dimension < 4, "Not equipped for higher dimensions");
    static_assert(sizeof...(bc_mom) == Dim, "Inconsistent boundary conditions for the momentum provided");
    static const std::size_t num_tsteps_momentum = PM::num_tsteps_momentum;
    auto tuple_bcmom = std::forward_as_tuple(bc_mom...);
    std::string m_names[] = {"u", "v", "w"};
    typename PM::GridType::Options opt;
    for (auto& o : opt)
        o = 0.;

    pm->InitializeContinuity("pressure",
                             grid->GetRepresentation(opt),
                             2,
                             bc_cont);

    if constexpr(Dim > 0) {
        auto opt_loc = opt;
        opt_loc[0] = 1;
        pm->InitializeMomentum(0,
                               m_names[0],
                               grid->GetRepresentation(opt_loc),
                               num_tsteps_momentum,
                               std::get<0>(tuple_bcmom));
    }
    if constexpr (Dim > 1) {
        auto opt_loc = opt;
        opt_loc[1] = 1;
        pm->InitializeMomentum(1,
                               m_names[1],
                               grid->GetRepresentation(opt_loc),
                               num_tsteps_momentum,
                               std::get<1>(tuple_bcmom));
    }
    if constexpr (Dim > 2) {
        auto opt_loc = opt;
        opt_loc[2] = 1;
        pm->InitializeMomentum(2,
                               m_names[2],
                               grid->GetRepresentation(opt_loc),
                               num_tsteps_momentum,
                               std::get<2>(tuple_bcmom));
    }
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
        for (auto& it : pm->GetMomentum(dir)->GetCustomMember()->beta_ex) {
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
    using CNB = dare::CartesianNeighbor;

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
        (*mblock) += free_pm_convection_Cartesian(pm, direction, *g_r, o_loc, epsilon, rho, velocities);

        // pressure force
        (*mblock) += free_pm_pressure_force_Cartesian(pm, direction, ind, epsilon);

        // viscous stress
        (*mblock) += free_pm_viscious_stress_Cartesian(pm, direction, *g_r, ind, epsilon_f * mu_f, velocities);

        // explicit forcing
        (*mblock) += free_pm_explicit_force_Cartesian(pm, direction, ind);

        // set initial guess
        free_pm_momentum_initialguess(pm, mblock->template Get<CNB::CENTER>(0, 0), mblock);

        // Boundary conditions
        free_pm_momentum_apply_boundary_conditions(pm, *pm->GetMomentum(direction), mblock);

        // normalize
        free_pm_apply_normalizer(pm, pm->GetMomentum(direction)->GetCustomMember()->normalizer, mblock);
    };

    pm->GetMomentum(dir)->Build(BuildStrategy);
}

template <typename PM>
typename PM::SC free_pm_defect_compressible_Cartesian(
    PM* pm,
    typename PM::LO ordinal_internal,
    const typename PM::Index& ind,
    const dare::Vector<PM::dimension, const dare::GridVector<typename PM::GridType, typename PM::SC, 1>*>& v,
    typename PM::DensityVariableType rho,
    typename PM::PorosityVariableType epsilon) {
    using Divergence = dare::Divergence<typename PM::GridType, dare::EULER_BACKWARD>;
    using TVD = dare::TVD<typename PM::GridType, typename PM::SC, typename PM::TVDScheme>;
    using FVStencil = dare::FaceValueStencil<typename PM::GridType, typename PM::SC, 1>;
    using DensityType = std::remove_cv_t<std::remove_pointer_t<typename PM::DensityVariableType>>;
    using PorosityType = std::remove_cv_t<std::remove_pointer_t<typename PM::PorosityVariableType>>;
    auto g_s = &pm->GetContinuity()->GetField()->GetGridRepresentation();
    TVD tvd(pm->GetContinuity()->GetField()->GetGridRepresentation(), ordinal_internal, v);
    Divergence div(pm->GetContinuity()->GetField()->GetGridRepresentation(), ordinal_internal);
    FVStencil fluxes = tvd.GetVelocityStencil();
    if constexpr (dare::is_field_v<DensityType>)
        fluxes *= tvd.Interpolate(dare::convert_to_ptr(rho)->GetDataVector());
    else if constexpr(!dare::is_none_v<DensityType>)
        fluxes *= tvd.Interpolate(rho);
    if constexpr (dare::is_field_v<PorosityType>)
        fluxes *= tvd.Interpolate(dare::convert_to_ptr(epsilon)->GetDataVector());
    else if constexpr (!dare::is_none_v<PorosityType>)
        fluxes *= tvd.Interpolate(epsilon);
    typename PM::SC defect = div(fluxes)[0];
    typename PM::SC eps_rho{1.}, eps_rho_old{1.};
    if constexpr(dare::is_field_v<DensityType>) {
        eps_rho *= pm->GetDensity()->GetDataVector().At(ind, 0);
        eps_rho_old *= pm->GetDensity()->GetDataVector(1).At(ind, 0);
    } else if constexpr (std::is_arithmetic_v<DensityType>) {
        eps_rho *= pm->GetDensity();
        eps_rho_old *= pm->GetDensity();
    } else if constexpr (!dare::is_none_v<DensityType>) {
        static_assert(dare::always_false<DensityType>, "no way defined to get the variable");
    }
    if constexpr (dare::is_field_v<PorosityType>) {
        eps_rho *= pm->GetPorosity()->GetDataVector().At(ind, 0);
        eps_rho_old *= pm->GetPorosity()->GetDataVector(1).At(ind, 0);
    } else if constexpr (std::is_arithmetic_v<PorosityType>) {
        eps_rho *= pm->GetPorosity();
        eps_rho_old *= pm->GetPorosity();
    } else if constexpr (!dare::is_none_v<PorosityType>) {
        static_assert(dare::always_false<PorosityType>, "no way defined to get the variable");
    }
    defect += (eps_rho - eps_rho_old) * g_s->GetCellVolume() / pm->GetTimeStepSize();
    return defect;
}

template <typename PM, typename PorosityVariableType>
typename PM::SC free_pm_defect_incompressible_Cartesian(
    PM* pm,
    const typename PM::LO ordinal_internal,
    const typename PM::Index& ind,
    const dare::Vector<PM::dimension, const dare::GridVector<typename PM::GridType, typename PM::SC, 1>*>& v,
    PorosityVariableType epsilon) {
    using Divergence = dare::Divergence<typename PM::GridType, dare::EULER_BACKWARD>;
    using TVD = dare::TVD<typename PM::GridType, typename PM::SC, dare::CDS>;
    using FVStencil = dare::FaceValueStencil<typename PM::GridType, typename PM::SC, 1>;
    auto g_s = &pm->GetContinuity()->GetField()->GetGridRepresentation();
    FVStencil eps_f = dare::InterpolateToFaceStencil(*g_s, ind, epsilon);
    Divergence div(pm->GetContinuity()->GetField()->GetGridRepresentation(), ordinal_internal);
    TVD tvd(pm->GetContinuity()->GetField()->GetGridRepresentation(), ordinal_internal, v);
    typename PM::SC defect = div(eps_f * tvd.GetVelocityStencil())[0];
    return defect;
}

template <typename PM, typename PorosityVariableType>
    requires dare::is_field_v<std::remove_cv_t<std::remove_pointer_t<PorosityVariableType>>>
typename PM::SC free_pm_defect_incompressible_Cartesian(
    PM* pm,
    const typename PM::LO ordinal_internal,
    const typename PM::Index& ind,
    const dare::Vector<PM::dimension, const dare::GridVector<typename PM::GridType, typename PM::SC, 1>*>& v,
    PorosityVariableType epsilon) {
    return free_pm_defect_incompressible_Cartesian(pm,
                                                   ordinal_internal,
                                                   ind,
                                                   v,
                                                   dare::convert_to_ptr(epsilon)->GetDataVector());
}

template <typename PM>
dare::CenterMatrixStencil<typename PM::GridType, typename PM::SC, 1>
free_pm_continuity_Jacobian_Cartesian(
    PM* pm,
    typename PM::LO ordinal_internal,
    const typename PM::Index& ind,
    typename PM::DensityVariableType rho,
    typename PM::PorosityVariableType epsilon) {
    using GridType = typename PM::GridType;
    using SC = PM::SC;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CMStencil = dare::CenterMatrixStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;
    using ForceType = std::remove_cv_t<std::remove_pointer_t<typename PM::ImplicitForceVariableType>>;
    using PorosityType = std::remove_cv_t<std::remove_pointer_t<decltype(epsilon)>>;
    using DensityType = std::remove_cv_t<std::remove_pointer_t<decltype(rho)>>;
    using Divergence = dare::Divergence<GridType, dare::EULER_BACKWARD>;
    using Gradient = dare::Gradient<GridType>;

    auto g_s = &pm->GetContinuity()->GetField()->GetGridRepresentation();
    typename PM::SC dt = pm->GetTimeStepSize();
    Divergence div(*g_s, ordinal_internal);
    Gradient grad(*g_s, ordinal_internal);
    FVStencil eps_f, rho_f;
    if constexpr (dare::is_field_v<DensityType>)
        rho_f = dare::InterpolateToFaceStencil(*g_s, ind, dare::convert_to_ptr(rho)->GetDataVector());
    else
        rho_f = dare::InterpolateToFaceStencil(*g_s, ind, rho);
    if constexpr (dare::is_field_v<PorosityType>)
        eps_f = dare::InterpolateToFaceStencil(*g_s, ind, dare::convert_to_ptr(epsilon)->GetDataVector());
    else
        eps_f = dare::InterpolateToFaceStencil(*g_s, ind, epsilon);
    FVStencil beta_f;
    beta_f.SetAll(0.);
    if constexpr (dare::is_field_v<ForceType>) {
        for (auto b : pm->GetContinuity()->GetCustomMember()->beta_im)
            beta_f += dare::InterpolateToFaceStencil(*g_s, ind, b->GetDataVector());
    } else if constexpr (!dare::is_none_v<ForceType>) {
        static_assert(dare::always_false<ForceType>, "This type is not supported for forces");
    }

    CMStencil s;
    if constexpr (PM::compressible) {
        FVStencil coef_f = -1. * eps_f * dt / (1. + beta_f * dt / (eps_f * rho_f));
        s = div(coef_f, grad(dare::ONE));

        // Main diagonal with density-derivative
        SC dd_dp = pm->GetDensityDerivative(ind);
        SC eps_c{1.};
        if constexpr (dare::is_field_v<PorosityType>) {
            eps_c = epsilon->GetDataVector().At(ind, 0);
        } else if constexpr (std::is_arithmetic_v<PorosityType>) {
            eps_c = epsilon;
        }
        dd_dp *= eps_c;
        dd_dp *= g_s->GetCellVolume() / dt;
        s.GetValue(CNB::CENTER, 0) += dd_dp;
    } else {
        FVStencil coef_f = -1. * eps_f * dt / (rho_f + beta_f / eps_f * dt);
        s = div(coef_f, grad(dare::ONE));
    }

    return s;
}

template <typename PM>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
void free_pm_compute_defect(PM* pm) {
    using LO = typename PM::LO;
    using Index = typename PM::Index;
    using GridVectorType = dare::GridVector<typename PM::GridType, typename PM::SC, 1>;
    using DefectType = GridVectorType;
    auto g_s = &pm->GetPressure()->GetGridRepresentation();
    DefectType* defect = &pm->GetContinuity()->GetDefect()->GetDataVector();
    dare::Vector<PM::dimension, const GridVectorType*> velocities;
    for (std::size_t d{0}; d < PM::dimension; d++) {
        velocities[d] = &pm->GetMomentum(d)->GetField()->GetDataVector(0);
    }

#pragma omp parallel for
    for (LO n = 0; n < g_s->GetNumberLocalCellsInternal(); n++) {
        Index ind_internal{g_s->MapOrdinalToIndexLocalInternal(n)};
        Index ind{g_s->MapInternalToLocal(ind_internal)};
        if constexpr (PM::compressible) {
            defect->At(ind, 0) = free_pm_defect_compressible_Cartesian(pm,
                                                                       n,
                                                                       ind,
                                                                       velocities,
                                                                       pm->GetDensity(),
                                                                       pm->GetPorosity());
        } else {
            defect->At(ind, 0) = free_pm_defect_incompressible_Cartesian(pm,
                                                                         n,
                                                                         ind,
                                                                         velocities,
                                                                         pm->GetPorosity());
        }
    }
}

template <typename PM>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
void free_pm_build_continuity(PM* pm, int iteration) {
    static_assert(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>, "Inconsistent dimensions");
    using GridType = typename PM::GridType;
    using LO = typename GridType::LocalOrdinalType;
    using SC = typename GridType::ScalarType;
    using DensityType = typename PM::DensityVariableType;
    using PorosityType = typename PM::PorosityVariableType;
    using IndexLocal = typename GridType::Index;
    using GridVectorType = dare::GridVector<GridType, SC, 1>;

    dare::Vector<PM::dimension, const GridVectorType*> velocities;
    for (std::size_t d{0}; d < PM::dimension; d++) {
        velocities[d] = &pm->GetMomentum(d)->GetField()->GetDataVector(1);
    }

    if (iteration == 0) {
        auto BuildStrategy = [=](auto mblock) {
            // const typename GridType::Representation* g_r{mblock->GetRepresentation()};
            LO o_loc{mblock->GetLocalOrdinal()};  // this refers to the internal one without ghost/halo cells
            IndexLocal ind{mblock->GetIndex()};

            const DensityType rho{pm->GetDensity()};
            const PorosityType epsilon{pm->GetPorosity()};

            (*mblock) = free_pm_continuity_Jacobian_Cartesian(pm, o_loc, ind, rho, epsilon);
            mblock->GetRhs(0) = -1. * pm->GetContinuity()->GetDefect()->GetDataVector().At(ind, 0);

            // Apply Boundary conditions
            pm->GetContinuity()->GetBoundaryStrategy()->Apply(mblock);
        };
        pm->GetContinuity()->Build(BuildStrategy);
    } else {
        auto BuildStrategy = [=](auto mblock) {
            // const typename GridType::Representation* g_r{mblock->GetRepresentation()};
            // LO o_loc{mblock->GetLocalOrdinal()};  // this refers to the internal one without ghost/halo cells
            IndexLocal ind{mblock->GetIndex()};

            mblock->GetRhs(0) = -1. * pm->GetContinuity()->GetDefect()->GetDataVector().At(ind, 0);

            // Apply Boundary conditions
            pm->GetContinuity()->GetBoundaryStrategy()->Apply(mblock);
        };
        pm->GetContinuity()->UpdateRhs(BuildStrategy);
    }
}

template <typename PM>
    requires(std::is_same_v<typename PM::GridType, dare::Cartesian<PM::dimension>>)
void free_pm_update_velocity(PM* pm, int iteration) {
    using GridType = typename PM::GridType;
    using LO = typename PM::LO;
    using SC = typename PM::SC;
    using Index = typename PM::Index;
    using GridVector = dare::GridVector<GridType, SC, 1>;
    using GridRepresentation = typename GridType::Representation;
    using ImplicitForceType = std::remove_cv_t<std::remove_pointer_t<typename PM::ImplicitForceVariableType>>;
    using CNB = dare::CartesianNeighbor;
    // using TVD = dare::TVD<typename PM::TVDScheme>;

    auto GetBeta = [=](const GridRepresentation& grep, Index ind, CNB face) -> SC {
        SC beta{0.};
        if constexpr (dare::is_field_v<ImplicitForceType>) {
            for (const auto& f : pm->GetContinuity()->GetCustomMember()->beta_im)
                beta += dare::InterpolateToFace(grep, ind, face, f->GetDataVector())[0];
        } else if constexpr (!dare::is_none_v<ImplicitForceType>) {
            static_assert(dare::always_false<ImplicitForceType>, "Cannot deal with provided type");
        }
        return beta;
    };

    const GridVector* dp = &pm->GetContinuity()->GetdP()->GetDataVector();
    const SC dt = pm->GetTimeStepSize();
    const GridRepresentation* grep_p = &pm->GetContinuity()->GetField()->GetGridRepresentation();
    for (std::size_t d{0}; d < pm->GetDimension(); d++) {
        GridVector* v = &pm->GetMomentum(d)->GetField()->GetDataVector();
        const GridRepresentation* grep = &v->GetGridRepresentation();
        SC dn_r = 1. / grep->GetDistances()[d];
        for (LO n_loc{0}; n_loc < grep->GetNumberLocalCellsInternal(); n_loc++) {
            Index ind = grep->MapOrdinalToIndexLocalInternal(n_loc);
            ind = grep->MapInternalToLocal(ind);
            if constexpr(PM::compressible) {
                // Do a good check for
                // - How beneficial is the use of a TVD based interpolation of the involved parameters?
                // - How beneficial is the use of the density of the previous iteration vs the iteration from
                //   the old timestep?
                // TVD tvd(grep, ind);  // Is it beneficial to interpolate here with a TVD scheme?
                SC v_prev = v->At(ind, 0);
                Index ind_lo(ind);
                ind_lo[d] -= 1;
                SC epsilon{0.5 * (pm->GetPorosity(ind_lo, 0) + pm->GetPorosity(ind, 0))};
                SC rho{0.5 * (pm->GetDensity(ind_lo, 0) + pm->GetDensity(ind, 0))};
                SC rho_prev{0.5 * (pm->GetDensityPreviousIteration()->At(ind_lo, 0)
                                    + pm->GetDensityPreviousIteration()->At(ind, 0))};
                SC beta{GetBeta(*grep_p, ind, dare::ToFace(d * 2)) * (iteration == 0)};
                SC dP_hi{dp->At(ind, 0)};
                SC dP_lo{dp->At(ind_lo, 0)};
                SC dP_dx = (dP_hi - dP_lo) * dn_r;
                SC v_new = epsilon * rho_prev / (epsilon * rho + beta * dt) * (v_prev - dt / rho_prev * dP_dx);
                v->At(ind, 0) = v_new;
            } else {
                SC v_prev = v->At(ind, 0);
                Index ind_lo(ind);
                ind_lo[d] -= 1;
                SC epsilon{0.5 * (pm->GetPorosity(ind_lo, 0) + pm->GetPorosity(ind, 0))};
                SC rho{0.5 * (pm->GetDensity(ind_lo, 0) + pm->GetDensity(ind, 0))};
                SC beta{GetBeta(*grep_p, ind, dare::ToFace(d * 2)) * (iteration == 0)};
                SC dP_hi{dp->At(ind, 0)};
                SC dP_lo{dp->At(ind_lo, 0)};
                SC dP_dx = (dP_hi - dP_lo) * dn_r;
                SC v_new = 1./ (1. + beta * dt / (epsilon * rho)) * (v_prev - dt / rho * dP_dx);
                v->At(ind, 0) = v_new;
            }
        }
        pm->GetMomentum(d)->UpdateBoundaries();
    }
}

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_CARTESIAN_H_
