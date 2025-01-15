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

#ifndef ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_
#define ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_

#include <concepts>
#include <array>
#include <memory>
#include <algorithm>
#include <bitset>

#include "PM_Information.h"
#include "PM_details.h"
#include "PM_Momentum.h"
#include "PM_Continuity.h"
#include "Algorithm/AlgorithmTraits.h"
#include "Data/Field.h"
#include "Equations/FluxLimiter.h"
#include "Utilities/Errors.h"
#include "Equations/TimeDiscretizationSchemes.h"
#include "Utilities/InitializationTracker.h"

namespace dare::algorithm {

template <typename Grid>
struct PMPropertyInfoDefault {
    using FieldType = dare::Data::Field<Grid, typename Grid::ScalarType, 1>;
    using density = PMDensityInfo<FieldType>;
    using viscosity = PMViscosityInfo<FieldType>;
    using porosity = PMPorosityInfo<dare::utils::None>;
    using implicit_force = PMImplicitForceInfo<dare::utils::None>;
    using explicit_force = PMExplicitForceInfo<dare::utils::None>;
    using compressible = PMCompressibleInfo<false>;
};

struct PMNumericalInfoDefault {
    using tvd = PMTVDInfo<dare::Matrix::MINMOD>;
    using viscous_stress = PMViscousStressInfo<PMDefaultStressTensor>;
    using momentum_iterations = PMMomentumIterationInfo<dare::algorithm::FixedPoint>;
    using continuity_iterations = PMContinuityIterationInfo<dare::algorithm::Newton>;
    using time_scheme_convective = PMTimeSchemeConvectiveInfo<dare::Matrix::EULER_BACKWARD>;
};

/*!
 * @brief holds all relevant entities to compute the flow field via a two-step projection method
 * @tparam Grid the grid type
 * @tparam AlgorithmInfo compile time information about the algorithm
 * @tparam PropertyInfo compile time information about the properties
 *
 * EXPLAIN THE COMPILE TIME INFORMATION APPROACH
 * Compile time information concerning the properties:
 * - density:
 * - viscosity:
 * - porosity:
 * - implicit_force:
 * - explicit_force:
 *
 *
 */
template <typename Grid, typename BoundaryStrategy, typename PropertyInfo, typename NumericalInfo>
class ProjectionMethod : public dare::utils::InitializationTracker {
public:
    enum {
        rho_init = 0b0000001,
        mu_init = 0b0000010,
        epsilon_init = 0b0000100,
        beta_im_init = 0b0001000,
        beta_ex_init = 0b0010000
    };
    // general types based on the grid
    using GridType = Grid;
    using BoundaryStrategyType = BoundaryStrategy;
    using SC = typename Grid::ScalarType;
    using LO = typename Grid::LocalOrdinalType;
    using GO = typename Grid::GlobalOrdinalType;
    using Index = typename Grid::Index;
    using IndexGlobal = typename Grid::IndexGlobal;
    using FieldType = Data::Field<GridType, SC, 1>;

    // properties determined from the PropertyInfo type
    using PropertyTypeInfo = detail::PMAssembledPropertyInfoWithDefaults<
                                        PropertyInfo,
                                        PMPropertyInfoDefault<Grid>>;
    using DensityInfo = typename PropertyTypeInfo::density;
    using ViscosityInfo = typename PropertyTypeInfo::viscosity;
    using PorosityInfo = typename PropertyTypeInfo::porosity;
    using ImplicitForceInfo = typename PropertyTypeInfo::implicit_force;
    using ExplicitForceInfo = typename PropertyTypeInfo::explicit_force;
    using CompressibilityInfo = typename PropertyTypeInfo::compressible;
    using DensityVariableType = detail::determine_density_variable_type_t<DensityInfo>;
    using ViscosityVariableType = detail::determine_viscosity_variable_type_t<ViscosityInfo>;
    using PorosityVariableType = detail::determine_porosity_variable_type_t<PorosityInfo>;
    using ImplicitForceVariableType = detail::determine_implicit_force_variable_type_t<ImplicitForceInfo>;
    using ExplicitForceVariableType = detail::determine_explicit_force_variable_type_t<ExplicitForceInfo>;
    using ImplicitForceMemberType = detail::determine_implicit_force_member_variable_type_t<ImplicitForceInfo>;
    using ExplicitForceMemberType = detail::determine_explicit_force_member_variable_type_t<ExplicitForceInfo>;
    static const bool compressible = CompressibilityInfo::flag;
    static const std::size_t dimension = Grid::Dimension;

    // algorithm and discretization properties
    using NumericalTypeInfo = detail::PMAssembledNumericalInfoWithDefaults<
                                        NumericalInfo,
                                        PMNumericalInfoDefault>;
    using TVDInfo = typename NumericalTypeInfo::tvd;
    using ViscousStressInfo = typename NumericalTypeInfo::viscous_stress;
    using MomentumIterationInfo = typename NumericalTypeInfo::momentum_iterations;
    using ContinuityIterationInfo = typename NumericalTypeInfo::continuity_iterations;
    using TVDScheme = typename TVDInfo::type;
    using ViscousStressTreatment = typename ViscousStressInfo::type;
    using MomentumIterationType = typename MomentumIterationInfo::type;
    using ContinuityIterationType = typename ContinuityIterationInfo::type;
    using ConvectiveTimeSchemeType = typename ContinuityIterationInfo::type;
    static const std::size_t num_tsteps_momentum = std::max(ConvectiveTimeSchemeType::NUM_TSTEPS + 1, 2);

    struct MomentumMembers{
        ExplicitForceMemberType beta_ex;
    };
    struct ContinuityMembers {
        ExplicitForceMemberType beta_im;
    };
    using MomentumType = PMMomentum<GridType, BoundaryStrategyType, MomentumMembers>;
    using ContinuityType = PMContinuity<GridType, BoundaryStrategyType, ContinuityMembers>;

    ProjectionMethod()
        : ex_man(nullptr),
          status(0), status_finalized(rho_init | mu_init | epsilon_init | beta_im_init | beta_ex_init) {
        if constexpr (dare::utils::is_none_v<PorosityVariableType>)
            status |= epsilon_init;
        if constexpr (dare::utils::is_none_v<ImplicitForceVariableType>)
            status |= beta_im_init;
        if constexpr (dare::utils::is_none_v<ExplicitForceVariableType>)
            status |= beta_ex_init;
    }

    template<typename... Args>
    void Initialize(const GridType& grid, Args&&... bc_args) {
        ex_man = grid.GetExecutionManager();
        free_pm_initialize(this, grid, bc_args...);
        this->Initialize();
    }

    // for access in the free functions
    template <typename... Args>
    void IntializeMomentum(std::size_t dim, Args&&... args) {
        momentum[dim] = std::make_unique<MomentumType>(args...);
    }

    void SolveFlowField() {
        if (!CheckStatus()) {
            ex_man->Terminate(__func__, "Projection method was not fully finalized!");
        }
        free_pm_solve(this);
    }

    constexpr bool IsCompressible() const { return compressible; }
    constexpr std::size_t GetDimension() const { return dimension; }

    std::unique_ptr<MomentumType>& GetMomentum(std::size_t dim) { return momentum[dim]; }
    const std::unique_ptr<MomentumType>& GetMomentum(std::size_t dim) const { return momentum[dim]; }

    ContinuityType* GetContinuity() { return &continuity; }
    const ContinuityType& GetContinuity() const { return continuity; }

    void SetDensity(DensityVariableType d) {
        rho = d;
        status |= rho_init;
    }
    void SetViscosity(ViscosityVariableType v) {
        mu = v;
        status |= mu_init;
    }
    void SetPorosity(PorosityVariableType p) {
        epsilon = p;
        status |= epsilon_init;
    }

    void AddImplicitForce(ImplicitForceVariableType f) {
        // add to continuity
        if (!IsInitialized()) {
            ex_man->Terminate(__func__, "Cannot add force terms prior to initialization");
        }
        if constexpr (!dare::utils::is_none_v<ImplicitForceVariableType>) {
            continuity.GetCustomMember()->beta_im.emplace(f);
            status |= beta_im_init;
        }
        status |= beta_im_init;
    }

    void AddExplicitForce(ExplicitForceVariableType f, std::size_t dim) {
        // add to continuity
        if (!IsInitialized()) {
            ex_man->Terminate(__func__, "Cannot add force terms prior to initialization");
        }
        if constexpr (!dare::utils::is_none_v<ExplicitForceVariableType>) {
            if ((dim < 1) || (dim >= dimension)) {
                ex_man->Terminate(__func__, "Invalid dimension choses for the force");
            }
            momentum[dim]->GetCustomMember()->beta_ex.emplace(f);

            // check if all explicit force members were set
            bool all_init{true};
            for (auto& m : momentum)
                all_init &= m->GetCustomMember()->beta_ex.empty();
            status |= (beta_ex_init & all_init);
        }
    }

    bool CheckStatus() const {
        return (status == status_finalized) && this->IsInitialized();
    }

private:
    dare::mpi::ExecutionManager* ex_man;
    DensityVariableType rho;
    ViscosityVariableType mu;
    PorosityVariableType epsilon;

    ContinuityType continuity;
    std::array<std::unique_ptr<MomentumType>, dimension> momentum;

    char status;
    char status_finalized;
};

}  // namespace dare::algorithm

#include "ProjectionMethod.inl"

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_
