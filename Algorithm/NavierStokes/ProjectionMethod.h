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
#include <utility>

#include "PM_Information.h"
#include "PM_details.h"
#include "PM_Continuity.h"
#include "Algorithm/AlgorithmTraits.h"
#include "Data/Field.h"
#include "Equations/FluxLimiter.h"
#include "Utilities/Errors.h"
#include "Equations/TimeDiscretizationSchemes.h"
#include "Utilities/InitializationTracker.h"
#include "Equations/GenericEquation.h"
#include "ProjectionMethod_freefunc.h"

namespace dare {

template <typename Grid>
struct PMPropertyInfoDefault {
    using FieldType = dare::Field<Grid, typename Grid::ScalarType, 1>;
    using density = PMDensityInfo<FieldType>;
    using viscosity = PMViscosityInfo<FieldType>;
    using porosity = PMPorosityInfo<dare::None>;
    using implicit_force = PMImplicitForceInfo<dare::None>;
    using explicit_force = PMExplicitForceInfo<dare::None>;
    using compressible = PMCompressibleInfo<false>;
};

struct PMNumericalInfoDefault {
    using tvd = PMTVDInfo<dare::MINMOD>;
    using viscous_stress = PMViscousStressInfo<PMDefaultStressTensor>;
    using momentum_iterations = PMMomentumIterationInfo<dare::FixedPoint>;
    using continuity_iterations = PMContinuityIterationInfo<dare::Newton>;
    using time_scheme_convective = PMTimeSchemeConvectiveInfo<dare::EULER_BACKWARD>;
    using momentum_normalizer = PMMomentumNormalizerInfo<dare::None>;
    using continuity_normalizer = PMContinuityNormalizerInfo<dare::None>;
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
class ProjectionMethod : public dare::InitializationTracker {
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
    using SC = typename GridType::ScalarType;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using FieldType = Field<GridType, SC, 1>;

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
    static const std::size_t dimension = GridType::Dimension;

    // algorithm and discretization properties
    using NumericalTypeInfo = detail::PMAssembledNumericalInfoWithDefaults<
                                        NumericalInfo,
                                        PMNumericalInfoDefault>;
    using TVDInfo = typename NumericalTypeInfo::tvd;
    using ViscousStressInfo = typename NumericalTypeInfo::viscous_stress;
    using MomentumIterationInfo = typename NumericalTypeInfo::momentum_iterations;
    using ContinuityIterationInfo = typename NumericalTypeInfo::continuity_iterations;
    using ConvectiveTimeSchemeInfo = typename NumericalTypeInfo::time_scheme_convective;
    using MomentumNormalizerInfo = typename NumericalTypeInfo::momentum_normalizer;
    using ContinuityNormalizerInfo = typename NumericalTypeInfo::continuity_normalizer;
    using TVDScheme = typename TVDInfo::type;
    using ViscousStressTreatment = typename ViscousStressInfo::type;
    using MomentumIterationType = typename MomentumIterationInfo::type;
    using ContinuityIterationType = typename ContinuityIterationInfo::type;
    using ConvectiveTimeSchemeType = typename ConvectiveTimeSchemeInfo::type;
    using MomentumNormalizerType = typename MomentumNormalizerInfo::type;
    using ContinuityNormalizerType = typename ContinuityNormalizerInfo::type;
    static const std::size_t num_tsteps_momentum
        = std::max(ConvectiveTimeSchemeType::NUM_TIMESTEPS + 1,
            static_cast<decltype(ConvectiveTimeSchemeType::NUM_TIMESTEPS)>(2));

    struct MomentumMembers{
        ExplicitForceMemberType beta_ex;
        MomentumNormalizerType normalizer;
    };
    struct ContinuityMembers {
        ExplicitForceMemberType beta_im;
        SC defect_max;
        ContinuityNormalizerType normalizer;
    };
    using MomentumType = dare::GenericEquation<GridType, BoundaryStrategyType, MomentumMembers>;
    using ContinuityType = PMContinuity<GridType, BoundaryStrategyType, ContinuityMembers>;

    ProjectionMethod()
        : ex_man(nullptr),
          dt(0.),
          max_iterations(100),
          status(0), status_finalized(rho_init | mu_init | epsilon_init | beta_im_init | beta_ex_init) {
        if constexpr (dare::is_none_v<PorosityVariableType>)
            status |= epsilon_init;
        if constexpr (dare::is_none_v<ImplicitForceVariableType>)
            status |= beta_im_init;
        if constexpr (dare::is_none_v<ExplicitForceVariableType>)
            status |= beta_ex_init;

        // method specific check for consistent compile time information
        free_compile_time_check(this);
    }

    template<TimeStepper T, typename... Args>
    void Initialize(GridType* grid, T* tstep, Args&&... bc_args) {
        ex_man = grid->GetExecutionManager();
        auto dt_obs_func = [&](const T& stepper, typename T::StateChange tag) {
            this->dt = stepper.GetTimeStepSize();
        };
        pimpl_dt_obs = dare::make_observer_handle<T>(dt_obs_func);
        dt = tstep->GetTimeStepSize();
        free_pm_initialize(this, grid, bc_args...);
        if constexpr(std::is_arithmetic_v<MomentumNormalizerType>) {
            for (auto& e : momentum)
                e->GetCustomMember()->normalizer = 1;
        }
        if constexpr (std::is_arithmetic_v<MomentumNormalizerType>) {
            continuity->GetCustomMember()->normalizer = 1;
        }
        this->dare::InitializationTracker::Initialize();
    }

    template<TimeStepper T, typename... Args>
    void Initialize(std::unique_ptr<GridType>& grid, T* tstep, Args&&... bc_args) {   // NOLINT
        Initialize(grid.get(), tstep, bc_args...);
    }

    // for access in the free functions
    template <typename... Args>
    void InitializeMomentum(std::size_t dim, Args&&... args) {
        momentum[dim] = std::make_unique<MomentumType>(args...);
    }

    template <typename... Args>
    void InitializeContinuity(Args&&... args) {
        continuity = std::make_unique<ContinuityType>(args...);
    }

    void SolveFlowField() {
        // here we could add switch between different approaches, but
        // for now there is only one
        if (!CheckStatus()) {
            ex_man->Terminate(__func__, "Projection method was not fully finalized!");
        }
        // build momentum and solve subsequently
        // a bit more verbose, but easier to debug
        // put it in a scope to limit variable lifetime
        {
            BuildMomentum(dare::ZERO);
            auto [success, iter] = SolveMomentum(0);
        }
        // add some output here
        if constexpr (dimension > 1) {
            BuildMomentum(dare::ONE);
            auto [success, iter] = SolveMomentum(1);
            // add some output here
        }
        if constexpr (dimension > 2) {
            BuildMomentum(dare::TWO);
            auto [success, iter] = SolveMomentum(2);
            // add some output here
        }

        // enforce continuity
        // first iteration we add the implicit force term
        // or in the case of no additional term this is just
        // the very first iteration
        int iteration = 0;
        for (; iteration < max_iterations; iteration++) {
            BuildContinuity(iteration);
            auto [success, iter] = SolveContinuity(iteration);
            UpdatePressure(iteration);
            UpdateVelocity(iteration);
            if (ContinuityConvergence(iteration))
                break;
        }

        // check if iterations == max_iterations
        if (iteration == (max_iterations - 1)) {
            // add warning for unconverged solution
        }

        ERROR << "Implementation not finished" << ERROR_CLOSE;
    }

    constexpr bool IsCompressible() const { return compressible; }
    constexpr std::size_t GetDimension() const { return dimension; }

    std::unique_ptr<MomentumType>& GetMomentum(std::size_t dim) { return momentum[dim]; }
    const std::unique_ptr<MomentumType>& GetMomentum(std::size_t dim) const { return momentum[dim]; }

    std::unique_ptr<ContinuityType>& GetContinuity() { return continuity; }
    const ContinuityType& GetContinuity() const { return *continuity; }

    FieldType* GetPressure() { return continuity->GetPressure(); }
    const FieldType& GetPressure() const { return continuity->GetPressure(); }

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

    DensityVariableType GetDensity() const {
        return rho;
    }

    ViscosityVariableType GetViscosity() const {
        return rho;
    }

    PorosityVariableType GetPorosity() const {
        return epsilon;
    }

    void AddImplicitForce(ImplicitForceVariableType f) {
        // add to continuity
        if (!IsInitialized()) {
            ex_man->Terminate(__func__, "Cannot add force terms prior to initialization");
        }
        if constexpr (!dare::is_none_v<ImplicitForceVariableType>) {
            continuity->GetCustomMember()->beta_im.emplace(f);
        }
        status |= beta_im_init;
    }

    void AddExplicitForce(ExplicitForceVariableType f, std::size_t dim) {
        // add to continuity
        if (!IsInitialized()) {
            ex_man->Terminate(__func__, "Cannot add force terms prior to initialization");
        }
        if constexpr (!dare::is_none_v<ExplicitForceVariableType>) {
            if (dim >= dimension) {
                ex_man->Terminate(__func__, "Invalid dimension choses for the force");
            }
            momentum[dim]->GetCustomMember()->beta_ex.emplace(f);

            // check if all explicit force members were set
            char all_init{beta_ex_init};
            for (auto& m : momentum)
                all_init &= !m->GetCustomMember()->beta_ex.empty()? beta_ex_init : 0;
            status |= (beta_ex_init & all_init);
        }
    }

    bool CheckStatus() const {
        return (status == status_finalized) && this->IsInitialized();
    }

    void SetTimeStepSize(SC _dt) {
        // This should be a callback
        dt = _dt;
    }

    SC GetTimeStepSize() const {
        return dt;
    }

    dare::ExecutionManager* GetExecutionManager() const {
        return ex_man;
    }

private:
    template <dare::NaturalNumber Direction>
    void BuildMomentum(Direction dir) {
        free_pm_build_momentum(this, dir);
    }

    template <dare::NaturalNumber Direction>
    std::pair<bool, int> SolveMomentum(Direction dir) {
        return free_pm_solve_momentum(this, dir);
    }

    void BuildContinuity(int iteration) {
        free_pm_build_continuity(this, iteration);
    }

    void SolveContinuity(int iteration) {
        free_pm_solve_continuity(this, iteration);
    }

    void UpdatePressure(int iteration) {
        free_pm_update_pressure(this, iteration);
    }

    void UpdateVelocity(int iteration) {
        free_pm_update_velocity(this, iteration);
    }

    bool ContinuityConvergence(int iteration) {
        return free_pm_continuity_convergence(this, iteration);
    }

    dare::ExecutionManager* ex_man;
    DensityVariableType rho;
    ViscosityVariableType mu;
    PorosityVariableType epsilon;

    std::unique_ptr<ContinuityType> continuity;
    std::array<std::unique_ptr<MomentumType>, dimension> momentum;

    SC dt;  // for now this is temporary, work with observer here!
    int max_iterations;
    char status;
    char status_finalized;
    UniqueObserverHandle pimpl_dt_obs;
};

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_
