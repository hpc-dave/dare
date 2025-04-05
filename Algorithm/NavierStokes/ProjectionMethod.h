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
#include <set>

#include "PM_Information.h"
#include "PM_details.h"
#include "PM_Continuity.h"
#include "Algorithm/AlgorithmTraits.h"
#include "Data/Field.h"
#include "Equations/FluxLimiter.h"
#include "Utilities/Errors.h"
#include "Utilities/Observer.h"
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
    using density_derivative = PMDensityDerivativeInfo<dare::None>;
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
        beta_ex_init = 0b0010000,
        density_derivative_init = 0b0100000
    };

    /*!
     * @brief an enum class for potential observers
     */
    enum class StateChange {
        VelocityUpdated,
        PressureUpdated
    };

    // general types based on the grid
    using GridType = Grid;
    using BoundaryStrategyType = BoundaryStrategy;
    using SelfType = ProjectionMethod<Grid, BoundaryStrategy, PropertyInfo, NumericalInfo>;
    using SC = typename GridType::ScalarType;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using FieldType = dare::Field<GridType, SC, 1>;
    using GridVectorType = dare::GridVector<GridType, SC, 1>;
    using ObserverType = dare::Observer<SelfType, StateChange>;

    // properties determined from the PropertyInfo type
    using PropertyTypeInfo = detail::PMAssembledPropertyInfoWithDefaults<
                                        PropertyInfo,
                                        PMPropertyInfoDefault<Grid>>;
    using DensityInfo = typename PropertyTypeInfo::density;
    using ViscosityInfo = typename PropertyTypeInfo::viscosity;
    using PorosityInfo = typename PropertyTypeInfo::porosity;
    using ImplicitForceInfo = typename PropertyTypeInfo::implicit_force;
    using ExplicitForceInfo = typename PropertyTypeInfo::explicit_force;
    // using CompressibilityInfo = typename PropertyTypeInfo::compressible;
    using DensityDerivativeInfo = typename PropertyTypeInfo::density_derivative;
    using DensityVariableType = detail::determine_density_variable_type_t<DensityInfo>;
    using ViscosityVariableType = detail::determine_viscosity_variable_type_t<ViscosityInfo>;
    using PorosityVariableType = detail::determine_porosity_variable_type_t<PorosityInfo>;
    using ImplicitForceVariableType = detail::determine_implicit_force_variable_type_t<ImplicitForceInfo>;
    using ExplicitForceVariableType = detail::determine_explicit_force_variable_type_t<ExplicitForceInfo>;
    using ImplicitForceMemberType = detail::determine_implicit_force_member_variable_type_t<ImplicitForceInfo>;
    using ExplicitForceMemberType = detail::determine_explicit_force_member_variable_type_t<ExplicitForceInfo>;
    using DensityDerivativeMemberType = detail::determine_density_derivative_member_variable_type_t<DensityDerivativeInfo>; // NOLINT
    // consider removing this boolean and check simply for equation of state
    static const bool compressible = !dare::is_none_v<DensityDerivativeMemberType>;
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
        ImplicitForceMemberType beta_im;
        SC defect_max;
        ContinuityNormalizerType normalizer;
    };
    using MomentumType = dare::GenericEquation<GridType, BoundaryStrategyType, MomentumMembers>;
    using ContinuityType = PMContinuity<GridType, BoundaryStrategyType, ContinuityMembers>;

    ProjectionMethod()
        : ex_man{nullptr},
          dt{0.},
          continuity_tolerance{1e-14},
          max_iterations{100},
          status{0},
          status_finalized{rho_init | mu_init | epsilon_init | beta_im_init | beta_ex_init | density_derivative_init} {
        if constexpr (dare::is_none_v<PorosityVariableType>)
            status |= epsilon_init;
        if constexpr (dare::is_none_v<ImplicitForceVariableType>)
            status |= beta_im_init;
        if constexpr (dare::is_none_v<ExplicitForceVariableType>)
            status |= beta_ex_init;
        if constexpr(dare::is_none_v<DensityDerivativeMemberType>)
            status |= density_derivative_init;

        // method specific check for consistent compile time information
        free_compile_time_check(this);
    }

    template<TimeStepper T, typename... Args>
    void Initialize(GridType* grid, T* tstep, BoundaryStrategyType bc_continuity, Args... bc_momentum) {
        ex_man = grid->GetExecutionManager();
        auto dt_obs_func = [&](const T& stepper, typename T::StateChange tag) {
            this->dt = stepper.GetTimeStepSize();
        };
        pimpl_dt_obs = dare::make_observer_handle<T>(dt_obs_func);
        dt = tstep->GetTimeStepSize();
        free_pm_initialize(this, grid, bc_continuity, bc_momentum...);
        if constexpr(std::is_arithmetic_v<MomentumNormalizerType>) {
            for (auto& e : momentum)
                e->GetCustomMember()->normalizer = 1;
        }
        if constexpr (std::is_arithmetic_v<MomentumNormalizerType>) {
            continuity->GetCustomMember()->normalizer = 1;
        }
        this->dare::InitializationTracker::Initialize();
    }

    template <TimeStepper T, typename... Args>
    void Initialize(std::unique_ptr<GridType>& grid, T* tstep, BoundaryStrategyType bc_continuity, Args... bc_momentum) {  // NOLINT
        Initialize(grid.get(), tstep, bc_continuity, bc_momentum...);
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
            auto [success, iter] = SolveMomentum(dare::ZERO);
            ERROR << "this is a placeholder - ignore for now: " << iter << " "
                << (success? "success": "fail") << ERROR_CLOSE;
        }
        // add some output here
        if constexpr (dimension > 1) {
            BuildMomentum(dare::ONE);
            auto [success, iter] = SolveMomentum(dare::ONE);
            // add some output here
            ERROR << "this is a placeholder - ignore for now: " << iter << " "
                << (success ? "success" : "fail") << ERROR_CLOSE;
        }
        if constexpr (dimension > 2) {
            BuildMomentum(dare::TWO);
            auto [success, iter] = SolveMomentum(dare::TWO);
            // add some output here
            ERROR << "this is a placeholder - ignore for now: " << iter << " "
                << (success ? "success" : "fail") << ERROR_CLOSE;
        }

        // enforce continuity
        // first iteration we add the implicit force term
        // or in the case of no additional term this is just
        // the very first iteration
        int iteration = 0;
        ComputeDefect();
        for (; iteration < max_iterations; iteration++) {
            if constexpr (compressible)
                *rho_prev = rho->GetGridVector(0);

            BuildContinuity(iteration);
            auto [success, iter] = SolveContinuity(iteration);

            if (!success) {
                ex_man->Print(dare::Verbosity::Low) << "Continuity system failed to converge after "
                                             << iter << " matrix-solver iterations";
            }

            UpdatePressure();
            Notify(StateChange::PressureUpdated);

            UpdateVelocity(iteration);
            Notify(StateChange::VelocityUpdated);

            ComputeDefect();

            // print update to terminal
            if (ContinuityConvergence())
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
        if constexpr(compressible) {
            static_assert(dare::is_field_v<std::remove_cv_t<std::remove_pointer_t<DensityVariableType>>>,
            "In the compressible case, the density needs to be a field!");
            rho_prev = std::make_unique<GridVectorType>("rho_prev", rho->GetGridRepresentation());
            *rho_prev = rho->GetDataVector();
        }
    }
    void SetViscosity(ViscosityVariableType v) {
        mu = v;
        status |= mu_init;
    }
    void SetPorosity(PorosityVariableType p) {
        epsilon = p;
        status |= epsilon_init;
    }

    void SetDensityDerivative(DensityDerivativeMemberType dd) {
        density_derivative = std::move(dd);
        status |= density_derivative_init;
    }

    const DensityDerivativeMemberType& GetDensityDerivative() const {
        return density_derivative;
    }

    SC GetDensityDerivative(Index ind) const {
        return density_derivative(ind);
    }

    DensityVariableType GetDensity() const {
        return rho;
    }

    const GridVectorType* GetDensityPreviousIteration() const {
        if constexpr(!compressible) {
            static_assert(dare::always_false<decltype(this)>, "In the incompressible case this should not be accessed");
        }
#ifndef DARE_NDEBUG
        if (!rho_prev) {
            this->ex_man->Terminate(__func__, "The array for the previous iteration was not allocated!");
        }
#endif
        return rho_prev.get();
    }

    SC GetDensity(Index ind, std::size_t time_level = 0) const {
        using DType = std::remove_cv_t<std::remove_pointer_t<DensityVariableType>>;
        if constexpr (dare::is_field_v<DType>) {
            return rho->GetDataVector(time_level).At(ind, 0);
        } else if constexpr(std::is_arithmetic_v<DType>) {
            return rho;
        } else {
            static_assert(dare::always_false<DType>, "Cannot handle this type");
        }
    }

    ViscosityVariableType GetViscosity() const {
        return rho;
    }

    SC GetViscosity(Index ind, std::size_t time_level = 0) const {
        using VType = std::remove_cv_t<std::remove_pointer_t<ViscosityVariableType>>;
        if constexpr (dare::is_field_v<VType>) {
            return mu->GetDataVector(time_level).At(ind, 0);
        } else if constexpr (std::is_arithmetic_v<VType>) {
            return mu;
        } else {
            static_assert(dare::always_false<VType>, "Cannot handle this type");
        }
    }

    PorosityVariableType GetPorosity() const {
        return epsilon;
    }

    SC GetPorosity(Index ind, std::size_t time_level = 0) const {
        using PType = std::remove_cv_t<std::remove_pointer_t<PorosityVariableType>>;
        if constexpr (dare::is_field_v<PType>) {
            return epsilon->GetDataVector(time_level).At(ind, 0);
        } else if constexpr (std::is_arithmetic_v<PType>) {
            return epsilon;
        } else if constexpr (dare::is_none_v<PType>) {
            return 1.;
        } else {
            static_assert(dare::always_false<PType>, "Cannot handle this type");
        }
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

    SC GetTimeStepSize() const {
        return dt;
    }

    dare::ExecutionManager* GetExecutionManager() const {
        return ex_man;
    }

    bool Attach(ObserverType* o) {
        auto [pos, success] = observers.emplace(o);
        return success;
    }

    bool Detach(ObserverType* o) {
        return (observers.erase(o) > 0U);
    }

    void Notify(StateChange property) {
        if constexpr (compressible) {
            if (property == StateChange::PressureUpdated && (observers.size() == 0)) {
                ERROR << "Notifying others of pressure update, however no observers were attached!"
                << " In the compressible case the means that the density is not adapted!" << ERROR_CLOSE;
            }
        }
        for (auto iter = observers.begin(); iter != observers.end();) {
            auto const pos = iter++;
            (*pos)->Update(*this, property);
        }
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

    std::pair<bool, int> SolveContinuity(int iteration) {
        return free_pm_solve_continuity(this, iteration);
    }

    void UpdatePressure() {
        free_pm_update_pressure(this);
    }

    void UpdateVelocity(int iteration) {
        free_pm_update_velocity(this, iteration);
    }

    void ComputeDefect() {
        free_pm_compute_defect(this);
    }

    SC DetermineMaxContinuityDefect() {
        return free_pm_determine_max_continuity_defect(this);
    }

    bool ContinuityConvergence() {
        if constexpr (uses_newton_iterations_v<ContinuityIterationType>) {
            return max_continuity_defect < continuity_tolerance;
        } else {
            static_assert(dare::always_false<decltype(this)>, "At the moment, only newton iterations are allowed for the continuity");  // NOLINT
        }
    }

    dare::ExecutionManager* ex_man;
    DensityVariableType rho;
    ViscosityVariableType mu;
    PorosityVariableType epsilon;
    DensityDerivativeMemberType density_derivative;
    std::unique_ptr<GridVectorType> rho_prev;

    std::unique_ptr<ContinuityType> continuity;
    std::array<std::unique_ptr<MomentumType>, dimension> momentum;

    SC dt;  // for now this is temporary, work with observer here!
    SC continuity_tolerance;
    SC max_continuity_defect;
    int max_iterations;
    char status;
    char status_finalized;
    UniqueObserverHandle pimpl_dt_obs;
    std::set<ObserverType*> observers;
};

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_
