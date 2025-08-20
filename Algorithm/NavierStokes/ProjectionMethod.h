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

#include <algorithm>
#include <array>
#include <bitset>
#include <concepts>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Algorithm/AlgorithmTraits.h"
#include "Data/Field.h"
#include "Equations/FluxLimiter.h"
#include "Equations/GenericEquation.h"
#include "Equations/TimeDiscretizationSchemes.h"
#include "IO/TerminalOutput.h"
#include "PMTerminalOutputDefault.h"
#include "PM_Continuity.h"
#include "PM_Information.h"
#include "PM_details.h"
#include "ProjectionMethod_freefunc.h"
#include "Utilities/Errors.h"
#include "Utilities/InitializationTracker.h"
#include "Utilities/Observer.h"
#include "Utilities/Timer.h"

namespace dare {

/*!
 * @brief a default boundary strategy utilizing type erasure for variability
 */
template<typename Grid, typename SC>
class PMDefaultBoundaryType {
public:
    using GridType = Grid;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using FieldType = dare::Field<GridType, SC, 1>;
    using MBTypeLO = dare::MatrixBlock<GridType, LO, SC, 1>;
    using MBTypeGO = dare::MatrixBlock<GridType, GO, SC, 1>;

    template<typename T>
    explicit PMDefaultBoundaryType(const T& bc)
        : pimpl(std::make_unique<Model<T>>(bc)) {}

    ~PMDefaultBoundaryType() = default;

    PMDefaultBoundaryType(const PMDefaultBoundaryType& other)
        : pimpl(other.pimpl->clone()) {}

    PMDefaultBoundaryType& operator=(PMDefaultBoundaryType other) {
        std::swap(pimpl, other.pimpl);
        return *this;
    }

    PMDefaultBoundaryType(PMDefaultBoundaryType&& other) = default;

    PMDefaultBoundaryType& operator=(PMDefaultBoundaryType&& other) = default;

    /*!
     * @brief apply the boundary condition to the local matrix block or field
     */
    void Apply(MBTypeLO* mb) const { pimpl->Apply(mb); }

    /*!
     * @brief apply the boundary condition to the global matrix block or field
     */
    void Apply(MBTypeGO* mb) const { pimpl->Apply(mb); }

    /*!
     * @brief adapt the field according to the boundary conditions
     */
    void Apply(FieldType* f) const { pimpl->Apply(f); }

    /*!
     * @brief apply the boundary condition to the local matrix block or field
     */
    void operator()(MBTypeLO* mb) const {
        pimpl->Apply(mb);
    }

    /*!
     * @brief apply the boundary condition to the global matrix block or field
     */
    void operator()(MBTypeGO* mb) const {
        pimpl->Apply(mb);
    }

    /*!
     * @brief adapt the field according to the boundary conditions
     */
    void operator()(FieldType* f) const {
        pimpl->Apply(f);
    }

private:
    /*!
     * @brief abstract virtual class for the concept of the type erasure defining the type signature
     */
    class Concept {
    public:
        virtual ~Concept() = default;
        virtual void Apply(MBTypeLO*) const = 0;
        virtual void Apply(MBTypeGO*) const = 0;
        virtual void Apply(FieldType*) const = 0;

        /*!
         * @brief prototype signature clone
         * @return unique base pointer of copied object
         */
        virtual std::unique_ptr<Concept> clone() const = 0;
    };

    /*!
     * @brief wrapper for the type erasure containing the original object
     * @tparam TBoundaryCondition type of the wrapped object
     */
    template<typename TBoundaryCondition>
    class Model : public Concept {
    public:
        explicit Model(const TBoundaryCondition& bc) : boundary_condition(bc) {}

        void Apply(MBTypeLO* mb) const override {
            boundary_condition.Apply(mb);
        }
        void Apply(MBTypeGO* mb) const override {
            boundary_condition.Apply(mb);
        }
        void Apply(FieldType* f) const override {
            boundary_condition.Apply(f);
        }

        /*!
         * @brief prototype clone instantiation
         * @return unique pointer of the copied object
         */
        std::unique_ptr<Concept> clone() const override {
            return std::make_unique<Model<TBoundaryCondition>>(*this);
        }

    private:
        TBoundaryCondition boundary_condition;  //!< original object
    };

    std::unique_ptr<Concept> pimpl;     //!< implementation of the concept
};

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
template <typename Grid, typename BoundaryInfo, typename PropertyInfo, typename NumericalInfo>
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
    // using BoundaryStrategyType = BoundaryStrategy;
    using SelfType = ProjectionMethod<Grid, BoundaryInfo, PropertyInfo, NumericalInfo>;
    using SC = typename GridType::ScalarType;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using FieldType = dare::Field<GridType, SC, 1>;
    using GridVectorType = dare::GridVector<GridType, SC, 1>;
    using ObserverType = dare::Observer<SelfType, StateChange>;
    using TerminalOutput = dare::PMTerminalOutputDefault;

    // boundary information for the different equations
    using BoundaryInfoType = BoundaryInfo;

    // properties determined from the PropertyInfo type
    using PropertyTypeInfo = detail::PMAssembledPropertyInfoWithDefaults<
                                        PropertyInfo,
                                        PMPropertyInfoDefault<Grid>>;
    using DensityInfo = typename PropertyTypeInfo::density;
    using ViscosityInfo = typename PropertyTypeInfo::viscosity;
    using PorosityInfo = typename PropertyTypeInfo::porosity;
    using ImplicitForceInfo = typename PropertyTypeInfo::implicit_force;
    using ExplicitForceInfo = typename PropertyTypeInfo::explicit_force;
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

    // boundary conditions (for now only the boundary info type, adapt at later point to be more flexible)
    using BoundaryStrategyContinuity = BoundaryInfoType;
    using BoundaryStrategyMomentum = BoundaryInfoType;

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
    using MomentumType = dare::GenericEquation<GridType, BoundaryStrategyMomentum, MomentumMembers>;
    using ContinuityType = PMContinuity<GridType, BoundaryStrategyContinuity, ContinuityMembers>;

    /*!
     * @brief default constructor
     */
    ProjectionMethod()
        : ex_man{nullptr},
          dt{0.},
          time{0.},
          tstep{0},
          continuity_tolerance{1e-14},
          max_iterations{100},
          status{0},
          status_finalized{rho_init | mu_init | epsilon_init | beta_im_init | beta_ex_init | density_derivative_init},
          terminal_output(dimension) {
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

    /*!
     * @brief allocates memory and instantiates the separate equations
     * @tparam ...Args input arguments for the momentum equation
     * @tparam BCContinuity boundary condition input type for the continuity
     * @tparam T a time step object
     * @param grid pointer to the grid
     * @param tstep pointer to the time stepper
     * @param bc_continuity boundary arguments for continuity
     * @param ...bc_momentum boundary arguments for the momentum
     */
    template<TimeStepper T, typename BCContinuity, typename... Args>
    void Initialize(GridType* grid, T* tstep, BCContinuity bc_continuity, Args... bc_momentum) {
        ex_man = grid->GetExecutionManager();
        auto dt_obs_func = [&](const T& stepper, typename T::StateChange tag) {
            using State = typename T::StateChange;
            switch (tag) {
                case State::AdvanceTimeStep:
                    this->time = stepper.GetTime();
                    this->tstep = stepper.GetTimeStepCounter();
                    break;
                case State::UpdateTimeStepSize:
                    this->dt = stepper.GetTimeStepSize();
                    break;
                }
        };
        pimpl_dt_obs = dare::make_observer_handle(tstep, dt_obs_func);
        dt = tstep->GetTimeStepSize();
        time = tstep->GetTime();
        free_pm_initialize(this, grid, bc_continuity, bc_momentum...);

        // Provide default solver settings
        free_pm_solver_settings_default(this);

        if constexpr (std::is_arithmetic_v<MomentumNormalizerType>) {
            for (auto& e : momentum)
                e->GetCustomMember()->normalizer = 1;
        }
        if constexpr (std::is_arithmetic_v<MomentumNormalizerType>) {
            continuity->GetCustomMember()->normalizer = 1;
        }
        this->dare::InitializationTracker::Initialize();
    }

        /*!
         * @brief overload for the case the the grid is provided as unique pointer
         * @tparam ...Args input arguments for the momentum equation
         * @tparam BCContinuity boundary condition input type for the continuity
         * @tparam T a time step object
         * @param grid pointer to the grid
         * @param tstep pointer to the time stepper
         * @param bc_continuity boundary arguments for continuity
         * @param ...bc_momentum boundary arguments for the momentum
         */
        template <TimeStepper T, typename BCContinuity, typename... Args>
        void Initialize(std::unique_ptr<GridType> & grid, T * tstep, BCContinuity bc_continuity, Args... bc_momentum) {  // NOLINT
            Initialize(grid.get(), tstep, bc_continuity, bc_momentum...);
        }

    /*!
     * @brief initializing a dediacted momentum equation
     * @tparam ...Args input arguments for the momentum equations
     * @param dim the dimension for which the momentum should be initialized
     * @param ...args actual arguments for the momentum
     * @warning only for access in free functions!
     */
    template <typename... Args>
    void InitializeMomentum(std::size_t dim, Args&&... args) {
        momentum[dim] = std::make_unique<MomentumType>(args...);
    }

    /*!
     * @brief initializes the continuity
     * @tparam ...Args arguments for the continuity
     * @param ...args arguments for the continuity
     * @warning only for access in free functions!
     */
    template <typename... Args>
    void InitializeContinuity(Args&&... args) {
        continuity = std::make_unique<ContinuityType>(args...);
    }

    /*!
     * @brief advance the flow field by a 2 step projection, based on Chorin's projection method
     */
    void SolveFlowField() {
        // here we could add switch between different approaches, but
        // for now there is only one
        if (!CheckStatus()) {
            ex_man->Terminate(__func__, "Projection method was not fully finalized!");
        }

        terminal_output.PrintHeader(*this);
        StartProfiling("FlowStep");
        // build momentum and solve subsequently
        // a bit more verbose, but easier to debug
        // put it in a scope to limit variable lifetime
        {
            BuildMomentum(dare::ZERO);
            auto [success, iter] = SolveMomentum(dare::ZERO);
            terminal_output.PrintMomentum(0, iter, success, this);
        }
        // add some output here
        if constexpr (dimension > 1) {
            BuildMomentum(dare::ONE);
            auto [success, iter] = SolveMomentum(dare::ONE);
            // add some output here
            terminal_output.PrintMomentum(1, iter, success, *this);
        }
        if constexpr (dimension > 2) {
            BuildMomentum(dare::TWO);
            auto [success, iter] = SolveMomentum(dare::TWO);
            terminal_output.PrintMomentum(2, iter, success, *this);
        }

        // enforce continuity
        // first iteration we add the implicit force term
        // or in the case of no additional term this is just
        // the very first iteration
        StartProfiling("Enforce_Continuity");
        int iteration = 0;
        ComputeDefect();
        terminal_output.PrintInitialDefect(*this);
        for (; iteration < max_iterations; iteration++) {
            if constexpr (compressible)
                *rho_prev = rho->GetGridVector(0);

            BuildContinuity(iteration);
            auto [success, iter] = SolveContinuity(iteration);

            if (!success) {
                dare::Print(dare::Verbosity::Low) << "Continuity system failed to converge after "
                                                  << iter << " matrix-solver iterations";
            }

            UpdatePressure();
            Notify(StateChange::PressureUpdated);

            UpdateVelocity(iteration);
            Notify(StateChange::VelocityUpdated);

            ComputeDefect();

            // print update to terminal
            terminal_output.PrintContinuity(iteration, iter, success, *this);
            if (ContinuityConvergence() || !success)
                break;
        }
        StopProfiling("Enforce_Continuity");
        StopProfiling("FlowStep");

        // check if iterations == max_iteration
        if (iteration == max_iterations) {
            Print(dare::Verbosity::Low) << "The continuity could not be conserved within "
                << std::to_string(max_iterations) << " iterations!" << std::endl;
            iteration--;
        }
        PrintProfiling(iteration);
    }

    /*!
     * @brief convenient constexpr access to compressibility
     * @return true if compressible
     */
    constexpr bool IsCompressible() const { return compressible; }

    /*!
     * @brief convenient constexpr access to the dimension of the system
     */
    constexpr std::size_t GetDimension() const { return dimension; }

    /*!
     * @brief provides momentum equation
     * @param dim direction of the momentum
     * @return reference to the unique_ptr holding the momentum
     */
    std::unique_ptr<MomentumType>& GetMomentum(std::size_t dim) { return momentum[dim]; }
    const std::unique_ptr<MomentumType>& GetMomentum(std::size_t dim) const { return momentum[dim]; }

    /*!
     * @brief provides continuity equation
     * @return reference to the unique_ptr holding the continuity
     */
    std::unique_ptr<ContinuityType>& GetContinuity() { return continuity; }
    const ContinuityType& GetContinuity() const { return *continuity; }

    /*!
     * @brief convenient access to the pressure field
     */
    FieldType* GetPressure() { return continuity->GetPressure(); }
    const FieldType& GetPressure() const { return continuity->GetPressure(); }

    /*!
     * @brief setter for the density
     * @param d density variable
     * Internally, the status is updated to track the density
     */
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

    /*!
     * @brief setter for the viscosity
     * @param v viscosity variable
     * Internally, the status is updated to track the viscosity
     */
    void SetViscosity(ViscosityVariableType v) {
        mu = v;
        status |= mu_init;
    }

    /*!
     * @brief setter for the porosity
     * @param p porosity variable
     * Internally, the status is updated to track the porosity
     */
    void SetPorosity(PorosityVariableType p) {
        epsilon = p;
        status |= epsilon_init;
    }

    /*!
     * @brief setter for the density derivative
     * @param dd density derivative variable
     * Internally, the status is updated to track the density derivative
     */
    void SetDensityDerivative(DensityDerivativeMemberType dd) {
        density_derivative = std::move(dd);
        status |= density_derivative_init;
    }

    /*!
     * @brief provides access to the underlying density derivative
     */
    const DensityDerivativeMemberType& GetDensityDerivative() const {
        return density_derivative;
    }

    /*!
     * @brief provides the density derivative at a certain index
     * @param ind local index
     * @return value at the cell center
     */
    SC GetDensityDerivative(Index ind) const {
        return density_derivative(ind);
    }

    /*!
     * @brief provides the underlying density variable
     */
    DensityVariableType GetDensity() const {
        return rho;
    }

    /*!
     * @brief provides the field with porosity at the previous iteration
     * @return pointer to the field
     * This function can only be used in the compressible case and won't compile
     * for incompressible flow!
     */
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

    /*!
     * @brief returns the density at a certain index
     * @param ind local index at which to provide the density
     * @param time_level time step at which to access
     * @return value at the cell center
     */
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

    /*!
     * @brief provides the underlying viscosity variable
     */
    ViscosityVariableType GetViscosity() const {
        return mu;
    }

    /*!
     * @brief returns the viscosity at a certain index
     * @param ind local index at which to provide the viscosity
     * @param time_level time step at which to access
     * @return value at the cell center
     */
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

    /*!
     * @brief provides the underlying porosity variable
     */
    PorosityVariableType GetPorosity() const {
        return epsilon;
    }

    /*!
     * @brief returns the porosity at a certain index
     * @param ind local index at which to provide the porosity
     * @param time_level time step at which to access
     * @return value at the cell center
     */
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

    /*!
     * @brief Adds an imiplicit force to the flow
     * @param f implicit force variable type
     * The implicit force is applied to the continuity equation and thus has to be
     * provided for the scalar grid
     */
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

    /*!
     * @brief Adds an explicit force to the flow
     * @param f explicit force variable type
     * @param dim staggered grid direction at which the explicit force is applied
     * The explict force is applied to the momentum equation and thus has to be
     * supplied in the staggered configuration
     */
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

    /*!
     * @brief Checks if the object is ready for flow step computation
     * @return true if prepared
     */
    bool CheckStatus() const {
        return (status == status_finalized) && this->IsInitialized();
    }

    /*!
     * @brief current time step size
     * @return 
     */
    SC GetTimeStepSize() const {
        return dt;
    }

    /*!
     * @brief returns the current simulation time
     */
    SC GetTime() const {
        return time;
    }

    /*!
     * @brief returns the current time step number
     * @return time step number
     */
    TimeStepCounter GetTimeStepCounter() const {
        return tstep;
    }

    /*!
     * @brief Sets the maximum number of loop iterations
     * @param max_loops max loops
     */
    void SetMaxLoopIterations(int max_loops) {
        if (max_loops < 0) {
            ERROR << "maximum continuity iterations may not be negative!" << ERROR_CLOSE;
            return;
        }
        max_iterations = max_loops;
    }

    /*!
     * @brief maximum number of continuity loops
     */
    int GetMaxLoopIterations() const {
        return max_iterations;
    }

    /*!
     * @brief provides the last determined maximum continuity defect
     * @return max continuity defect
     */
    SC GetMaxContinuityDefect() const {
        return max_continuity_defect;
    }

    /*!
     * @brief provides reference to execution manager
     * @return address of execution manager
     */
    dare::ExecutionManager* GetExecutionManager() const {
        return ex_man;
    }

    /*!
     * @brief attached an observer
     * @param o address of the observer
     * @return true, if the observer could be attached
     */
    bool Attach(ObserverType* o) {
        auto [pos, success] = observers.emplace(o);
        return success;
    }

    /*!
     * @brief detaches a previous attached observer
     * @param o address of the observer
     * @return true, if observer was found and deattached
     */
    bool Detach(ObserverType* o) {
        return (observers.erase(o) > 0U);
    }

    /*!
     * @brief notifies the attached observers of a state change
     * @param property state change
     */
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

    /*!
     * @brief Convenient overload for copying all internal fields to the old timestep
     */
    void CopyToOld() {
        for (std::size_t i{0}; i < dimension; i++)
            GetMomentum(i)->GetField()->CopyDataVectorsToOldTimeStep();
        GetContinuity()->GetField()->CopyDataVectorsToOldTimeStep();
    }

    /*!
     * @brief provides a timer and starts profiling
     * @param t instance of the timer
     */
    void SetTimer(Timer t) {
        timer = std::make_unique<Timer>(std::move(t));
    }

private:
    /*!
     * @brief builds the matrix system for the momentum equation in specific direction
     * @tparam Direction Natural number type
     * @param dir direction of the momentum equation
     */
    template <dare::NaturalNumber Direction>
    void BuildMomentum(Direction dir) {
        std::string id = std::string {"Build_Momentum_"} + std::to_string(dir);
        StartProfiling(id);
        free_pm_build_momentum(this, dir);
        StopProfiling(id);
    }

    /*!
     * @brief solves a specific momentum equation
     * @tparam Direction Natural number type
     * @param dir direction of the momentum equation
     * @return pair of <bool, int>, indicating convergence success and number of solver loops
     */
    template <dare::NaturalNumber Direction>
    std::pair<bool, int> SolveMomentum(Direction dir) {
        std::string id = std::string {"Solve_Momentum_"} + std::to_string(dir);
        StartProfiling(id);
        auto ret = free_pm_solve_momentum(this, dir);
        StopProfiling(id);
        return ret;
    }

    /*!
     * @brief builds the continuity matrix system
     * @param iteration continuity loop iteration
     */
    void BuildContinuity(int iteration) {
        std::string id = std::string {"Build_Continuity_"} + std::to_string(iteration);
        StartProfiling(id);
        free_pm_build_continuity(this, iteration);
        StopProfiling(id);
    }

    /*!
     * @brief calls the solver
     * @param iteration continuity loop iteration
     * @return pair of <bool, int>, indicating convergence success and number of solver loops
     */
    std::pair<bool, int> SolveContinuity(int iteration) {
        std::string id = std::string {"Solve_Continuity_"} + std::to_string(iteration);
        StartProfiling(id);
        auto ret = free_pm_solve_continuity(this, iteration);
        StopProfiling(id);
        return ret;
    }

    /*!
     * @brief updates the pressure with the dP field
     */
    void UpdatePressure() {
        free_pm_update_pressure(this);
    }

    /*!
     * @brief updates the velocity with the dP field
     * @param iteration continuity loop iteration
     * \note implementation is provided as an overload of the
     * function free_pm_update_velocity
     */
    void UpdateVelocity(int iteration) {
        free_pm_update_velocity(this, iteration);
    }

    /*!
     * @brief computes the defect according to the method specific function
     * Also determines the maximum defect
     */
    void ComputeDefect() {
        free_pm_compute_defect(this);
        max_continuity_defect = DetermineMaxContinuityDefect();
    }

    /*!
     * @brief determines maximum defect on the defect field
     * @return maximum defect
     * \note implementation can be customized depending
     * on the template parameters of the class by providing
     * a dedicated overload of the underlying function
     * free_pm_determine_max_continuity_defect
     */
    SC DetermineMaxContinuityDefect() {
        return free_pm_determine_max_continuity_defect(this);
    }

    /*!
     * @brief checks if continuity has been achieved within tolerance
     * @return true, if converged
     */
    bool ContinuityConvergence() {
        if constexpr (uses_newton_iterations_v<ContinuityIterationType>) {
            return max_continuity_defect < continuity_tolerance;
        } else {
            static_assert(dare::always_false<decltype(this)>, "At the moment, only newton iterations are allowed for the continuity");  // NOLINT
        }
    }

    /*!
     * @brief starts profiling, if timer was allocated
     * @param id profile step id
     */
    void StartProfiling(std::string id) {
        if (timer)
            timer->Tic(id);
    }

    /*!
     * @brief stops profiling, if timer was allocated
     * @param id profile step id
     * @return elapsed time
     */
    Timer::ValueType StopProfiling(std::string id) {
        if (timer)
            return timer->Toc(id);
        return 0.;
    }

    /*!
     * @brief pretty printing of the profiling data
     * @param c_loops numer of loops for continuity
     * 
     * \note if no timer was allocated, this function returns without doing anything
     */
    void PrintProfiling(int c_loops) const {
        if (!timer)
            return;

        char mom_id[] = { 'X', 'Y', 'Z' };

        const int prec{3};  //!< precision of the output
        Timer::ValueType t_flowstep = timer->GetElapsedTime("FlowStep");
        std::array<Timer::ValueType, dimension> t_build_mom, t_solve_mom;
        for (std::size_t d{0}; d < dimension; d++) {
            std::string id = std::string {"Build_Momentum_"} + std::to_string(d);
            t_build_mom[d] = timer->GetElapsedTime(id);
            id = std::string {"Solve_Momentum_"} + std::to_string(d);
            t_solve_mom[d] = timer->GetElapsedTime(id);
        }
        Timer::ValueType t_enforce_continuity = timer->GetElapsedTime("Enforce_Continuity");
        std::vector<Timer::ValueType> t_build_cont(c_loops+1);
        std::vector<Timer::ValueType> t_solve_cont(c_loops+1);

        for (int l{0}; l <= c_loops; l++) {
            std::string id = std::string {"Build_Continuity_"} + std::to_string(l);
            t_build_cont[l] = timer->GetElapsedTime(id);
            id = std::string("Solve_Continuity_") + std::to_string(l);
            t_solve_cont[l] = timer->GetElapsedTime(id);
        }
        Timer::ValueType t_build{0.}, t_solve{0.};
        Timer::ValueType t_build_cont_tot{0.}, t_solve_cont_tot{0.};
        for (std::size_t d{0}; d < dimension; d++) {
            t_build += t_build_mom[d];
            t_solve += t_solve_mom[d];
        }
        for (int l{0}; l <= c_loops; l++) {
            t_build += t_build_cont[l];
            t_solve += t_solve_cont[l];
            t_solve_cont_tot += t_solve_cont[l];
            t_build_cont_tot += t_build_cont[l];
        }
        Print(dare::Verbosity::Low)
            << "Time estimates\n"
            << "--------------\n"
            << "ID                -> elapsed in s (% of step)\n";
        // Print assembly times
        Print(dare::Verbosity::Low)
            << "Assembly Time     -> "
            << std::setprecision(prec) << t_build
            << " (" << std::setprecision(prec) << t_build / t_flowstep * 100 << " %)\n";
        for (std::size_t d{0}; d < dimension; d++)
            Print(dare::Verbosity::Medium)
                << "   Momentum " << mom_id[d] << "     -> "
                << std::setprecision(prec) << t_build_mom[d]
                << " (" << std::setprecision(prec) << t_build_mom[d] / t_flowstep * 100 << " %)\n";
        Print(dare::Verbosity::Medium)
            << "   Continuity     -> " << std::setprecision(prec) << t_build_cont_tot
            << " (" << std::setprecision(prec) << t_build_cont_tot / t_flowstep * 100 << " %)\n";
        if (c_loops > 1) {
            for (int l{0}; l <= c_loops; l++)
                Print(dare::Verbosity::High)
                    << "      it " << std::to_string(l) << "       -> "
                    << std::setprecision(prec) << t_build_cont[l] << " ("
                    << std::setprecision(prec) << t_build_cont[l] / t_flowstep * 100 << " %)\n";
        }
        // Print solving times
        Print(dare::Verbosity::Low)
            << "Solving Time      -> "
            << std::setprecision(prec) << t_solve
            << "(" << std::setprecision(prec) << t_solve / t_flowstep * 100 << " %)\n";
        for (std::size_t d{0}; d < dimension; d++)
            Print(dare::Verbosity::Medium)
                << "   Momentum " << mom_id[d] << "     -> "
                << std::setprecision(prec) << t_solve_mom[d]
                << " (" << std::setprecision(prec) << t_solve_mom[d] / t_flowstep * 100 << " %)\n";
        Print(dare::Verbosity::Medium)
            << "   Continuity     -> " << std::setprecision(prec) << t_solve_cont_tot
            << " (" << std::setprecision(prec) << t_solve_cont_tot / t_flowstep * 100 << " %)\n";
        if (c_loops > 1) {
            for (int l{0}; l <= c_loops; l++)
                Print(dare::Verbosity::High)
                    << "      it " << std::to_string(l) << "       -> "
                    << std::setprecision(prec) << t_solve_cont[l] << " ("
                    << std::setprecision(prec) << t_solve_cont[l] / t_flowstep * 100 << " %)\n";
        }
        Print(dare::Verbosity::Low)
            << "Momentum Total    -> " << std::setprecision(prec) << t_flowstep - t_enforce_continuity << " ("
            << std::setprecision(prec) << (1.-t_enforce_continuity / t_flowstep) * 100 << " %)\n";
        Print(dare::Verbosity::Low)
            << "Continuity Total  -> " << std::setprecision(prec) << t_enforce_continuity << " ("
            << std::setprecision(prec) << t_enforce_continuity / t_flowstep * 100 << " %)\n";
        Print(dare::Verbosity::Low)
            << "Total Time        -> " << std::setprecision(prec) << t_flowstep << " s" << std::endl;
    }

    dare::ExecutionManager* ex_man;                  //!< reference to execution manager
    DensityVariableType rho;                         //!< mass density
    ViscosityVariableType mu;                        //!< dynamic viscosity
    PorosityVariableType epsilon;                    //!< porosity
    DensityDerivativeMemberType density_derivative;  //!< density derivative (compressible)
    std::unique_ptr<GridVectorType> rho_prev;        //!< density at previous iteration (compressible)

    std::unique_ptr<ContinuityType> continuity;      //!< instance of the continuity equation
    std::array<std::unique_ptr<MomentumType>, dimension> momentum;  //!< momentum equations

    SC dt;                      //!< current time step size
    SC time;                    //!< current simulation time
    TimeStepCounter tstep;      //!< current time step
    SC continuity_tolerance;    //!< convergence tolerance for the continuity loop
    SC max_continuity_defect;   //!< last determined maximum continuity defect
    int max_iterations;         //!< maximum iterations for the continuity loop
    char status;                //!< bitflag for checking the object status
    char status_finalized;      //!< expected status of object when computing a flow step
    UniqueObserverHandle pimpl_dt_obs;  //!< observer handle for updating time and timesteps
    std::set<ObserverType*> observers;  //!< attached observers
    TerminalOutput terminal_output;     //!< prints data to the terminal (should probably be a logger)
    std::unique_ptr<Timer> timer;       //!< timing instance for profiling the execution
};

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_
