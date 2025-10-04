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

namespace dare {

template <typename G, typename BInf, typename PInf, typename NInf>
ProjectionMethod<G, BInf, PInf, NInf>::ProjectionMethod()
    : ex_man{nullptr},
      dt{0.},
      time{0.},
      tstep{0},
      continuity_tolerance{1e-14},
      status{0},
      status_finalized{rho_init | mu_init | epsilon_init | beta_im_init | beta_ex_init | density_derivative_init},
      status_build{0},
      terminal_output(dimension),
      params(detail::GetDefaultParameterListPM()) {
    if constexpr (dare::is_none_v<PorosityVariableType>)
        status |= epsilon_init;
    if constexpr (dare::is_none_v<ImplicitForceVariableType>)
        status |= beta_im_init;
    if constexpr (dare::is_none_v<ExplicitForceVariableType>)
        status |= beta_ex_init;
    if constexpr (dare::is_none_v<DensityDerivativeMemberType>)
        status |= density_derivative_init;

    // method specific check for consistent compile time information
    free_compile_time_check(this);
}

template <typename G, typename BInf, typename PInf, typename NInf>
ProjectionMethod<G, BInf, PInf, NInf>::~ProjectionMethod() {}

template <typename G, typename BInf, typename PInf, typename NInf>
template <TimeStepper T, typename BCContinuity, typename... Args>
void ProjectionMethod<G, BInf, PInf, NInf>::Initialize(
    G* grid, T* tstep, BCContinuity bc_continuity, Args... bc_momentum) {
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

template <typename G, typename BInf, typename PInf, typename NInf>
template <TimeStepper T, typename BCContinuity, typename... Args>
void ProjectionMethod<G, BInf, PInf, NInf>::Initialize(std::unique_ptr<G>& grid, T* tstep, BCContinuity bc_continuity, Args... bc_momentum) {  // NOLINT
    Initialize(grid.get(), tstep, bc_continuity, bc_momentum...);
}

template <typename G, typename BInf, typename PInf, typename NInf>
template <typename... Args>
void ProjectionMethod<G, BInf, PInf, NInf>::InitializeMomentum(std::size_t dim, Args&&... args) {
    momentum[dim] = std::make_unique<MomentumType>(args...);
}

template <typename G, typename BInf, typename PInf, typename NInf>
template <typename... Args>
void ProjectionMethod<G, BInf, PInf, NInf>::InitializeContinuity(Args&&... args) {
    continuity = std::make_unique<ContinuityType>(args...);
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::AddImplicitForce(ImplicitForceVariableType f) {
    // add to continuity
    if (!IsInitialized()) {
        ex_man->Terminate(__func__, "Cannot add force terms prior to initialization");
    }
    if constexpr (!dare::is_none_v<ImplicitForceVariableType>) {
        continuity->GetCustomMember()->beta_im.emplace(f);
    }
    status |= beta_im_init;
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::AddExplicitForce(ExplicitForceVariableType f, std::size_t dim) {
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
            all_init &= !m->GetCustomMember()->beta_ex.empty() ? beta_ex_init : 0;
        status |= (beta_ex_init & all_init);
    }
}

template <typename G, typename BInf, typename PInf, typename NInf>
bool ProjectionMethod<G, BInf, PInf, NInf>::Attach(ObserverType* o) {
    auto [pos, success] = observers.emplace(o);
    return success;
}

template <typename G, typename BInf, typename PInf, typename NInf>
bool ProjectionMethod<G, BInf, PInf, NInf>::Detach(ObserverType* o) {
    return (observers.erase(o) > 0U);
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::Notify(StateChange property) {
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

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::CopyToOld() {
    for (std::size_t i{0}; i < dimension; i++)
        GetMomentum(i)->GetField()->CopyDataVectorsToOldTimeStep();
    GetContinuity()->GetField()->CopyDataVectorsToOldTimeStep();
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::SolveFlowField() {
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
        GetMomentum(dare::ZERO)->PreStep();
        BuildMomentum(dare::ZERO);
        status_build |= A_mom_0_build;
        auto [success, iter] = SolveMomentum(dare::ZERO);
        GetMomentum(dare::ZERO)->PostStep();
        terminal_output.PrintMomentum(0, iter, success, this);
    }
    // add some output here
    if constexpr (dimension > 1) {
        GetMomentum(dare::ONE)->PreStep();
        BuildMomentum(dare::ONE);
        status_build |= A_mom_1_build;
        auto [success, iter] = SolveMomentum(dare::ONE);
        GetMomentum(dare::ONE)->PostStep();
        // add some output here
        terminal_output.PrintMomentum(1, iter, success, *this);
    }
    if constexpr (dimension > 2) {
        GetMomentum(dare::TWO)->PreStep();
        BuildMomentum(dare::TWO);
        status_build |= A_mom_2_build;
        auto [success, iter] = SolveMomentum(dare::TWO);
        GetMomentum(dare::TWO)->PostStep();
        terminal_output.PrintMomentum(2, iter, success, *this);
    }

    // enforce continuity
    // first iteration we add the implicit force term
    // or in the case of no additional term this is just
    // the very first iteration
    StartProfiling("Enforce_Continuity");
    int iteration = 0;
    int max_iterations = params.get<int>("continuity: Newton iterations max");
    ComputeDefect();
    terminal_output.PrintInitialDefect(*this);
    for (; iteration < max_iterations; iteration++) {
        if constexpr (compressible)
            *rho_prev = rho->GetGridVector(0);

        BuildContinuity(iteration);
        status_build |= A_p_build;
        auto [success, iter] = SolveContinuity(iteration);

        if (!success) {
            dare::Print(dare::Verbosity::Low) << "Continuity system failed to converge after "
                                              << iter << " matrix-solver iterations" << std::endl;
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

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::SetDensity(DensityVariableType d) {
    rho = d;
    status |= rho_init;
    if constexpr (compressible) {
        static_assert(dare::is_field_v<std::remove_cv_t<std::remove_pointer_t<DensityVariableType>>>,
                      "In the compressible case, the density needs to be a field!");
        rho_prev = std::make_unique<GridVectorType>("rho_prev", rho->GetGridRepresentation());
        *rho_prev = rho->GetDataVector();
    }
}


template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::SetViscosity(ViscosityVariableType v) {
    mu = v;
    status |= mu_init;
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::SetPorosity(PorosityVariableType p) {
    epsilon = p;
    status |= epsilon_init;
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::SetDensityDerivative(DensityDerivativeMemberType dd) {
    density_derivative = std::move(dd);
    status |= density_derivative_init;
}

template <typename G, typename BInf, typename PInf, typename NInf>
const typename ProjectionMethod<G, BInf, PInf, NInf>::GridVectorType*
ProjectionMethod<G, BInf, PInf, NInf>::GetDensityPreviousIteration() const {
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

template <typename G, typename BInf, typename PInf, typename NInf>
typename ProjectionMethod<G, BInf, PInf, NInf>::SC
ProjectionMethod<G, BInf, PInf, NInf>::GetDensity(Index ind, std::size_t time_level) const {
    using DType = std::remove_cv_t<std::remove_pointer_t<DensityVariableType>>;
    if constexpr (dare::is_field_v<DType>) {
        return rho->GetDataVector(time_level).At(ind, 0);
    } else if constexpr(std::is_arithmetic_v<DType>) {
        return rho;
    } else {
        static_assert(dare::always_false<DType>, "Cannot handle this type");
    }
}

template <typename G, typename BInf, typename PInf, typename NInf>
typename ProjectionMethod<G, BInf, PInf, NInf>::SC
ProjectionMethod<G, BInf, PInf, NInf>::GetViscosity(Index ind, std::size_t time_level) const {
    using VType = std::remove_cv_t<std::remove_pointer_t<ViscosityVariableType>>;
    if constexpr (dare::is_field_v<VType>) {
        return mu->GetDataVector(time_level).At(ind, 0);
    } else if constexpr (std::is_arithmetic_v<VType>) {
        return mu;
    } else {
        static_assert(dare::always_false<VType>, "Cannot handle this type");
    }
}

template <typename G, typename BInf, typename PInf, typename NInf>
typename ProjectionMethod<G, BInf, PInf, NInf>::SC
ProjectionMethod<G, BInf, PInf, NInf>::GetPorosity(Index ind, std::size_t time_level) const {
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

template <typename G, typename BInf, typename PInf, typename NInf>
template <dare::NaturalNumber Direction>
void ProjectionMethod<G, BInf, PInf, NInf>::BuildMomentum(Direction dir) {
    std::string id = std::string("Build_Momentum_") + std::to_string(dir);
    StartProfiling(id);
    free_pm_build_momentum(this, dir);
    StopProfiling(id);
}

template <typename G, typename BInf, typename PInf, typename NInf>
template <dare::NaturalNumber Direction>
std::pair<bool, int> ProjectionMethod<G, BInf, PInf, NInf>::SolveMomentum(Direction dir) {
    using SProp = typename MomentumType::SolverNumericalPropertiesType;
    std::string id = std::string("Solve_Momentum_") + std::to_string(dir);
    std::string sprop = "momentum[" + std::to_string(dir) + "]: solver properties";

    if (!params.isParameter(sprop) || params.isType<dare::None>(sprop)) {
        Print(dare::Verbosity::High) << "Did not find following key in the properties: " << sprop;
        sprop = "momentum: solver properties";
        Print(dare::Verbosity::High) << " Testing for " << sprop << " instead" << std::endl;
    }
    if (!params.isParameter(sprop) || params.isType<dare::None>(sprop))
        ex_man->Terminate(__func__,
            "Did not find the solver properties for the momentum in direction " + std::to_string(dir)
            + ". Possible missing key: 'momentum: solver properties'");
    GetMomentum(dir)->SetSolverNumericalProperties(params.template get<SProp>(sprop));

    StartProfiling(id);
    auto ret = free_pm_solve_momentum(this, dir);
    StopProfiling(id);
    return ret;
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::BuildContinuity(int iteration) {
    std::string id = std::string("Build_Continuity_") + std::to_string(iteration);
    StartProfiling(id);
    free_pm_build_continuity(this, iteration);
    StopProfiling(id);
}


template <typename G, typename BInf, typename PInf, typename NInf>
std::pair<bool, int> ProjectionMethod<G, BInf, PInf, NInf>::SolveContinuity(int iteration) {
    using SProp = typename ContinuityType::SolverNumericalPropertiesType;
    std::string id = std::string("Solve_Continuity_") + std::to_string(iteration);
    std::string sprop = "continuity: solver properties";
    if (!params.isParameter(sprop) || params.isType<dare::None>(sprop))
        ex_man->Terminate(__func__,
            "Did not find the solver properties for the continuity with the key 'continuity: solver properties'.");
    GetContinuity()->SetSolverNumericalProperties(params.template get<SProp>(sprop));

    StartProfiling(id);
    auto ret = free_pm_solve_continuity(this, iteration);
    StopProfiling(id);
    return ret;
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::StartProfiling(std::string id) {
    if (timer)
        timer->Tic(id);
}

template <typename G, typename BInf, typename PInf, typename NInf>
typename dare::Timer::ValueType
ProjectionMethod<G, BInf, PInf, NInf>::StopProfiling(std::string id) {
    if (timer)
        return timer->Toc(id);
    return 0.;
}

template <typename G, typename BInf, typename PInf, typename NInf>
void ProjectionMethod<G, BInf, PInf, NInf>::PrintProfiling(int c_loops) const {
    if (!timer)
        return;

    char mom_id[] = {'X', 'Y', 'Z'};

    const int prec{3};  //!< precision of the output
    Timer::ValueType t_flowstep = timer->GetElapsedTime("FlowStep");
    std::array<Timer::ValueType, dimension> t_build_mom, t_solve_mom;
    for (std::size_t d{0}; d < dimension; d++) {
        std::string id = std::string("Build_Momentum_") + std::to_string(d);
        t_build_mom[d] = timer->GetElapsedTime(id);
        id = std::string("Solve_Momentum_") + std::to_string(d);
        t_solve_mom[d] = timer->GetElapsedTime(id);
    }
    Timer::ValueType t_enforce_continuity = timer->GetElapsedTime("Enforce_Continuity");
    std::vector<Timer::ValueType> t_build_cont(c_loops + 1);
    std::vector<Timer::ValueType> t_solve_cont(c_loops + 1);

    for (int l{0}; l <= c_loops; l++) {
        std::string id = std::string("Build_Continuity_") + std::to_string(l);
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
        << std::setprecision(prec) << (1. - t_enforce_continuity / t_flowstep) * 100 << " %)\n";
    Print(dare::Verbosity::Low)
        << "Continuity Total  -> " << std::setprecision(prec) << t_enforce_continuity << " ("
        << std::setprecision(prec) << t_enforce_continuity / t_flowstep * 100 << " %)\n";
    Print(dare::Verbosity::Low)
        << "Total Time        -> " << std::setprecision(prec) << t_flowstep << " s" << std::endl;
}

}  // namespace dare
