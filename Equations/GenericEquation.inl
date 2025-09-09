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

template <typename Grid, typename BS, typename CM>
template <typename BC>
GenericEquation<Grid, BS, CM>::GenericEquation(const std::string& name,
                                               GridRepresentation grid,
                                               std::size_t num_tsteps,
                                               BC bc_strat)
    : grep(std::move(grid)),
      exec_man(nullptr),
      field(name, grid, num_tsteps),
      boundary_strategy(bc_strat) {
    exec_man = grep.GetExecutionManager();
    if (!exec_man) {
        ERROR << "Execution Manager is nullptr!" << ERROR_CLOSE;
    }
    matrix_system.Initialize(exec_man);
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::PreStep() {
    for (auto& f : pre_step_strategy)
        f(this);
    Notify(StateChange::PreStep);
}

template <typename Grid, typename BS, typename CM>
template <typename BuildStrategy>
void GenericEquation<Grid, BS, CM>::Build(BuildStrategy build_lambda) {
    const bool rebuild{false};
    Notify(StateChange::PreBuild);
    matrix_system.Build(grep, field.GetDataVector(), build_lambda, rebuild);
    Notify(StateChange::PostBuild);
}

template <typename Grid, typename BS, typename CM>
template <typename BuildStrategy>
void GenericEquation<Grid, BS, CM>::UpdateRhs(BuildStrategy build_lambda) {
    Notify(StateChange::PreBuildRhs);
    matrix_system.SetB(grep, field.GetDataVector(), build_lambda);
    Notify(StateChange::PostBuildRhs);
}

template <typename Grid, typename BS, typename CM>
template <typename UpdateStrategy>
std::pair<bool, int> GenericEquation<Grid, BS, CM>::Solve(UpdateStrategy strat_update, bool build_prec) {
    Notify(StateChange::PreSolve);
    MatrixSolverType solver;
    if (build_prec || matrix_system.GetM().is_null())
        matrix_system.GetM() = solver.BuildPreconditioner(solver_prop, matrix_system.GetA());

    auto ret = solver.Solve(solver_prop,
                            matrix_system.GetM(),
                            matrix_system.GetA(),
                            matrix_system.GetX(),
                            matrix_system.GetB());

    strat_update(*this, matrix_system);
    Notify(StateChange::PostSolve);
    return {ret == Belos::ReturnType::Converged, solver.GetNumIterations()};
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::UpdateBoundaries() {
    boundary_strategy(&field);
    field.ExchangeHaloCells();
    Notify(StateChange::UpdateBoundaries);
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::PostStep() {
    for (auto& f : post_step_strategy)
        f(this);
    Notify(StateChange::PostStep);
}

template <typename Grid, typename BS, typename CM>
typename GenericEquation<Grid, BS, CM>::GridRepresentation*
GenericEquation<Grid, BS, CM>::GetGridRepresentation() {
    return &grep;
}

template <typename Grid, typename BS, typename CM>
const typename GenericEquation<Grid, BS, CM>::GridRepresentation&
GenericEquation<Grid, BS, CM>::GetGridRepresentation() const {
    return grep;
}


template <typename Grid, typename BS, typename CM>
BS* GenericEquation<Grid, BS, CM>::GetBoundaryStrategy() {
    return &boundary_strategy;
}

template <typename Grid, typename BS, typename CM>
const BS& GenericEquation<Grid, BS, CM>::GetBoundaryStrategy() const {
    return boundary_strategy;
}

template <typename Grid, typename BS, typename CM>
GenericEquation<Grid, BS, CM>::FieldType* GenericEquation<Grid, BS, CM>::GetField() {
    return &field;
}

template <typename Grid, typename BS, typename CM>
const GenericEquation<Grid, BS, CM>::FieldType& GenericEquation<Grid, BS, CM>::GetField() const {
    return field;
}

template <typename Grid, typename BS, typename CM>
CM* GenericEquation<Grid, BS, CM>::GetCustomMember() {
    return &custom_member;
}

template <typename Grid, typename BS, typename CM>
const CM& GenericEquation<Grid, BS, CM>::GetCustomMember() const {
    return custom_member;
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::SetSolverNumericalProperties(const SolverNumericalPropertiesType& prop) {
    solver_prop = prop;
}

template <typename Grid, typename BS, typename CM>
const GenericEquation<Grid, BS, CM>::SolverNumericalPropertiesType&
GenericEquation<Grid, BS, CM>::GetSolverNumericalProperties() const {
    return solver_prop;
}

template <typename Grid, typename BS, typename CM>
GenericEquation<Grid, BS, CM>::MatrixSystemType*
GenericEquation<Grid, BS, CM>::GetMatrixSystem() {
    return &matrix_system;
}

template <typename Grid, typename BS, typename CM>
const GenericEquation<Grid, BS, CM>::MatrixSystemType&
GenericEquation<Grid, BS, CM>::GetMatrixSystem() const {
    return matrix_system;
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::AddPreStepStrategy(std::function<void(SelfType*)> f) {
    pre_step_strategy.insert(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::SetPreStepStrategy(std::function<void(SelfType*)> f) {
    ClearPreStepStrategy();
    AddPreStepStrategy(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::ClearPreStepStrategy() {
    pre_step_strategy.clear();
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::AddPostStepStrategy(std::function<void(SelfType*)> f) {
    post_step_strategy.insert(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::SetPostStepStrategy(std::function<void(SelfType*)> f) {
    ClearPostStepStrategy();
    AddPostStepStrategy(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::ClearPostStepStrategy() {
    post_step_strategy.clear();
}

template <typename Grid, typename BS, typename CM>
bool GenericEquation<Grid, BS, CM>::Attach(ObserverType* o) {
    auto [pos, success] = observers.emplace(o);
    return success;
}

template <typename Grid, typename BS, typename CM>
bool GenericEquation<Grid, BS, CM>::Detach(ObserverType* o) {
    return (observers.erase(o) > 0U);
}

template <typename Grid, typename BS, typename CM>
void GenericEquation<Grid, BS, CM>::Notify(StateChange property) {
    for (auto iter = observers.begin(); iter != observers.end();) {
        auto const pos = iter++;
        (*pos)->Update(*this, property);
    }
}
}  // namespace dare
