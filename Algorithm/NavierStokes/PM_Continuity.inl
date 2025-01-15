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

namespace dare::algorithm {

template <typename Grid, typename BS, typename CM>
PMContinuity<Grid, BS, CM>::PMContinuity(const std::string& name,
                                         GridRepresentation grid,
                                         dare::mpi::ExecutionManager* ex_man,
                                         std::size_t num_tsteps,
                                         BS bc_strat)
    : grep(grid),
      exec_man(ex_man),
      pressure(name, grid, num_tsteps),
      defect("defect", grid, 1),
      dP("dP", grid, 1),
      boundary_strategy(std::move(bc_strat)) {
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::PreStep() {
    for (auto it : pre_step_strategy)
        (*it)(this);
}

template <typename Grid, typename BS, typename CM>
template <typename BuildStrategy>
void PMContinuity<Grid, BS, CM>::Build(BuildStrategy build) {}

template <typename Grid, typename BS, typename CM>
std::pair<bool, int> PMContinuity<Grid, BS, CM>::Solve(SolverPropertyType sprop,
                                                       PreconditionerPropertyType mprop) {
    return {false, -1};
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::UpdateBoundaries() {
    boundary_strategy(&data);
    data.ExchangeHaloCells();
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::PostStep() {
    for (auto it : post_step_strategy)
        (*it)(this);
}

template <typename Grid, typename BS, typename CM>
BS* PMContinuity<Grid, BS, CM>::GetBoundaryStrategy() {
    return &boundary_strategy;
}

template <typename Grid, typename BS, typename CM>
const BS& PMContinuity<Grid, BS, CM>::GetBoundaryStrategy() const {
    return boundary_strategy;
}

template <typename Grid, typename BS, typename CM>
PMContinuity<Grid, BS, CM>::FieldType* PMContinuity<Grid, BS, CM>::GetPressure() {
    return &pressure;
}

template <typename Grid, typename BS, typename CM>
const PMContinuity<Grid, BS, CM>::FieldType& PMContinuity<Grid, BS, CM>::GetPressure() const {
    return pressure;
}

template <typename Grid, typename BS, typename CM>
PMContinuity<Grid, BS, CM>::FieldType* PMContinuity<Grid, BS, CM>::GetDefect() {
    return &defect;
}

template <typename Grid, typename BS, typename CM>
const PMContinuity<Grid, BS, CM>::FieldType& PMContinuity<Grid, BS, CM>::GetDefect() const {
    return defect;
}

template <typename Grid, typename BS, typename CM>
PMContinuity<Grid, BS, CM>::FieldType* PMContinuity<Grid, BS, CM>::GetdP() {
    return &dP;
}

template <typename Grid, typename BS, typename CM>
const PMContinuity<Grid, BS, CM>::FieldType& PMContinuity<Grid, BS, CM>::GetdP() const {
    return dP;
}
template <typename Grid, typename BS, typename CM>
CM* PMContinuity<Grid, BS, CM>::GetCustomMember() {
    return &custom_member;
}

template <typename Grid, typename BS, typename CM>
const CM& PMContinuity<Grid, BS, CM>::GetCustomMember() const {
    return custom_member;
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::AddPreStepStrategy(std::function<void(SelfType*)> f) {
    pre_step_strategy.insert(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::SetPreStepStrategy(std::function<void(SelfType*)> f) {
    ClearPreStepStrategy();
    AddPreStepStrategy(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::ClearPreStepStrategy() {
    pre_step_strategy.clear();
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::AddPostStepStrategy(std::function<void(SelfType*)> f) {
    post_step_strategy.insert(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::SetPostStepStrategy(std::function<void(SelfType*)> f) {
    ClearPostStepStrategy();
    AddPostStepStrategy(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::ClearPostStepStrategy() {
    post_step_strategy.clear();
}
}  // namespace dare::algorithm
