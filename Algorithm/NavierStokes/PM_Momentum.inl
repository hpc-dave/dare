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
PMMomentum<Grid, BS, CM>::PMMomentum(const std::string& name,
                                     GridRepresentation grid,
                                     dare::mpi::ExecutionManager* ex_man,
                                     std::size_t num_tsteps,
                                     BS bc_strat)
    : grep(grid),
      exec_man(ex_man),
      data(name, grid, num_tsteps),
      boundary_strategy(std::move(bc_strat)) {
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::PreStep() {
    for (auto it : pre_step_strategy)
        (*it)(this);
}

template <typename Grid, typename BS, typename CM>
template <typename BuildStrategy>
void PMMomentum<Grid, BS, CM>::Build(BuildStrategy build) {}

template <typename Grid, typename BS, typename CM>
std::pair<bool, int> PMMomentum<Grid, BS, CM>::Solve(SolverPropertyType sprop, PreconditionerPropertyType mprop) {
    return {false, -1};
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::UpdateBoundaries() {
    boundary_strategy(&data);
    data.ExchangeHaloCells();
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::PostStep() {
    for (auto it : post_step_strategy)
        (*it)(this);
}

template <typename Grid, typename BS, typename CM>
BS* PMMomentum<Grid, BS, CM>::GetBoundaryStrategy() {
    return &boundary_strategy;
}

template <typename Grid, typename BS, typename CM>
const BS& PMMomentum<Grid, BS, CM>::GetBoundaryStrategy() const {
    return boundary_strategy;
}

template <typename Grid, typename BS, typename CM>
PMMomentum<Grid, BS, CM>::FieldType* PMMomentum<Grid, BS, CM>::GetField() {
    return &data;
}

template <typename Grid, typename BS, typename CM>
const PMMomentum<Grid, BS, CM>::FieldType& PMMomentum<Grid, BS, CM>::GetField() const {
    return data;
}

template <typename Grid, typename BS, typename CM>
CM* PMMomentum<Grid, BS, CM>::GetCustomMember() {
    return &custom_member;
}

template <typename Grid, typename BS, typename CM>
const CM& PMMomentum<Grid, BS, CM>::GetCustomMember() const {
    return custom_member;
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::AddPreStepStrategy(std::function<void(SelfType*)> f) {
    pre_step_strategy.insert(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::SetPreStepStrategy(std::function<void(SelfType*)> f) {
    ClearPreStepStrategy();
    AddPreStepStrategy(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::ClearPreStepStrategy() {
    pre_step_strategy.clear();
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::AddPostStepStrategy(std::function<void(SelfType*)> f) {
    post_step_strategy.insert(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::SetPostStepStrategy(std::function<void(SelfType*)> f) {
    ClearPostStepStrategy();
    AddPostStepStrategy(std::move(f));
}

template <typename Grid, typename BS, typename CM>
void PMMomentum<Grid, BS, CM>::ClearPostStepStrategy() {
    post_step_strategy.clear();
}
}  // namespace dare::algorithm
