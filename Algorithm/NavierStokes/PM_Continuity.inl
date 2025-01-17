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
    : Matrix::GenericEquation<Grid, BS, CM>(name, std::move(grid), ex_man, num_tsteps, std::move(bc_strat)),
      defect("defect", grid, 1),
      dP("dP", grid, 1) {
}

template <typename Grid, typename BS, typename CM>
PMContinuity<Grid, BS, CM>::FieldType* PMContinuity<Grid, BS, CM>::GetPressure() {
    return this->GetField();
}

template <typename Grid, typename BS, typename CM>
const PMContinuity<Grid, BS, CM>::FieldType& PMContinuity<Grid, BS, CM>::GetPressure() const {
    return this->GetField();
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

}  // namespace dare::algorithm
