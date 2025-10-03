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
template<typename BC>
PMContinuity<Grid, BS, CM>::PMContinuity(const std::string& name,
                                         GridRepresentation grid,
                                         std::size_t num_tsteps,
                                         BC bc_strat)
    : dare::GenericEquation<Grid, BS, CM>(name, std::move(grid), num_tsteps, bc_strat),
      defect("defect", grid, 1),
      dP("dP", grid, 1) {
}

template <typename Grid, typename BS, typename CM>
typename PMContinuity<Grid, BS, CM>::FieldType* PMContinuity<Grid, BS, CM>::GetPressure() {
    return this->GetField();
}

template <typename Grid, typename BS, typename CM>
const typename PMContinuity<Grid, BS, CM>::FieldType& PMContinuity<Grid, BS, CM>::GetPressure() const {
    return this->GetField();
}

template <typename Grid, typename BS, typename CM>
typename PMContinuity<Grid, BS, CM>::FieldType* PMContinuity<Grid, BS, CM>::GetDefect() {
    return &defect;
}

template <typename Grid, typename BS, typename CM>
const typename PMContinuity<Grid, BS, CM>::FieldType& PMContinuity<Grid, BS, CM>::GetDefect() const {
    return defect;
}

template <typename Grid, typename BS, typename CM>
typename PMContinuity<Grid, BS, CM>::FieldType* PMContinuity<Grid, BS, CM>::GetdP() {
    return &dP;
}

template <typename Grid, typename BS, typename CM>
const typename PMContinuity<Grid, BS, CM>::FieldType& PMContinuity<Grid, BS, CM>::GetdP() const {
    return dP;
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::UpdatedPBoundaries() {
    (*this->GetBoundaryStrategy())(&dP);
    dP.ExchangeHaloCells();
}

template <typename Grid, typename BS, typename CM>
void PMContinuity<Grid, BS, CM>::UpdatePressureBoundaries() {
    this->GetField()->ExchangeHaloCells();
}

}  // namespace dare
