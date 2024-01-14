/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder

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

#ifndef MATH_INTERPOLATION_H_
#define MATH_INTERPOLATION_H_

#include "../Data/GridVector.h"
#include "Pow.h"

namespace dare::math {

/*!
 * @brief dummy for SFINAE
 * @tparam GridType type of grid
 * @tparam SC type of scalar
 * @tparam N number of components
 * @param target target grid
 * @param ind_target indices of the target cell
 * @param field 
 * @param n
 * @return 
 */
template<typename GridType, typename SC, std::size_t N>
SC InterpolateToFace(const typename GridType::GridRepresentation& target,
                     const typename GridType::Index& ind_target,
                     const typename GridType::NeighborID face,
                     const Data::GridVector<GridType, SC, N>& field,
                     std::size_t n) {
}

}  // end namespace dare::math

#endif  // MATH_INTERPOLATION_H_
