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
#ifndef MATRIXSYSTEM_OPERATORS_H_
#define MATRIXSYSTEM_OPERATORS_H_

namespace dare::Matrix {

/*!
 * @brief basic dummy for SFINAE for a gradient operator
 * @tparam Grid type of grid
 *
 * The purpose of this class is to enable easy translation of mathematical formulae and
 * equations into any matrix building procedure. Any specialization of this class should allow
 * following code examples to work
 *
 * @code{.cpp}
 * // assuming an instance 'grid' of the Grid exists
 *
 * // Initialization may deviate per specialization, however the maximum requirement
 * // should be the one where you specify an Index of the relevant cell
 * Gradient<Grid> grad(grid, Index ind);
 *
 * // determine the gradient values at the faces from a field with N components
 * FaceValueStencil<Grid, N> gradients_at_faces = grad(field);
 *
 * // if a field has multiple components, the gradient should also be able to return
 * // the gradient of only one component
 * FaceValueStencil<Grid, 1> gradient_single_at_faces = grad(field, n);  // n here is the component ID
 *
 * // For adding to a matrix, using MatrixBlock instance mb, following interface must exist
 * FaceMatrixStencil<Grid, N> gradient_for_matrix = grad(mb);
 *
 * @endcode
 */
template <typename Grid>
class Gradient {
};

/*!
 * @brief basic dummy for SFINAE for a divergence operator
 * @tparam Grid type of grid
 *
 * The purpose of this class is to enable easy translation of mathematical formulae and
 * equations into any matrix building procedure. Any specialization of this class should allow
 * following code examples to work
 *
 * @code{.cpp}
 * // assuming an instance 'grid' of the Grid exists
 *
 * // Initialization may deviate per specialization, however the maximum requirement
 * // should be the one where you specify an Index of the relevant cell
 * Divergence<Grid> div(grid, Index ind);
 *
 * // determine the divergence from values at the faces with N components, SC here is the scalar type
 * FaceValueStencil values_at_faces;
 * dare::utils::Vector<N, SC> div_values = div(values_at_faces);
 *
 * // For adding to a matrix, following interface must exist
 * FaceMatrixValues matrix_entries_at_faces;
 * CenterMatrixStencil<Grid, N> div_for_matrix = div(matrix_entries_at_faces);
 *
 * // Finally, the whole thing should be able to look as follows with the gradient operator
 * dare::utils::Vector<N, SC> div_values = div(rho * grad(field));
 *
 * // or
 * MatrixBlock<Grid, N> mb;
 * mb += div(rho * grad());
 *
 * @endcode
 */
template<typename Grid>
class Divergence {
};

}  // end namespace dare::Matrix

#endif  // MATRIXSYSTEM_OPERATORS_H_