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
#ifndef MATRIXSYSTEM_OPERATORS_CARTESIAN_H_
#define MATRIXSYSTEM_OPERATORS_CARTESIAN_H_

namespace dare::Matrix {

template <std::size_t Dim>
Divergence<dare::Grid::Cartesian<Dim>>::Divergence(
    const GridRepresentation& grid, LO ordinal_internal)
    : A(grid.GetFaceArea()) {
    // the ordinal is here for the standard interface, nothing more
    // it's not really required when using the Cartesian grid instance
}

template <std::size_t Dim>
template <typename SC, std::size_t N>
dare::utils::Vector<N, SC>
Divergence<dare::Grid::Cartesian<Dim>>::operator()(
    const dare::Data::FaceValueStencil<GridType, SC, N>& s) const {
    dare::utils::Vector<N, SC> div_v;
    for (std::size_t n{0}; n < N; n++) {
        // Divergence in X
        div_v[n] = A[0] * s.GetValue(Positions::EAST, n);
        div_v[n] -= A[0] * s.GetValue(Positions::WEST, n);

        // Divergence in Y
        if constexpr (Dim > 1) {
            div_v[n] += A[1] * s.GetValue(Positions::NORTH, n);
            div_v[n] -= A[1] * s.GetValue(Positions::SOUTH, n);
        }

        // Divergence in Z
        if constexpr (Dim > 2) {
            div_v[n] += A[2] * s.GetValue(Positions::TOP, n);
            div_v[n] -= A[2] * s.GetValue(Positions::BOTTOM, n);
        }
    }
    return div_v;
}

template <std::size_t Dim>
template <typename SC, std::size_t N>
dare::Data::CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
Divergence<dare::Grid::Cartesian<Dim>>::operator()(
    const dare::Data::FaceMatrixStencil<GridType, SC, N>& s) const {
    dare::Data::CenterMatrixStencil<GridType, SC, N> s_c;
    for (std::size_t n{0}; n < N; n++) {
        // Divergence in X
        SC coef_c = s.GetValueCenter(Positions::EAST, n);
        SC coef_f = s.GetValueNeighbor(Positions::EAST, n);
        s_c.GetValue(Positions::EAST, n) = A[0] * coef_f;
        s_c.GetValue(Positions::CENTER, n) = A[0] * coef_c;
        s_c.GetRHS(n) = A[0] * s.GetRHS(Positions::EAST, n);

        coef_c = s.GetValueCenter(Positions::WEST, n);
        coef_f = s.GetValueNeighbor(Positions::WEST, n);
        s_c.GetValue(Positions::WEST, n) = -A[0] * coef_f;
        s_c.GetValue(Positions::CENTER, n) -= A[0] * coef_c;
        s_c.GetRHS(n) -= A[0] * s.GetRHS(Positions::WEST, n);

        // Divergence in Y
        if constexpr (Dim > 1) {
            SC coef_c = s.GetValueCenter(Positions::NORTH, n);
            SC coef_f = s.GetValueNeighbor(Positions::NORTH, n);
            s_c.GetValue(Positions::NORTH, n) = A[1] * coef_f;
            s_c.GetValue(Positions::CENTER, n) += A[1] * coef_c;
            s_c.GetRHS(n) += A[1] * s.GetRHS(Positions::NORTH, n);

            coef_c = s.GetValueCenter(Positions::SOUTH, n);
            coef_f = s.GetValueNeighbor(Positions::SOUTH, n);
            s_c.GetValue(Positions::SOUTH, n) = -A[1] * coef_f;
            s_c.GetValue(Positions::CENTER, n) -= A[1] * coef_c;
            s_c.GetRHS(n) -= A[1] * s.GetRHS(Positions::SOUTH, n);
        }

        // Divergence in Z
        if constexpr (Dim > 2) {
            SC coef_c = s.GetValueCenter(Positions::TOP, n);
            SC coef_f = s.GetValueNeighbor(Positions::TOP, n);
            s_c.GetValue(Positions::TOP, n) = A[2] * coef_f;
            s_c.GetValue(Positions::CENTER, n) += A[2] * coef_c;
            s_c.GetRHS(n) += A[2] * s.GetRHS(Positions::TOP, n);

            coef_c = s.GetValueCenter(Positions::BOTTOM, n);
            coef_f = s.GetValueNeighbor(Positions::BOTTOM, n);
            s_c.GetValue(Positions::BOTTOM, n) = -A[2] * coef_f;
            s_c.GetValue(Positions::CENTER, n) -= A[2] * coef_c;
            s_c.GetRHS(n) -= A[2] * s.GetRHS(Positions::BOTTOM, n);
        }
    }
    return s_c;
}

}  // end namespace dare::Matrix


#endif  // MATRIXSYSTEM_OPERATORS_CARTESIAN_H_
