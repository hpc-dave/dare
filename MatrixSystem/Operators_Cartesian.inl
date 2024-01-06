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

namespace dare::Matrix {

template <std::size_t Dim, typename LO, typename GO, typename SC>
Gradient<dare::Grid::Cartesian<Dim, LO, GO, SC>>::Gradient(const GridRepresentation& _grid, LO ordinal_int)
    : grid(&_grid), ordinal_internal(ordinal_int) {
    for (auto& e : dn_r)
        e = static_cast<SC>(1.);
    dn_r /= grid->GetDistances();
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
template <typename O, std::size_t N>
dare::Data::FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
Gradient<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator() (
    const MatrixBlock<GridType, O, SC, N>&) const {
    dare::Data::FaceMatrixStencil<GridType, N> s;
    s.SetValues(Positions::WEST, 0, -dn_r[0],  dn_r[0]);
    s.SetValues(Positions::EAST, 0,  dn_r[0], -dn_r[0]);
    if constexpr (Dim > 1) {
        s.SetValues(Positions::SOUTH, 0, -dn_r[1],  dn_r[1]);
        s.SetValues(Positions::NORTH, 0,  dn_r[1], -dn_r[1]);
    }
    if constexpr (Dim > 2) {
        s.SetValues(Positions::BOTTOM, 0, -dn_r[2],  dn_r[2]);
        s.SetValues(Positions::TOP,    0,  dn_r[2], -dn_r[2]);
    }

    for (std::size_t n{1}; n < N; n++) {
        s.GetDataCenter()[n] = s.GetDataCenter()[n - 1];
        s.GetDataNeighbor()[n] = s.GetDataNeighbor()[n - 1];
    }
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
template <std::size_t N>
dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
Gradient<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator()(
    const dare::Data::GridVector<GridType, SC, N>& field) const {
    dare::Data::FaceValueStencil<GridType, N> s_f;
    const LO ordinal_local = grid->MapInternalToLocal(ordinal_internal);
    Index ind = grid->MapOrdinalToIndexLocal(ordinal_local);
    for (std::size_t n{0}; n < N; n++) {
        Index ind_nb(ind);
        ind_nb.i() = ind.i() - 1;
        s_f.SetValue(Positions::WEST, n,
                     (field.At(ind, n) - field.At(ind_nb, n)) * dn_r[0]);
        ind_nb.i() = ind.i() + 1;
        s_f.SetValue(Positions::EAST, n,
                     (field.At(ind_nb, n) - field.At(ind, n)) * dn_r[0]);
        if constexpr (Dim > 1) {
            ind_nb = ind;
            ind_nb.j() = ind.j() - 1;
            s_f.SetValue(Positions::SOUTH, n,
                         (field.At(ind, n) - field.At(ind_nb, n)) * dn_r[1]);
            ind_nb.j() = ind.j() + 1;
            s_f.SetValue(Positions::NORTH, n,
                         (field.At(ind_nb, n) - field.At(ind, n)) * dn_r[1]);
        }
        if constexpr (Dim > 2) {
            ind_nb = ind;
            ind_nb.k() = ind.k() - 1;
            s_f.SetValue(Positions::BOTTOM, n,
                         (field.At(ind, n) - field.At(ind_nb, n)) * dn_r[2]);
            ind_nb.k() = ind.k() + 1;
            s_f.SetValue(Positions::TOP, n,
                         (field.At(ind_nb, n) - field.At(ind, n)) * dn_r[2]);
        }
    }
    return s_f;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
template <std::size_t N>
dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
Gradient<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator()(
    const dare::Data::CenterValueStencil<GridType, N>& s_c) const {
    dare::Data::FaceValueStencil<GridType, N> s_f;
    for (std::size_t n{0}; n < N; n++) {
        s_f.SetValue(Positions::WEST, n,
            (s_c.GetValue(Positions::CENTER, n) - s_c.GetValue(Positions::WEST, n)) * dn_r[0]);
        s_f.SetValue(Positions::EAST, n,
                     (s_c.GetValue(Positions::EAST, n) - s_c.GetValue(Positions::CENTER, n)) * dn_r[0]);
        if constexpr(Dim > 1) {
            s_f.SetValue(Positions::SOUTH, n,
                         (s_c.GetValue(Positions::CENTER, n) - s_c.GetValue(Positions::SOUTH, n)) * dn_r[1]);
            s_f.SetValue(Positions::EAST, n,
                         (s_c.GetValue(Positions::NORTH, n) - s_c.GetValue(Positions::CENTER, n)) * dn_r[1]);
        }
        if constexpr (Dim > 2) {
            s_f.SetValue(Positions::BOTTOM, n,
                         (s_c.GetValue(Positions::CENTER, n) - s_c.GetValue(Positions::BOTTOM, n)) * dn_r[2]);
            s_f.SetValue(Positions::TOP, n,
                         (s_c.GetValue(Positions::TOP, n) - s_c.GetValue(Positions::CENTER, n)) * dn_r[2]);
        }
    }
    return s_f;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
template <std::size_t N>
dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, 1>
Gradient<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator() (
    const dare::Data::GridVector<GridType, SC, N>& field, std::size_t n) const {
    dare::Data::FaceValueStencil<GridType, 1> s_f;
    const LO ordinal_local = grid->MapInternalToLocal(ordinal_internal);
    Index ind = grid->MapOrdinalToIndexLocal(ordinal_local);

    Index ind_nb(ind);
    ind_nb.i() = ind.i() - 1;
    s_f.SetValue(Positions::WEST, 1,
                 (field.At(ind, n) - field.At(ind_nb, n)) * dn_r[0]);
    ind_nb.i() = ind.i() + 1;
    s_f.SetValue(Positions::EAST, 1,
                 (field.At(ind_nb, n) - field.At(ind, n)) * dn_r[0]);
    if constexpr (Dim > 1) {
        ind_nb = ind;
        ind_nb.j() = ind.j() - 1;
        s_f.SetValue(Positions::SOUTH, 1,
                     (field.At(ind, n) - field.At(ind_nb, n)) * dn_r[1]);
        ind_nb.j() = ind.j() + 1;
        s_f.SetValue(Positions::NORTH, 1,
                     (field.At(ind_nb, n) - field.At(ind, n)) * dn_r[1]);
    }
    if constexpr (Dim > 2) {
        ind_nb = ind;
        ind_nb.k() = ind.k() - 1;
        s_f.SetValue(Positions::BOTTOM, 1,
                     (field.At(ind, n) - field.At(ind_nb, n)) * dn_r[2]);
        ind_nb.k() = ind.k() + 1;
        s_f.SetValue(Positions::TOP, 1,
                     (field.At(ind_nb, n) - field.At(ind, n)) * dn_r[2]);
    }
    return s_f;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
template <std::size_t N>
dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, 1>
Gradient<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator()(
    const dare::Data::CenterValueStencil<GridType, N>& s_c, std::size_t n) const {
    dare::Data::FaceValueStencil<GridType, 1> s_f;

    s_f.SetValue(Positions::WEST, 1,
                 (s_c.GetValue(Positions::CENTER, n) - s_c.GetValue(Positions::WEST, n)) * dn_r[0]);
    s_f.SetValue(Positions::EAST, 1,
                 (s_c.GetValue(Positions::EAST, n) - s_c.GetValue(Positions::CENTER, n)) * dn_r[0]);
    if constexpr (Dim > 1) {
        s_f.SetValue(Positions::SOUTH, 1,
                     (s_c.GetValue(Positions::CENTER, n) - s_c.GetValue(Positions::SOUTH, n)) * dn_r[1]);
        s_f.SetValue(Positions::EAST, 1,
                     (s_c.GetValue(Positions::NORTH, n) - s_c.GetValue(Positions::CENTER, n)) * dn_r[1]);
    }
    if constexpr (Dim > 2) {
        s_f.SetValue(Positions::BOTTOM, 1,
                     (s_c.GetValue(Positions::CENTER, n) - s_c.GetValue(Positions::BOTTOM, n)) * dn_r[2]);
        s_f.SetValue(Positions::TOP, 1,
                     (s_c.GetValue(Positions::TOP, n) - s_c.GetValue(Positions::CENTER, n)) * dn_r[2]);
    }
    return s_f;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
Divergence<dare::Grid::Cartesian<Dim, LO, GO, SC>>::Divergence(
    const GridRepresentation& grid, LO ordinal_internal)
    : A(grid.GetFaceArea()) {
    // the ordinal is here for the standard interface, nothing more
    // it's not really required when using the Cartesian grid instance
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
template <std::size_t N>
dare::utils::Vector<N, SC>
Divergence<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator()(
    const dare::Data::FaceValueStencil<GridType, N>& s) const {
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
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
template <std::size_t N>
dare::Data::CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
Divergence<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator()(
    const dare::Data::FaceMatrixStencil<GridType, N>& s) const {
    dare::Data::CenterMatrixStencil<GridType, N> s_c;
    for (std::size_t n{0}; n < N; n++) {
        // Divergence in X
        SC coef_c = s.GetValueCenter(Positions::EAST, n);
        SC coef_f = s.GetValueNeighbor(Positions::EAST, n);
        s_c.GetValue(Positions::EAST, n)   = A[0] * coef_f;
        s_c.GetValue(Positions::CENTER, n) = A[0] * coef_c;

        coef_c = s.GetValueCenter(Positions::WEST, n);
        coef_f = s.GetValueNeighbor(Positions::WEST, n);
        s_c.GetValue(Positions::WEST, n)    = -A[0] * coef_f;
        s_c.GetValue(Positions::CENTER, n) -=  A[0] * coef_c;

        // Divergence in Y
        if constexpr(Dim > 1) {
            SC coef_c = s.GetValueCenter(Positions::NORTH, n);
            SC coef_f = s.GetValueNeighbor(Positions::NORTH, n);
            s_c.GetValue(Positions::NORTH, n)   = A[1] * coef_f;
            s_c.GetValue(Positions::CENTER, n) += A[1] * coef_c;

            coef_c = s.GetValueCenter(Positions::SOUTH, n);
            coef_f = s.GetValueNeighbor(Positions::SOUTH, n);
            s_c.GetValue(Positions::SOUTH, n)   = -A[1] * coef_f;
            s_c.GetValue(Positions::CENTER, n) -=  A[1] * coef_c;
        }

        // Divergence in Z
        if constexpr (Dim > 1) {
            SC coef_c = s.GetValueCenter(Positions::TOP, n);
            SC coef_f = s.GetValueNeighbor(Positions::TOP, n);
            s_c.GetValue(Positions::TOP, n)     = A[2] * coef_f;
            s_c.GetValue(Positions::CENTER, n) += A[2] * coef_c;

            coef_c = s.GetValueCenter(Positions::BOTTOM, n);
            coef_f = s.GetValueNeighbor(Positions::BOTTOM, n);
            s_c.GetValue(Positions::BOTTOM, n)  = -A[2] * coef_f;
            s_c.GetValue(Positions::CENTER, n) -=  A[2] * coef_c;
        }
    }
}

}  // end namespace dare::Matrix