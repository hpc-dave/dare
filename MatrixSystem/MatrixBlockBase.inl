/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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
template <typename O, typename SC, std::size_t N>
MatrixBlockBase<O, SC, N>::MatrixBlockBase()
    : MatrixBlockBase<O, SC, N>(0, dare::utils::Vector<N, std::size_t>()) {
}

template <typename O, typename SC, std::size_t N>
MatrixBlockBase<O, SC, N>::MatrixBlockBase(const O& _node, const dare::utils::Vector<N, std::size_t>& size_hint)
    : initial_guess("initial guess", N), rhs("RHS", N), node(_node) {
    ProvideSizeHint(size_hint);
}

template <typename O, typename SC, std::size_t N>
MatrixBlockBase<O, SC, N>::MatrixBlockBase(const MatrixBlockBase<O, SC, N>& other)
    : node(other.node) {
    for (std::size_t n{0}; n < N; n++) {
        Kokkos::resize(ordinals[n], other.ordinals[n].size());
        Kokkos::resize(coefficients[n], other.coefficients[n].size());
        Kokkos::deep_copy(ordinals[n], other.ordinals[n]);
        Kokkos::deep_copy(coefficients[n], other.coefficients[n]);
    }
    Kokkos::resize(rhs, other.rhs.size());
    Kokkos::resize(initial_guess, other.initial_guess.size());
    Kokkos::deep_copy(rhs, other.rhs);
    Kokkos::deep_copy(initial_guess, other.initial_guess);
}

template <typename O, typename SC, std::size_t N>
MatrixBlockBase<O, SC, N>& MatrixBlockBase<O, SC, N>::operator=(const MatrixBlockBase<O, SC, N>& other) {
    if (this == &other)
        return *this;
    node = other.node;
    for (std::size_t n{0}; n < N; n++) {
        Kokkos::resize(ordinals[n], other.ordinals[n].size());
        Kokkos::resize(coefficients[n], other.coefficients[n].size());
        Kokkos::deep_copy(ordinals[n], other.ordinals[n]);
        Kokkos::deep_copy(coefficients[n], other.coefficients[n]);
    }
    Kokkos::resize(rhs, other.rhs.size());
    Kokkos::resize(initial_guess, other.initial_guess.size());
    Kokkos::deep_copy(rhs, other.rhs);
    Kokkos::deep_copy(initial_guess, other.initial_guess);
    return *this;
}

template <typename O, typename SC, std::size_t N>
MatrixBlockBase<O, SC, N>::~MatrixBlockBase() {}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::Initialize(O _node, const dare::utils::Vector<N, std::size_t>& size_hint) {
    node = _node;
    ProvideSizeHint(size_hint);
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::ProvideSizeHint(const dare::utils::Vector<N, std::size_t>& size_hint) {
    for (std::size_t n{0}; n < N; n++) {
        ordinals[n] = OrdinalArray("ordinals", size_hint[n]);
        coefficients[n] = ScalarArray("scalars", size_hint[n]);
        for (std::size_t i{0}; i < size_hint[n]; i++)
            ordinals[n][i] = -1;
    }
}

template <typename O, typename SC, std::size_t N>
O MatrixBlockBase<O, SC, N>::GetRow(std::size_t n) const {
    return node * N + n;
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetRhs(std::size_t n) const {
    return rhs[n];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetRhs(std::size_t n) {
    return rhs[n];
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetInitialGuess(std::size_t n) const {
    return initial_guess[n];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetInitialGuess(std::size_t n) {
    return initial_guess[n];
}

template <typename O, typename SC, std::size_t N>
std::size_t MatrixBlockBase<O, SC, N>::GetNumEntries(std::size_t n) const {
    return ordinals[n].size();
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetRhs(std::size_t n, SC value) {
    rhs[n] = value;
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetInitialGuess(std::size_t n, SC value) {
    initial_guess[n] = value;
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetRhs(const SC* rhs_array) {
    for (std::size_t n{0}; n < N; n++) {
        SetRhs(n, *(rhs_array + n));
    }
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetInitialGuess(const SC* x_array) {
    for (std::size_t n{0}; n < N; n++) {
        SetInitialGuess(n, *(x_array + n));
    }
}

template <typename O, typename SC, std::size_t N>
const typename MatrixBlockBase<O, SC, N>::OrdinalArray&
MatrixBlockBase<O, SC, N>::GetColumnOrdinals(std::size_t n) const {
    return ordinals[n];
}

template <typename O, typename SC, std::size_t N>
const typename MatrixBlockBase<O, SC, N>::ScalarArray& MatrixBlockBase<O, SC, N>::GetColumnValues(std::size_t n) const {
    return coefficients[n];
}

template <typename O, typename SC, std::size_t N>
O& MatrixBlockBase<O, SC, N>::GetOrdinalByPosition(std::size_t n, std::size_t pos) {
    return ordinals[n][pos];
}

template <typename O, typename SC, std::size_t N>
O MatrixBlockBase<O, SC, N>::GetOrdinalByPosition(std::size_t n, std::size_t pos) const {
    return ordinals[n][pos];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetCoefficientByPosition(std::size_t n, std::size_t pos) {
    return coefficients[n][pos];
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetCoefficientByPosition(std::size_t n, std::size_t pos) const {
    return coefficients[n][pos];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetCoefficientByOrdinal(std::size_t n, O ordinal) {
    O pos{-1};
    for (std::size_t i{0}; i < ordinals[n].size(); i++) {
        if (ordinals[n][i] == ordinal) {
            pos = i;
            break;
        }
    }
    if (pos == -1) {
        pos = ordinals[n].size();
        Kokkos::resize(coefficients[n], coefficients[n].size() + 1);
        Kokkos::resize(ordinals[n], ordinals[n].size() + 1);
        ordinals[n][pos] = ordinal;
        coefficients[n][pos] = 0;
    }

    return GetCoefficientByPosition(n, static_cast<std::size_t>(pos));
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetCoefficientByOrdinal(std::size_t n, O ordinal) const {
    std::size_t pos{0};
    for (std::size_t i{0}; i < ordinals[n].size(); i++) {
        if (ordinals[n][i] == ordinal)
            pos = i;
    }
    return GetCoefficientByPosition(n, pos);
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetCoefficient(std::size_t n, O id_col, SC value) {
    for (std::size_t count{0}; count < ordinals[n].size(); count++) {
        if (ordinals[n][count] == id_col || ordinals[n][count] == -1) {
            ordinals[n][count] = id_col;
            coefficients[n][count] = value;
            break;
        }
    }
}

template <typename O, typename SC, std::size_t N>
template <typename Array1, typename Array2>
void MatrixBlockBase<O, SC, N>::SetCoefficients(std::size_t n_row,
                                                std::size_t size,
                                                const Array1& id_col,
                                                const Array2& values) {
#ifndef DARE_NDEBUG
    if (size > ordinals.size())
        std::cerr << "SetCoefficients received a size value larger than the allocated storage "
                  << "(" << size << " > " << ordinals.size() << ")! "
                  << "Expect Segmentation fault!" << std::endl;
#endif

    for (std::size_t n{0}; n < size; n++) {
        ordinals[n_row][n] = id_col[n];
        coefficients[n_row][n] = values[n];
    }
}
}  // namespace dare::Matrix
