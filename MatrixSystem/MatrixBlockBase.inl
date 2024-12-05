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
        ordinals[n].resize(other.ordinals[n].h_view.size());
        coefficients[n].resize(other.coefficients[n].h_view.size());
        Kokkos::deep_copy(ordinals[n], other.ordinals[n]);
        Kokkos::deep_copy(coefficients[n], other.coefficients[n]);
    }
    rhs.resize(other.rhs.h_view.size());
    initial_guess.resize(other.initial_guess.h_view.size());
    Kokkos::deep_copy(rhs, other.rhs);
    Kokkos::deep_copy(initial_guess, other.initial_guess);
    Synchronize<ExecutionSpace>();
}

template <typename O, typename SC, std::size_t N>
MatrixBlockBase<O, SC, N>& MatrixBlockBase<O, SC, N>::operator=(const MatrixBlockBase<O, SC, N>& other) {
    if (this == &other)
        return *this;
    node = other.node;
    for (std::size_t n{0}; n < N; n++) {
        ordinals[n].resize(other.ordinals[n].h_view.size());
        coefficients[n].resize(other.coefficients[n].h_view.size());
        Kokkos::deep_copy(ordinals[n], other.ordinals[n]);
        Kokkos::deep_copy(coefficients[n], other.coefficients[n]);
    }
    rhs.resize(other.rhs.h_view.size());
    initial_guess.resize(other.initial_guess.h_view.size());
    Kokkos::deep_copy(rhs, other.rhs);
    Kokkos::deep_copy(initial_guess, other.initial_guess);
    Synchronize<ExecutionSpace>();
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
        ordinals[n] = DualViewOrdinalArrayType("ordinals", size_hint[n]);
        coefficients[n] = DualViewScalarArrayType("scalars", size_hint[n]);
        for (std::size_t i{0}; i < size_hint[n]; i++) {
            ordinals[n].h_view[i] = -1;
            ordinals[n].d_view[i] = -1;
        }
    }
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::Resize(std::size_t n, std::size_t size) {
    std::size_t old_size{ordinals[n].h_view.size()};
    ordinals[n].resize(size);
    coefficients[n].resize(size);
    for (std::size_t i{old_size}; i < size; i++) {
        ordinals[n].d_view[i] = -1;
        ordinals[n].h_view[i] = -1;
    }
}

template <typename O, typename SC, std::size_t N>
template <typename Space>
std::size_t MatrixBlockBase<O, SC, N>::GetSize(std::size_t n) const {
    if constexpr (std::is_same_v<Space, HostSpace>)
        return ordinals[n].h_view.size();
    else
        return ordinals[n].d_view.size();
}

template <typename O, typename SC, std::size_t N>
constexpr std::size_t MatrixBlockBase<O, SC, N>::GetNumComponents() const {
    return N;
}

template <typename O, typename SC, std::size_t N>
O MatrixBlockBase<O, SC, N>::GetNode() const {
    return node;
}

template <typename O, typename SC, std::size_t N>
O MatrixBlockBase<O, SC, N>::GetRow(std::size_t n) const {
    return node * N + n;
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetRhs(std::size_t n) const {
    return rhs.h_view[n];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetRhs(std::size_t n) {
    return rhs.h_view[n];
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetInitialGuess(std::size_t n) const {
    return initial_guess.h_view[n];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetInitialGuess(std::size_t n) {
    return initial_guess.h_view[n];
}

template <typename O, typename SC, std::size_t N>
std::size_t MatrixBlockBase<O, SC, N>::GetNumEntries(std::size_t n) const {
    std::size_t s{0};
    for (std::size_t i{0}; i < ordinals[n].h_view.size(); i++) {
        s += ordinals[n].h_view[i] > -1;
    }
    return s;
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetRhs(std::size_t n, SC value) {
    rhs.h_view[n] = value;
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetInitialGuess(std::size_t n, SC value) {
    initial_guess.h_view[n] = value;
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
    return ordinals[n].h_view;
}

template <typename O, typename SC, std::size_t N>
const typename MatrixBlockBase<O, SC, N>::OrdinalArrayDevice&
MatrixBlockBase<O, SC, N>::GetColumnOrdinalsDevice(std::size_t n) const {
    return ordinals[n].d_view;
}

template <typename O, typename SC, std::size_t N>
const typename MatrixBlockBase<O, SC, N>::ScalarArray&
MatrixBlockBase<O, SC, N>::GetColumnValues(std::size_t n) const {
    return coefficients[n].h_view;
}

template <typename O, typename SC, std::size_t N>
const typename MatrixBlockBase<O, SC, N>::ScalarArrayDevice&
MatrixBlockBase<O, SC, N>::GetColumnValuesDevice(std::size_t n) const {
    return coefficients[n].d_view;
}

template <typename O, typename SC, std::size_t N>
O& MatrixBlockBase<O, SC, N>::GetOrdinalByPosition(std::size_t n, std::size_t pos) {
    return ordinals[n].h_view[pos];
}

template <typename O, typename SC, std::size_t N>
O MatrixBlockBase<O, SC, N>::GetOrdinalByPosition(std::size_t n, std::size_t pos) const {
    return ordinals[n].h_view[pos];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetCoefficientByPosition(std::size_t n, std::size_t pos) {
    return coefficients[n].h_view[pos];
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetCoefficientByPosition(std::size_t n, std::size_t pos) const {
    return coefficients[n].h_view[pos];
}

template <typename O, typename SC, std::size_t N>
SC& MatrixBlockBase<O, SC, N>::GetCoefficientByOrdinal(std::size_t n, O ordinal) {
    O pos{-1};
    for (std::size_t i{0}; i < ordinals[n].h_view.size(); i++) {
        if (ordinals[n].h_view[i] == ordinal) {
            pos = i;
            break;
        }
    }
    if (pos == -1) {
        pos = ordinals[n].h_view.size();
        Kokkos::resize(coefficients[n].h_view, coefficients[n].h_view.size() + 1);
        Kokkos::resize(ordinals[n].h_view, ordinals[n].h_view.size() + 1);
        ordinals[n].h_view[pos] = ordinal;
        coefficients[n].h_view[pos] = 0;
    }

    return GetCoefficientByPosition(n, static_cast<std::size_t>(pos));
}

template <typename O, typename SC, std::size_t N>
SC MatrixBlockBase<O, SC, N>::GetCoefficientByOrdinal(std::size_t n, O ordinal) const {
    std::size_t pos{0};
    for (std::size_t i{0}; i < ordinals[n].h_view.size(); i++) {
        if (ordinals[n].h_view[i] == ordinal)
            pos = i;
    }
    return GetCoefficientByPosition(n, pos);
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::SetCoefficient(std::size_t n, O id_col, SC value) {
#ifndef DARE_NDEBUG
    bool substituted{false};
#endif
    for (std::size_t count{0}; count < ordinals[n].h_view.size(); count++) {
        if (ordinals[n].h_view[count] == id_col || ordinals[n].h_view[count] == -1) {
            ordinals[n].h_view[count] = id_col;
            coefficients[n].h_view[count] = value;
#ifndef DARE_NDEBUG
            substituted = true;
#endif
            break;
        }
    }
#ifndef DARE_NDEBUG
    if (!substituted) {
        std::cerr << "Coefficient " << id_col << " was not found in the array "
           << "and no additional space is available! Resize the array appropriately beforehand!" << std::endl;
    }
#endif
}

template <typename O, typename SC, std::size_t N>
template <typename Array1, typename Array2>
void MatrixBlockBase<O, SC, N>::SetCoefficients(std::size_t n_row,
                                                std::size_t size,
                                                const Array1& id_col,
                                                const Array2& values,
                                                bool allow_resize) {
    if (allow_resize) {
        this->Resize(n_row, this->GetNumEntries(n_row) + size);
    }
    bool is_empty{ordinals[n_row].h_view[0] < 0};
#ifndef DARE_NDEBUG
    if (size > ordinals[n_row].h_view.size())
        std::cerr << "In " << __func__ << ": received a size value larger than the allocated storage "
                  << "(" << size << " > " << ordinals[n_row].h_view.size() << ")! "
                  << "Expect Segmentation fault!" << std::endl;
    for (std::size_t n{0}; n < n_row; n++) {
        if (id_col[n] < 0) {
            std::cerr << "In " << __func__ << ": received an ordinal < 0 "
                      << "(" << id_col[n] << ") at row " << n_row << "! "
                      << "Expect error thrown by Trilinos!" << std::endl;
        }
        if (std::isnan(values[n]) || !std::isfinite(values[n])) {
            std::cerr << "In " << __func__ << ": received an unreal coefficient "
                      << "(" << values[n] << ") at row " << n_row << " and column " << id_col[n] << "! "
                      << "Expect nonsense solution!" << std::endl;
        }
    }
#endif

    if (is_empty) {
        for (std::size_t n{0}; n < size; n++) {
            ordinals[n_row].h_view[n] = id_col[n];
            coefficients[n_row].h_view[n] = values[n];
        }
    } else {
        for (std::size_t n{0}; n < size; n++) {
            SetCoefficient(n_row, id_col[n], values[n]);
        }
    }

    if (allow_resize) {
        this->Resize(n_row, this->GetNumEntries(n_row));
    }
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::RemoveCoefficientByPosition(std::size_t n, std::size_t pos) {
    const std::size_t range_new{ordinals[n].h_view.size() - 1};
    for (; pos < range_new; pos++) {
        ordinals[n].h_view[pos] = ordinals[n].h_view[pos + 1];
        coefficients[n].h_view[pos] = coefficients[n].h_view[pos + 1];
    }
    ordinals[n].resize(range_new);
    coefficients[n].resize(range_new);
}

template <typename O, typename SC, std::size_t N>
template <typename Array>
void MatrixBlockBase<O, SC, N>::RemoveCoefficientsByPositions(std::size_t n,
                                                              const Array& positions,
                                                              std::size_t num_entries) {
    OrdinalArray ordinals_new("ordinals", ordinals[n].h_view.size() - num_entries);
    ScalarArray coefficients_new("coefficients", coefficients[n].h_view.size() - num_entries);
    std::size_t q{0};
    for (std::size_t p{0}; p < ordinals[n].h_view.size(); p++) {
        bool skip{false};
        for (std::size_t i{0}; i < num_entries; i++) {
            skip |= positions[i] == p;
        }
        if (skip)
            continue;

        ordinals_new[q] = ordinals[n].h_view[p];
        coefficients_new[q] = coefficients[n].h_view[p];
        ++q;
    }
    ordinals[n].resize(ordinals[n].h_view.size() - num_entries);
    coefficients[n].resize(coefficients[n].h_view.size() - num_entries);
    Kokkos::deep_copy(ordinals[n].h_view, ordinals_new);
    Kokkos::deep_copy(coefficients[n].h_view, coefficients_new);
}

template <typename O, typename SC, std::size_t N>
void MatrixBlockBase<O, SC, N>::RemoveCoefficientByOrdinal(std::size_t n, O ordinal) {
    std::size_t pos{0};
    for (; pos < ordinals[n].h_view.size(); pos++) {
        if (ordinals[n].h_view[pos] == ordinal)
            break;
    }
    RemoveCoefficientByPosition(n, pos);
}

template <typename O, typename SC, std::size_t N>
template <typename Array>
void MatrixBlockBase<O, SC, N>::RemoveCoefficientsByOrdinals(std::size_t n,
                                                             const Array& col_ids,
                                                             std::size_t num_entries) {
    std::vector<std::size_t> pos(num_entries);
    for (std::size_t p{0}; p < num_entries; p++) {
        for (std::size_t q{0}; q < ordinals[n].h_view.size(); q++) {
            if (col_ids[p] == ordinals[n].h_view[q]) {
                pos[p] = q;
                break;
            }
        }
    }
    RemoveCoefficientsByPositions(n, pos, num_entries);
}

template <typename O, typename SC, std::size_t N>
template <typename TargetSpace>
void MatrixBlockBase<O, SC, N>::Synchronize() {
    if constexpr (std::is_same_v<TargetSpace, HostSpace>) {
        for (std::size_t n{0}; n < N; n++) {
            ordinals[n].template modify<ExecutionSpace>();
            coefficients[n].template modify<ExecutionSpace>();
            Kokkos::resize(ordinals[n].h_view, ordinals[n].d_view.size());
            Kokkos::resize(coefficients[n].h_view, coefficients[n].d_view.size());
        }
        rhs.template modify<ExecutionSpace>();
        initial_guess.template modify<ExecutionSpace>();
        Kokkos::resize(rhs.h_view, rhs.d_view.size());
        Kokkos::resize(initial_guess.h_view, initial_guess.d_view.size());
        for (std::size_t n{0}; n < N; n++) {
            ordinals[n].template sync<HostSpace>();
            coefficients[n].template sync<HostSpace>();
        }
        rhs.template sync<HostSpace>();
        initial_guess.template sync<HostSpace>();
    } else if constexpr (std::is_same_v<TargetSpace, ExecutionSpace>) {
        for (std::size_t n{0}; n < N; n++) {
            ordinals[n].template modify<HostSpace>();
            coefficients[n].template modify<HostSpace>();
            Kokkos::resize(ordinals[n].d_view, ordinals[n].h_view.size());
            Kokkos::resize(coefficients[n].d_view, coefficients[n].h_view.size());
        }
        rhs.template modify<HostSpace>();
        initial_guess.template modify<HostSpace>();
        Kokkos::resize(rhs.d_view, rhs.h_view.size());
        Kokkos::resize(initial_guess.d_view, initial_guess.h_view.size());
        for (std::size_t n{0}; n < N; n++) {
            ordinals[n].template sync<ExecutionSpace>();
            coefficients[n].template sync<ExecutionSpace>();
        }
        rhs.template sync<ExecutionSpace>();
        initial_guess.template sync<ExecutionSpace>();
    }
}
}  // namespace dare::Matrix
