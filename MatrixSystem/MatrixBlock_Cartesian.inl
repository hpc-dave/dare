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

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::MatrixBlock()
    : MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>(nullptr, 0, dare::utils::Vector<N, std::size_t>()) {
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::MatrixBlock(const GridRepresentation* _g_rep,
                                         O _node)
    : MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>(_g_rep, _node, dare::utils::Vector<N, std::size_t>()) {
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::MatrixBlock(const GridRepresentation* _g_rep,
                                                               O _node,
                                                               const dare::utils::Vector<N, std::size_t>& size_hint)
    : MatrixBlockBase<O, SC, N>(_node, size_hint), g_rep(_g_rep),
      neighbors("neighbors", N, N * 2 + 1), neighbor_set(N) {
    static_assert(std::is_same_v<O, LocalOrdinalType>
               || std::is_same_v<O, GlobalOrdinalType>,
                  "The ordinal type needs to be either a local or global ordinal type!");
    for (std::size_t n{0}; n < N; n++) {
        // neighbor_set.h_view(n) = 0;
        // neighbor_set.d_view(n) = 0;
        neighbor_set[n] = 0;
        neighbor_set[n] = 0;
    }
    if (g_rep != nullptr) {
        if constexpr (IsGlobal()) {
            ind_internal = g_rep->MapOrdinalToIndexGlobalInternal(_node);
            ind_full = g_rep->MapInternalToLocal(ind_internal);
        } else {
            ind_internal = g_rep->MapOrdinalToIndexLocalInternal(_node);
            ind_full = g_rep->MapInternalToLocal(ind_internal);
        }
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::MatrixBlock(const SelfType& other)
    : MatrixBlockBase<O, SC, N>(other), g_rep(other.g_rep) {
    static_assert(std::is_same_v<O, LocalOrdinalType>
               || std::is_same_v<O, GlobalOrdinalType>,
                  "The ordinal type needs to be either a local or global ordinal type!");
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::~MatrixBlock() {}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator=(const SelfType& other) {
    if (&other == this)
        return *this;
    MatrixBlockBase<O, SC, N>::operator=(other);
    g_rep = other.g_rep;
    return *this;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
void MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Initialize(const GridRepresentation* _g_rep, O _node) {
    Initialize(_g_rep, _node, dare::utils::Vector<N, std::size_t>());
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
void MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Initialize(const GridRepresentation* _g_rep,
                                             O _node,
                                             const dare::utils::Vector<N, std::size_t>& size_hint) {
    MatrixBlockBase<O, SC, N>::Initialize(_node, size_hint);
    g_rep = _g_rep;
    if constexpr (IsGlobal()) {
        ind_internal = g_rep->MapOrdinalToIndexGlobalInternal(_node);
        ind_full = g_rep->MapInternalToLocal(ind_internal);
    } else {
        ind_internal = g_rep->MapOrdinalToIndexLocalInternal(_node);
        ind_full = g_rep->MapInternalToLocal(ind_internal);
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
const typename MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GridRepresentation*
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetRepresentation() const {
    return g_rep;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <CartesianNeighbor CNB>
bool MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::IsSet(std::size_t n) const {
    if constexpr (IsSame<CartesianNeighbor::CENTER, CNB>())
        return GetNeighborBitSet<HostSpace>()[n] & static_cast<char>(CartesianNeighborBitSet::CENTER);
    else if constexpr (IsSame<CartesianNeighbor::WEST, CNB>())
        return GetNeighborBitSet<HostSpace>()[n] & static_cast<char>(CartesianNeighborBitSet::WEST);
    else if constexpr (IsSame<CartesianNeighbor::EAST, CNB>())
        return GetNeighborBitSet<HostSpace>()[n] & static_cast<char>(CartesianNeighborBitSet::EAST);
    else if constexpr (IsSame<CartesianNeighbor::SOUTH, CNB>() && (Dim > 1))
        return GetNeighborBitSet<HostSpace>()[n] & static_cast<char>(CartesianNeighborBitSet::SOUTH);
    else if constexpr (IsSame<CartesianNeighbor::NORTH, CNB>() && (Dim > 1))
        return GetNeighborBitSet<HostSpace>()[n] & static_cast<char>(CartesianNeighborBitSet::NORTH);
    else if constexpr (IsSame<CartesianNeighbor::BOTTOM, CNB>() && (Dim > 2))
        return GetNeighborBitSet<HostSpace>()[n] & static_cast<char>(CartesianNeighborBitSet::BOTTOM);
    else if constexpr (IsSame<CartesianNeighbor::TOP, CNB>() && (Dim > 2))
        return GetNeighborBitSet<HostSpace>()[n] & static_cast<char>(CartesianNeighborBitSet::TOP);
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
void MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Finalize() {
    using CN = CartesianNeighbor;
    O node = this->GetNode();
    Index ind = ind_internal;
    // if constexpr (this->IsGlobal()) {
    //     ind = g_rep->MapOrdinalToIndexGlobal(node);
    //     ind = g_rep->MapGlobalToInternal(ind);
    // } else {
    //     ind = g_rep->MapOrdinalToIndexLocal(node);
    //     ind = g_rep->MapLocalToInternal(ind);
    // }
    for (std::size_t n{0}; n < N; n++) {
        char num_entries{0};
        for (std::size_t p{0}; p < 8; p++) {
            num_entries += (GetNeighborBitSet<HostSpace>()[n] & (1 << p)) != 0;
        }
        std::vector<O> ordinals(num_entries);
        std::vector<SC> coef(num_entries);
        std::size_t pos{0};

        if (IsSet<CN::WEST>(n)) {
            O col{0};
            ind.i()--;
            if constexpr (this->IsGlobal()) {
                col = g_rep->MapIndexToOrdinalGlobalInternal(ind) * N + n;
            } else {
                col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
            }
            ind.i()++;
            ordinals[pos] = col;
            coef[pos] = Get<CN::WEST>(n);
            pos++;
        }
        if constexpr (Dim > 1) {
            if (IsSet<CN::SOUTH>(n)) {
                O col{0};
                ind.j()--;
                if constexpr (this->IsGlobal()) {
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind) * N + n;
                } else {
                    col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
                }
                ind.j()++;
                ordinals[pos] = col;
                coef[pos] = Get<CN::SOUTH>(n);
                pos++;
            }
        }
        if constexpr (Dim > 2) {
            if (IsSet<CN::BOTTOM>(n)) {
                O col{0};
                ind.k()--;
                if constexpr (this->IsGlobal()) {
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind) * N + n;
                } else {
                    col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
                }
                ind.k()++;
                ordinals[pos] = col;
                coef[pos] = Get<CN::BOTTOM>(n);
                pos++;
            }
        }
        if (IsSet<CN::CENTER>(n)) {
            O col{0};
            if constexpr (this->IsGlobal()) {
                col = g_rep->MapIndexToOrdinalGlobalInternal(ind) * N + n;
            } else {
                col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
            }
            ordinals[pos] = col;
            coef[pos] = Get<CN::CENTER>(n);
            pos++;
        }
        if constexpr (Dim > 2) {
            if (IsSet<CN::TOP>(n)) {
                O col{0};
                ind.k()++;
                if constexpr (this->IsGlobal()) {
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind) * N + n;
                } else {
                    col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
                }
                ind.k()--;
                ordinals[pos] = col;
                coef[pos] = Get<CN::TOP>(n);
                pos++;
            }
        }
        if constexpr (Dim > 1) {
            if (IsSet<CN::NORTH>(n)) {
                O col{0};
                ind.j()++;
                if constexpr (this->IsGlobal()) {
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind) * N + n;
                } else {
                    col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
                }
                ind.j()--;
                ordinals[pos] = col;
                coef[pos] = Get<CN::NORTH>(n);
                pos++;
            }
        }
        if (IsSet<CN::EAST>(n)) {
            O col{0};
            ind.i()++;
            if constexpr (this->IsGlobal()) {
                col = g_rep->MapIndexToOrdinalGlobalInternal(ind) * N + n;
            } else {
                col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
            }
            ind.i()--;
            ordinals[pos] = col;
            coef[pos] = Get<CN::EAST>(n);
            pos++;
        }
        this->SetCoefficients(n, num_entries, ordinals, coef);
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <CartesianNeighbor CNB>
SC& MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n) {
    if constexpr (IsSame<CartesianNeighbor::CENTER, CNB>()) {
        GetNeighborBitSet<HostSpace>()[n] |= static_cast<char>(CartesianNeighborBitSet::CENTER);
        return GetNeighbors<HostSpace>()(n, static_cast<char>(CartesianNeighbor::CENTER));
    } else if constexpr (IsSame<CartesianNeighbor::WEST, CNB>()) {
        GetNeighborBitSet<HostSpace>()[n] |= static_cast<char>(CartesianNeighborBitSet::WEST);
        return GetNeighbors<HostSpace>()(n, static_cast<char>(CartesianNeighbor::WEST));
    } else if constexpr (IsSame<CartesianNeighbor::EAST, CNB>()) {
        GetNeighborBitSet<HostSpace>()[n] |= static_cast<char>(CartesianNeighborBitSet::EAST);
        return GetNeighbors<HostSpace>()(n, static_cast<char>(CartesianNeighbor::EAST));
    } else if constexpr (IsSame<CartesianNeighbor::SOUTH, CNB>() && (Dim > 1)) {
        GetNeighborBitSet<HostSpace>()[n] |= static_cast<char>(CartesianNeighborBitSet::SOUTH);
        return GetNeighbors<HostSpace>()(n, static_cast<char>(CartesianNeighbor::SOUTH));
    } else if constexpr (IsSame<CartesianNeighbor::NORTH, CNB>() && (Dim > 1)) {
        GetNeighborBitSet<HostSpace>()[n] |= static_cast<char>(CartesianNeighborBitSet::NORTH);
        return GetNeighbors<HostSpace>()(n, static_cast<char>(CartesianNeighbor::NORTH));
    } else if constexpr (IsSame<CartesianNeighbor::BOTTOM, CNB>() && (Dim > 2)) {
        GetNeighborBitSet<HostSpace>()[n] |= static_cast<char>(CartesianNeighborBitSet::BOTTOM);
        return GetNeighbors<HostSpace>()(n, static_cast<char>(CartesianNeighbor::BOTTOM));
    } else if constexpr (IsSame<CartesianNeighbor::TOP, CNB>() && (Dim > 2)) {
        GetNeighborBitSet<HostSpace>()[n] |= static_cast<char>(CartesianNeighborBitSet::TOP);
        return GetNeighbors<HostSpace>()(n, static_cast<char>(CartesianNeighbor::TOP));
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
SC& MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n, CartesianNeighbor cnb) {
    switch (cnb) {
    case CartesianNeighbor::CENTER:
        return Get<CartesianNeighbor::CENTER>();
    case CartesianNeighbor::WEST:
        return Get<CartesianNeighbor::WEST>();
    case CartesianNeighbor::EAST:
        return Get<CartesianNeighbor::EAST>();
    case CartesianNeighbor::SOUTH:
        return Get<CartesianNeighbor::SOUTH>();
    case CartesianNeighbor::NORTH:
        return Get<CartesianNeighbor::NORTH>();
    case CartesianNeighbor::BOTTOM:
        return Get<CartesianNeighbor::BOTTOM>();
    case CartesianNeighbor::TOP:
        return Get<CartesianNeighbor::TOP>();
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <CartesianNeighbor CNB>
SC MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n) const {
    if constexpr (IsSame<CartesianNeighbor::CENTER, CNB>()) {
        return neighbors(n, CartesianNeighbor::CENTER);
    } else if constexpr (IsSame<CartesianNeighbor::WEST, CNB>()) {
        return neighbors(n, CartesianNeighbor::WEST);
    } else if constexpr (IsSame<CartesianNeighbor::EAST, CNB>()) {
        return neighbors(n, CartesianNeighbor::EAST);
    } else if constexpr (IsSame<CartesianNeighbor::SOUTH, CNB>() && (Dim > 1)) {
        return neighbors(n, CartesianNeighbor::SOUTH);
    } else if constexpr (IsSame<CartesianNeighbor::NORTH, CNB>() && (Dim > 1)) {
        return neighbors(n, CartesianNeighbor::NORTH);
    } else if constexpr (IsSame<CartesianNeighbor::BOTTOM, CNB>() && (Dim > 2)) {
        return neighbors(n, CartesianNeighbor::BOTTOM);
    } else if constexpr (IsSame<CartesianNeighbor::TOP, CNB>() && (Dim > 2)) {
        return neighbors(n, CartesianNeighbor::TOP);
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
SC MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n, CartesianNeighbor cnb) const {
    switch (cnb) {
    case CartesianNeighbor::CENTER:
        return Get<CartesianNeighbor::CENTER>();
    case CartesianNeighbor::WEST:
        return Get<CartesianNeighbor::WEST>();
    case CartesianNeighbor::EAST:
        return Get<CartesianNeighbor::EAST>();
    case CartesianNeighbor::SOUTH:
        return Get<CartesianNeighbor::SOUTH>();
    case CartesianNeighbor::NORTH:
        return Get<CartesianNeighbor::NORTH>();
    case CartesianNeighbor::BOTTOM:
        return Get<CartesianNeighbor::BOTTOM>();
    case CartesianNeighbor::TOP:
        return Get<CartesianNeighbor::TOP>();
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
constexpr bool MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::IsGlobal() {
    return std::is_same_v<O, GlobalOrdinalType>;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
typename MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Index
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetIndex() const {
    return ind_full;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
typename MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Index
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetIndexInternal() const {
    return ind_internal;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
bool MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::IsStencilLocal() const {
    if constexpr (std::is_same_v<LocalOrdinalType, O>) {
        return true;
    } else {
        bool is_local{true};
        for (std::size_t n{0}; n < N; n++) {
            for (std::size_t i{0}; i < this->GetNumEntries(n); i++) {
                is_local &= g_rep->IsLocalInternal(this->GetOrdinalByPosition(n, i) / N);
            }
        }
        return is_local;
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <typename Space>
typename Kokkos::DualView<SC**>::t_host&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighbors() {
        return neighbors.h_view;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <typename Space>
const typename Kokkos::DualView<SC**>::t_host&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighbors() const {
        return neighbors.h_view;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <typename Space>
// typename Kokkos::DualView<char*>::t_host&
std::vector<char>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighborBitSet() {
    // return neighbor_set.h_view;
    return neighbor_set;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <typename Space>
// const typename Kokkos::DualView<char*>::t_host&
const std::vector<char>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighborBitSet() const {
    // return neighbor_set.h_view;
    return neighbor_set;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator=(
    const dare::Data::CenterMatrixStencil<GridType, N>& s) {
    for (std::size_t n{0}; n < N; n++) {
        this->Get<CartesianNeighbor::CENTER>() = s.GetValue(CartesianNeighbor::CENTER, n);
        this->Get<CartesianNeighbor::WEST>() = s.GetValue(CartesianNeighbor::WEST, n);
        this->Get<CartesianNeighbor::EAST>() = s.GetValue(CartesianNeighbor::EAST, n);
        if constexpr (Dim > 1) {
            this->Get<CartesianNeighbor::SOUTH>() = s.GetValue(CartesianNeighbor::SOUTH, n);
            this->Get<CartesianNeighbor::NORTH>() = s.GetValue(CartesianNeighbor::NORTH, n);
        }
        if constexpr (Dim > 2) {
            this->Get<CartesianNeighbor::BOTTOM>() = s.GetValue(CartesianNeighbor::BOTTOM, n);
            this->Get<CartesianNeighbor::TOP>() = s.GetValue(CartesianNeighbor::TOP, n);
        }
    }
    return *this;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator+=(
    const dare::Data::CenterMatrixStencil<GridType, N>& s) {
    for (std::size_t n{0}; n < N; n++) {
        this->Get<CartesianNeighbor::CENTER>() += s.GetValue(CartesianNeighbor::CENTER, n);
        this->Get<CartesianNeighbor::WEST>() += s.GetValue(CartesianNeighbor::WEST, n);
        this->Get<CartesianNeighbor::EAST>() += s.GetValue(CartesianNeighbor::EAST, n);
        if constexpr (Dim > 1) {
            this->Get<CartesianNeighbor::SOUTH>() += s.GetValue(CartesianNeighbor::SOUTH, n);
            this->Get<CartesianNeighbor::NORTH>() += s.GetValue(CartesianNeighbor::NORTH, n);
        }
        if constexpr (Dim > 2) {
            this->Get<CartesianNeighbor::BOTTOM>() += s.GetValue(CartesianNeighbor::BOTTOM, n);
            this->Get<CartesianNeighbor::TOP>() += s.GetValue(CartesianNeighbor::TOP, n);
        }
    }
    return *this;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator-=(
    const dare::Data::CenterMatrixStencil<GridType, N>& s) {
    return (*this) += (-1. * s);
}

}  // end namespace dare::Matrix
