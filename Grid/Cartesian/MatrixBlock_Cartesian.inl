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
    : MatrixBlockBase<O, SC, N>(_node, size_hint), g_rep(_g_rep) {
    static_assert(std::is_same_v<O, LocalOrdinalType>
               || std::is_same_v<O, GlobalOrdinalType>,
                  "The ordinal type needs to be either a local or global ordinal type!");
    neighbor_set.SetAllValues(0);
    for (std::size_t n{0}; n < N; ++n) {
        neighbors[n].resize(Dim * 2 + 1);
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
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator=(SelfType other) {
    std::swap(*this, other);
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
        return GetNeighborBitSet()[n] & static_cast<char>(CartesianNeighborBitSet::CENTER);
    else if constexpr (IsSame<CartesianNeighbor::WEST, CNB>())
        return GetNeighborBitSet()[n] & static_cast<char>(CartesianNeighborBitSet::WEST);
    else if constexpr (IsSame<CartesianNeighbor::EAST, CNB>())
        return GetNeighborBitSet()[n] & static_cast<char>(CartesianNeighborBitSet::EAST);
    else if constexpr (IsSame<CartesianNeighbor::SOUTH, CNB>() && (Dim > 1))
        return GetNeighborBitSet()[n] & static_cast<char>(CartesianNeighborBitSet::SOUTH);
    else if constexpr (IsSame<CartesianNeighbor::NORTH, CNB>() && (Dim > 1))
        return GetNeighborBitSet()[n] & static_cast<char>(CartesianNeighborBitSet::NORTH);
    else if constexpr (IsSame<CartesianNeighbor::BOTTOM, CNB>() && (Dim > 2))
        return GetNeighborBitSet()[n] & static_cast<char>(CartesianNeighborBitSet::BOTTOM);
    else if constexpr (IsSame<CartesianNeighbor::TOP, CNB>() && (Dim > 2))
        return GetNeighborBitSet()[n] & static_cast<char>(CartesianNeighborBitSet::TOP);

    ERROR << "The specified cartesian neighbor (" << std::to_string(ToNum(CNB)) << ") is out of range!";
    return false;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
void MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Finalize() {
    using CN = CartesianNeighbor;
    const bool map_periodic{true};
    Index ind = ind_internal;
    for (std::size_t n{0}; n < N; n++) {
        std::size_t num_entries{0};
        for (std::size_t p{0}; p < 8; p++) {
            num_entries += (GetNeighborBitSet()[n] & (1 << p)) != 0;
        }
        std::vector<O> ordinals(num_entries);
        std::vector<SC> coef(num_entries);
        std::size_t pos{0};

        if (IsSet<CN::WEST>(n)) {
            O col{0};
            ind.i()--;
            if constexpr (this->IsGlobal()) {
                col = g_rep->MapIndexToOrdinalGlobalInternal(ind, map_periodic) * N + n;
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
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind, map_periodic) * N + n;
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
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind, map_periodic) * N + n;
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
                col = g_rep->MapIndexToOrdinalGlobalInternal(ind, map_periodic) * N + n;
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
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind, map_periodic) * N + n;
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
                    col = g_rep->MapIndexToOrdinalGlobalInternal(ind, map_periodic) * N + n;
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
                col = g_rep->MapIndexToOrdinalGlobalInternal(ind, map_periodic) * N + n;
            } else {
                col = g_rep->MapIndexToOrdinalLocalInternal(ind) * N + n;
            }
            ind.i()--;
            ordinals[pos] = col;
            coef[pos] = Get<CN::EAST>(n);
            pos++;
        }
        this->SetCoefficients(n, num_entries, ordinals, coef, true);
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <CartesianNeighbor CNB>
SC& MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n) {
    if constexpr (IsSame<CartesianNeighbor::CENTER, CNB>()) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::CENTER);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::CENTER)];
    } else if constexpr (IsSame<CartesianNeighbor::WEST, CNB>()) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::WEST);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::WEST)];
    } else if constexpr (IsSame<CartesianNeighbor::EAST, CNB>()) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::EAST);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::EAST)];
    } else if constexpr (IsSame<CartesianNeighbor::SOUTH, CNB>() && (Dim > 1)) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::SOUTH);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::SOUTH)];
    } else if constexpr (IsSame<CartesianNeighbor::NORTH, CNB>() && (Dim > 1)) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::NORTH);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::NORTH)];
    } else if constexpr (IsSame<CartesianNeighbor::BOTTOM, CNB>() && (Dim > 2)) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::BOTTOM);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::BOTTOM)];
    } else if constexpr (IsSame<CartesianNeighbor::TOP, CNB>() && (Dim > 2)) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::TOP);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::TOP)];
    } else if constexpr (IsSame<CartesianNeighbor::FOURD_LOW, CNB>() && (Dim > 3)) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::FOURD_LOW);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::FOURD_LOW)];
    } else if constexpr (IsSame<CartesianNeighbor::FOURD_UP, CNB>() && (Dim > 3)) {
        GetNeighborBitSet()[n] |= static_cast<char>(CartesianNeighborBitSet::FOURD_UP);
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::FOURD_UP)];
    } else {
        ERROR << "The neighbor indicator is out of bounds!" << ERROR_CLOSE;
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::CENTER)];
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
SC& MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n, CartesianNeighbor cnb) {
    switch (cnb) {
    case CartesianNeighbor::CENTER:
        return Get<CartesianNeighbor::CENTER>(n);
        break;
    case CartesianNeighbor::WEST:
        return Get<CartesianNeighbor::WEST>(n);
        break;
    case CartesianNeighbor::EAST:
        return Get<CartesianNeighbor::EAST>(n);
        break;
    case CartesianNeighbor::SOUTH:
        return Get<CartesianNeighbor::SOUTH>(n);
        break;
    case CartesianNeighbor::NORTH:
        return Get<CartesianNeighbor::NORTH>(n);
        break;
    case CartesianNeighbor::BOTTOM:
        return Get<CartesianNeighbor::BOTTOM>(n);
        break;
    case CartesianNeighbor::TOP:
        return Get<CartesianNeighbor::TOP>(n);
        break;
    case CartesianNeighbor::FOURD_LOW:
        return Get<CartesianNeighbor::FOURD_LOW>(n);
        break;
    case CartesianNeighbor::FOURD_UP:
        return Get<CartesianNeighbor::FOURD_UP>(n);
        break;
    }
    ERROR << "something horrible happened here: " << ERROR_CLOSE;
    return Get<CartesianNeighbor::CENTER>(n);
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <CartesianNeighbor CNB>
void MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Remove(std::size_t n) {
    if constexpr (IsSame<CartesianNeighbor::CENTER, CNB>()) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::CENTER);
    } else if constexpr (IsSame<CartesianNeighbor::WEST, CNB>()) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::WEST);
    } else if constexpr (IsSame<CartesianNeighbor::EAST, CNB>()) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::EAST);
    } else if constexpr (IsSame<CartesianNeighbor::SOUTH, CNB>() && (Dim > 1)) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::SOUTH);
    } else if constexpr (IsSame<CartesianNeighbor::NORTH, CNB>() && (Dim > 1)) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::NORTH);
    } else if constexpr (IsSame<CartesianNeighbor::BOTTOM, CNB>() && (Dim > 2)) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::BOTTOM);
    } else if constexpr (IsSame<CartesianNeighbor::TOP, CNB>() && (Dim > 2)) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::TOP);
    } else if constexpr (IsSame<CartesianNeighbor::FOURD_LOW, CNB>() && (Dim > 3)) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::FOURD_LOW);
    } else if constexpr (IsSame<CartesianNeighbor::FOURD_UP, CNB>() && (Dim > 3)) {
        GetNeighborBitSet()[n] &= ~static_cast<char>(CartesianNeighborBitSet::FOURD_UP);
    } else {
        ERROR << "The neighbor indicator is out of bounds!" << ERROR_CLOSE;
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
void MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Remove(std::size_t n, CartesianNeighbor cnb) {
    switch (cnb) {
    case CartesianNeighbor::CENTER:
        Remove<CartesianNeighbor::CENTER>(n);
        break;
    case CartesianNeighbor::WEST:
        Remove<CartesianNeighbor::WEST>(n);
        break;
    case CartesianNeighbor::EAST:
        Remove<CartesianNeighbor::EAST>(n);
        break;
    case CartesianNeighbor::SOUTH:
        Remove<CartesianNeighbor::SOUTH>(n);
        break;
    case CartesianNeighbor::NORTH:
        Remove<CartesianNeighbor::NORTH>(n);
        break;
    case CartesianNeighbor::BOTTOM:
        Remove<CartesianNeighbor::BOTTOM>(n);
        break;
    case CartesianNeighbor::TOP:
        Remove<CartesianNeighbor::TOP>(n);
        break;
    case CartesianNeighbor::FOURD_LOW:
        Remove<CartesianNeighbor::FOURD_LOW>(n);
        break;
    case CartesianNeighbor::FOURD_UP:
        Remove<CartesianNeighbor::FOURD_UP>(n);
        break;
    }
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
template <CartesianNeighbor CNB>
SC MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n) const {
    if constexpr (IsSame<CartesianNeighbor::CENTER, CNB>()) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::CENTER)];
    } else if constexpr (IsSame<CartesianNeighbor::WEST, CNB>()) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::WEST)];
    } else if constexpr (IsSame<CartesianNeighbor::EAST, CNB>()) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::EAST)];
    } else if constexpr (IsSame<CartesianNeighbor::SOUTH, CNB>() && (Dim > 1)) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::SOUTH)];
    } else if constexpr (IsSame<CartesianNeighbor::NORTH, CNB>() && (Dim > 1)) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::NORTH)];
    } else if constexpr (IsSame<CartesianNeighbor::BOTTOM, CNB>() && (Dim > 2)) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::BOTTOM)];
    } else if constexpr (IsSame<CartesianNeighbor::TOP, CNB>() && (Dim > 2)) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::TOP)];
    } else if constexpr (IsSame<CartesianNeighbor::FOURD_LOW, CNB>() && (Dim > 3)) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::FOURD_LOW)];
    } else if constexpr (IsSame<CartesianNeighbor::FOURD_UP, CNB>() && (Dim > 3)) {
        return GetNeighbors()[n][static_cast<char>(CartesianNeighbor::FOURD_UP)];
    } else {
       ERROR <<"The neighbor indicator is out of bounds!" << ERROR_CLOSE;
    }
    return std::numeric_limits<SC>::signaling_NaN();
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
SC MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Get(std::size_t n, CartesianNeighbor cnb) const {
    switch (cnb) {
    case CartesianNeighbor::CENTER:
        return Get<CartesianNeighbor::CENTER>(n);
    case CartesianNeighbor::WEST:
        return Get<CartesianNeighbor::WEST>(n);
    case CartesianNeighbor::EAST:
        return Get<CartesianNeighbor::EAST>(n);
    case CartesianNeighbor::SOUTH:
        return Get<CartesianNeighbor::SOUTH>(n);
    case CartesianNeighbor::NORTH:
        return Get<CartesianNeighbor::NORTH>(n);
    case CartesianNeighbor::BOTTOM:
        return Get<CartesianNeighbor::BOTTOM>(n);
    case CartesianNeighbor::TOP:
        return Get<CartesianNeighbor::TOP>(n);
    case CartesianNeighbor::FOURD_LOW:
        return Get<CartesianNeighbor::FOURD_LOW>(n);
    case CartesianNeighbor::FOURD_UP:
        return Get<CartesianNeighbor::FOURD_UP>(n);
    }
    return std::numeric_limits<SC>::signaling_NaN();
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
dare::utils::Vector<N, typename MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::ScalarArray>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighbors() {
    return neighbors;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
const dare::utils::Vector<N, typename MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::ScalarArray>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighbors() const {
    return neighbors;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
dare::utils::Vector<N, char>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighborBitSet() {
    return neighbor_set;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
const dare::utils::Vector<N, char>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetNeighborBitSet() const {
    return neighbor_set;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Set(const dare::Data::CenterMatrixStencil<GridType, SC, N>& s) {
    for (std::size_t n{0}; n < N; n++) {
        Get<CartesianNeighbor::CENTER>(n) = s.GetValue(CartesianNeighbor::CENTER, n);
        Get<CartesianNeighbor::WEST>(n) = s.GetValue(CartesianNeighbor::WEST, n);
        Get<CartesianNeighbor::EAST>(n) = s.GetValue(CartesianNeighbor::EAST, n);
        if constexpr (Dim > 1) {
            Get<CartesianNeighbor::SOUTH>(n) = s.GetValue(CartesianNeighbor::SOUTH, n);
            Get<CartesianNeighbor::NORTH>(n) = s.GetValue(CartesianNeighbor::NORTH, n);
        }
        if constexpr (Dim > 2) {
            Get<CartesianNeighbor::BOTTOM>(n) = s.GetValue(CartesianNeighbor::BOTTOM, n);
            Get<CartesianNeighbor::TOP>(n) = s.GetValue(CartesianNeighbor::TOP, n);
        }
        this->GetRhs(n) = s.GetRHS()[n];
    }
    return *this;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Set(std::size_t n,
                                                       const dare::Data::CenterMatrixStencil<GridType, SC, 1>& s) {
    Get<CartesianNeighbor::CENTER>(n) = s.GetValue(CartesianNeighbor::CENTER, 0);
    Get<CartesianNeighbor::WEST>(n) = s.GetValue(CartesianNeighbor::WEST, 0);
    Get<CartesianNeighbor::EAST>(n) = s.GetValue(CartesianNeighbor::EAST, 0);
    if constexpr (Dim > 1) {
        Get<CartesianNeighbor::SOUTH>(n) = s.GetValue(CartesianNeighbor::SOUTH, 0);
        Get<CartesianNeighbor::NORTH>(n) = s.GetValue(CartesianNeighbor::NORTH, 0);
    }
    if constexpr (Dim > 2) {
        Get<CartesianNeighbor::BOTTOM>(n) = s.GetValue(CartesianNeighbor::BOTTOM, 0);
        Get<CartesianNeighbor::TOP>(n) = s.GetValue(CartesianNeighbor::TOP, 0);
    }
    this->GetRhs(n) = s.GetRHS()[0];
    return *this;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator=(
    const dare::Data::CenterMatrixStencil<GridType, SC, N>& s) {
    return this->Set(s);
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Add(
    const dare::Data::CenterMatrixStencil<GridType, SC, N>& s) {
    for (std::size_t n{0}; n < N; n++) {
        Get<CartesianNeighbor::CENTER>(n) += s.GetValue(CartesianNeighbor::CENTER, n);
        Get<CartesianNeighbor::WEST>(n) += s.GetValue(CartesianNeighbor::WEST, n);
        Get<CartesianNeighbor::EAST>(n) += s.GetValue(CartesianNeighbor::EAST, n);
        if constexpr (Dim > 1) {
            Get<CartesianNeighbor::SOUTH>(n) += s.GetValue(CartesianNeighbor::SOUTH, n);
            Get<CartesianNeighbor::NORTH>(n) += s.GetValue(CartesianNeighbor::NORTH, n);
        }
        if constexpr (Dim > 2) {
            Get<CartesianNeighbor::BOTTOM>(n) += s.GetValue(CartesianNeighbor::BOTTOM, n);
            Get<CartesianNeighbor::TOP>(n) += s.GetValue(CartesianNeighbor::TOP, n);
        }
        this->GetRhs(n) += s.GetRHS()[n];
    }
    return *this;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Add(std::size_t n,
                                                       const dare::Data::CenterMatrixStencil<GridType, SC, 1>& s) {
    Get<CartesianNeighbor::CENTER>(n) += s.GetValue(CartesianNeighbor::CENTER, 0);
    Get<CartesianNeighbor::WEST>(n) += s.GetValue(CartesianNeighbor::WEST, 0);
    Get<CartesianNeighbor::EAST>(n) += s.GetValue(CartesianNeighbor::EAST, 0);
    if constexpr (Dim > 1) {
        Get<CartesianNeighbor::SOUTH>(n) += s.GetValue(CartesianNeighbor::SOUTH, 0);
        Get<CartesianNeighbor::NORTH>(n) += s.GetValue(CartesianNeighbor::NORTH, 0);
    }
    if constexpr (Dim > 2) {
        Get<CartesianNeighbor::BOTTOM>(n) += s.GetValue(CartesianNeighbor::BOTTOM, 0);
        Get<CartesianNeighbor::TOP>(n) += s.GetValue(CartesianNeighbor::TOP, 0);
    }
    this->GetRhs(n) += s.GetRHS()[0];
    return *this;
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator+=(
    const dare::Data::CenterMatrixStencil<GridType, SC, N>& s) {
    return this->Add(s);
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Subtract(
    const dare::Data::CenterMatrixStencil<GridType, SC, N>& s) {
    return this->Add(-1. * s);
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::Subtract(
    std::size_t n,
    const dare::Data::CenterMatrixStencil<GridType, SC, 1>& s) {
    return this->Add(n, -1. * s);
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>&
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::operator-=(
    const dare::Data::CenterMatrixStencil<GridType, SC, N>& s) {
    return (*this) += (-1. * s);
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
typename MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GO
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetGlobalOrdinal() const {
    if constexpr (std::is_same_v<O, GlobalOrdinalType>)
        return this->GetNode();
    else
        return g_rep->MapLocalToGlobalInternal(this->GetNode());
}

template <std::size_t Dim, typename O, typename SC, std::size_t N>
typename MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::LO
MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>::GetLocalOrdinal() const {
    if constexpr (std::is_same_v<O, LocalOrdinalType>)
        return this->GetNode();
    else
        return g_rep->MapGlobalToLocalInternal(this->GetNode());
}

}  // end namespace dare::Matrix

template <std::size_t Dim, typename O, typename SC, std::size_t N>
void std::swap(dare::Matrix::MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>& m1,
               dare::Matrix::MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>& m2) {
    std::swap(static_cast<dare::Matrix::MatrixBlockBase<O, SC, N>&>(m1),
              static_cast<dare::Matrix::MatrixBlockBase<O, SC, N>&>(m2));
    std::swap(m1.g_rep, m2.g_rep);
    std::swap(m1.neighbors, m2.neighbors);
    std::swap(m1.neighbor_set, m2.neighbor_set);
    std::swap(m1.ind_internal, m2.ind_internal);
    std::swap(m1.ind_full, m2.ind_full);
}
