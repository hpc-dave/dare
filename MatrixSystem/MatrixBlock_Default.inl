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
template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>::MatrixBlock()
    : MatrixBlock<Grid, O, SC, N>(nullptr, 0, dare::utils::Vector<N, std::size_t>()) {
}

template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>::MatrixBlock(const GridRepresentation* _g_rep,
                                         O _node)
    : MatrixBlock<Grid, O, SC, N>(_g_rep, _node, dare::utils::Vector<N, std::size_t>()) {
}

template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>::MatrixBlock(const GridRepresentation* _g_rep,
                                         O _node,
                                         const dare::utils::Vector<N, std::size_t>& size_hint)
    : MatrixBlockBase<O, SC, N>(_node, size_hint), g_rep(_g_rep) {
    static_assert(std::is_same_v<O, typename Grid::LocalOrdinalType>
               || std::is_same_v<O, typename Grid::GlobalOrdinalType>,
                  "The ordinal type needs to be either a local or global ordinal type!");
}

template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>::MatrixBlock(const SelfType& other)
    : MatrixBlockBase<O, SC, N>(other), g_rep(other.g_rep) {
    static_assert(std::is_same_v<O, typename Grid::LocalOrdinalType>
               || std::is_same_v<O, typename Grid::GlobalOrdinalType>,
                  "The ordinal type needs to be either a local or global ordinal type!");
}

template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>::~MatrixBlock() {}

template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>& MatrixBlock<Grid, O, SC, N>::operator=(const MatrixBlock<Grid, O, SC, N>& other) {
    if (&other == this)
        return *this;
    MatrixBlockBase<O, SC, N>::operator=(other);
    g_rep = other.g_rep;
    return *this;
}

template <typename Grid, typename O, typename SC, std::size_t N>
void MatrixBlock<Grid, O, SC, N>::Initialize(const GridRepresentation* _g_rep, O _node) {
    Initialize(_g_rep, _node, dare::utils::Vector<N, std::size_t>());
}

template <typename Grid, typename O, typename SC, std::size_t N>
void MatrixBlock<Grid, O, SC, N>::Initialize(const GridRepresentation* _g_rep,
                                             O _node,
                                             const dare::utils::Vector<N, std::size_t>& size_hint) {
    MatrixBlockBase<O, SC, N>::Initialize(_node, size_hint);
    g_rep = _g_rep;
}

template <typename Grid, typename O, typename SC, std::size_t N>
const typename MatrixBlock<Grid, O, SC, N>::GridRepresentation*
MatrixBlock<Grid, O, SC, N>::GetRepresentation() const {
    return g_rep;
}

template <typename Grid, typename O, typename SC, std::size_t N>
bool MatrixBlock<Grid, O, SC, N>::IsGlobal() const {
    return std::is_same_v<O, GlobalOrdinalType>;
}

template <typename Grid, typename O, typename SC, std::size_t N>
bool MatrixBlock<Grid, O, SC, N>::IsStencilLocal() const {
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

}  // end namespace dare::Matrix
