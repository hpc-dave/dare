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
MatrixBlock<Grid, O, SC, N>::MatrixBlock(GridRepresentation* _g_rep,
                                         Ordinal Node,
                                         const dare::utils::Vector<N, std::size_t>& size_hint)
    : MatrixBlockBase<O, SC, N>(node, size_hint), g_rep(_g_rep) {
    static_assert(std::is_same_v<O, tpyename GridRepresentation::LocalOrdinalType>
               || std::is_same_v<O, tpyename GridRepresentation::GlobalOrdinalType>,
                  "The ordinal type needs to be either a local or global ordinal type!");
}

template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>::MatrixBlock(const SelfType& other)
    : MatrixBlockBase<O, SC, N>(other), g_rep(other.g_rep) {
    static_assert(std::is_same_v<O, tpyename GridRepresentation::LocalOrdinalType>
               || std::is_same_v<O, tpyename GridRepresentation::GlobalOrdinalType>,
                  "The ordinal type needs to be either a local or global ordinal type!");
}

template <typename Grid, typename O, typename SC, std::size_t N>
MatrixBlock<Grid, O, SC, N>::~MatrixBlock() {}

template <typename Grid, typename O, typename SC, std::size_t N>
SelfType& MatrixBlock<Grid, O, SC, N>::operator=(const SelfType& other) {
    MatrixBase<O, SC, N>::operator=(other);
    g_rep = other.g_rep;
}

template <typename Grid, typename O, typename SC, std::size_t N>
typename MatrixBlock<Grid, O, SC, N>::GridRepresentation*
MatrixBlock<Grid, O, SC, N>::GetRepresentation() const {
    return g_rep;
}

template <typename Grid, typename O, typename SC, std::size_t N>
bool MatrixBlock<Grid, O, SC, N>::IsGlobal() const {
    return std::is_same_v<O, GlobalOrdinalType>;
}

}  // end namespace dare::Matrix
