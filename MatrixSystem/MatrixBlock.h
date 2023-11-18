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

#ifndef MATRIXSYSTEM_MATRIXBLOCK_H_
#define MATRIXSYSTEM_MATRIXBLOCK_H_

#include <type_traits>

#include "MatrixBlockBase.h"
#include "../Utilities/Vector.h"

namespace dare::Matrix {

template<typename Grid, typename O, typename SC, std::size_t N>
class MatrixBlock : public MatrixBlockBase<O, SC, N> {
public:
    using GridRepresentation = typename Grid::Representation;
    using LocalOrdinalType = typename Grid::LocalOrdinalType;
    using GlobalOrdinalType = typename Grid::GlobalOrdinalType;
    using SelfType = MatrixBlock<Grid O, SC, N>;
    MatrixBlock();
    MatrixBlock(GridRepresentation* g_rep, Ordinal Node, const dare::utils::Vector<N, std::size_t>& size_hint);
    explicit MatrixBlock(const SelfType& other);
    virtual ~MatrixBlock();

    SelfType& operator=(const SelfType& other);

    bool IsGlobal() const;

    GridRepresentation* GetRepresentation() const;

private:
    GridRepresentation* g_rep;
};
}  // end namespace dare::Matrix

#include "MatrixBlock.inl"

#endif  // MATRIXSYSTEM_MATRIXBLOCK_H_
