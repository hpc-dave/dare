/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder
 *
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

#include <gtest/gtest.h>

#include "../MatrixBlock.h"

namespace dare::Matrix::test {
class Grid {
public:
    using GlobalOrdinalType = int64_t;
    using LocalOrdinalType = int32_t;
    using ScalarType = double;
    class TestRepresentation {
        public:
        using GlobalOrdinalType = Grid::GlobalOrdinalType;
        using LocalOrdinalType = Grid::LocalOrdinalType;
        const GlobalOrdinalType offset{10};
        GlobalOrdinalType MapLocalToGlobalInternal(LocalOrdinalType node) {
            return node + offset;
        }

        LocalOrdinalType MapGlobalToLocalInternal(GlobalOrdinalType node) {
            return node - offset;
        }
    };
    using Representation = TestRepresentation;

    Representation GetRepresentation() { return Representation(); }
};

}  // end namespace dare::Matrix::test

class MatrixBlockTest : public testing::Test {
public:
    using GridType = dare::Matrix::test::Grid;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = double;
    using GridRepresentation = typename GridType::Representation;
    const std::size_t num_ghost = 2;
    static const std::size_t N = 3;
    GridType grid;
};

TEST_F(MatrixBlockTest, Initialization) {
    GridRepresentation g_rep = grid.GetRepresentation();
    LO node = 11;
    dare::utils::Vector<N, std::size_t> size_hint;
    for (auto& e : size_hint)
        e = 3;
    dare::Matrix::MatrixBlock<GridType, LO, SC, N> mblock(&g_rep, node, size_hint);
    dare::Matrix::MatrixBlock<GridType, LO, SC, N> mblock_copy_construct(mblock);
    dare::Matrix::MatrixBlock<GridType, LO, SC, N> mblock_copy_assign;
    mblock_copy_assign = mblock;

    EXPECT_FALSE(mblock.IsGlobal());
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mblock.GetNumEntries(n), mblock_copy_construct.GetNumEntries(n));
        EXPECT_EQ(mblock.GetNumEntries(n), mblock_copy_assign.GetNumEntries(n));
    }
}

TEST_F(MatrixBlockTest, Convert) {
    GridRepresentation g_rep = grid.GetRepresentation();
    LO node_local = 11;
    GO node_global = node_local + g_rep.offset;
    dare::utils::Vector<N, std::size_t> size_hint;
    for (auto& e : size_hint)
        e = 3;
    dare::Matrix::MatrixBlock<GridType, LO, SC, N> mblock_local(&g_rep, node_local, size_hint);
    dare::Matrix::MatrixBlock<GridType, GO, SC, N> mblock_global(&g_rep, node_global, size_hint);

    for (std::size_t n{0}; n < N; n++) {
        for (std::size_t i{0}; i < size_hint[n]; i++) {
            mblock_local.SetCoefficient(n, node_local * N + n + i, 1.);
            mblock_global.SetCoefficient(n, node_global * N + n + i, 1.);
        }
    }

    dare::Matrix::MatrixBlock<GridType, LO, SC, N> mblock_loc_conv = dare::Matrix::Convert<LO>(mblock_global);
    dare::Matrix::MatrixBlock<GridType, GO, SC, N> mblock_glob_conv = dare::Matrix::Convert<GO>(mblock_local);
    dare::Matrix::MatrixBlock<GridType, LO, SC, N> mblock_no_conv = dare::Matrix::Convert<LO>(mblock_local);
    EXPECT_EQ(mblock_loc_conv.GetNode(), node_local);
    EXPECT_EQ(mblock_glob_conv.GetNode(), node_global);
    EXPECT_EQ(mblock_no_conv.GetNode(), node_local);

    for (std::size_t n{0}; n < N; n++) {
        for (std::size_t i{0}; i < size_hint[n]; i++) {
            EXPECT_EQ(mblock_loc_conv.GetOrdinalByPosition(n, i), node_local * N + n + i);
            EXPECT_EQ(mblock_glob_conv.GetOrdinalByPosition(n, i), node_global * N + n + i);
            EXPECT_EQ(mblock_no_conv.GetOrdinalByPosition(n, i), node_local * N + n + i);
        }
    }
}
