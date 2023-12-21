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

#include <gtest/gtest.h>
#include "../Grid/Cartesian.h"
#include "../MatrixSystem/MatrixBlock.h"

template <std::size_t Dim, typename GO>
dare::utils::Vector<Dim, GO> GetResolutionTestCartesianGrid() {
    dare::utils::Vector<Dim, GO> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 100 + n;
    return res;
}
template <std::size_t Dim, typename SC>
dare::utils::Vector<Dim, SC> GetSizeTestCartesianGrid() {
    dare::utils::Vector<Dim, SC> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

TEST(IntegrationCartesianMatrixBlock, Initialization) {
    static const std::size_t Dim{3};
    static const std::size_t N{3};
    const std::size_t num_ghost{2};
    using GridType = dare::Grid::Cartesian<Dim>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using SC = double;
    dare::mpi::ExecutionManager exec_man;
    GridType grid(&exec_man, GetResolutionTestCartesianGrid<Dim, GO>(), GetSizeTestCartesianGrid<Dim, SC>(), num_ghost);

    dare::utils::Vector<N, std::size_t> size_hint;
    for (auto& e : size_hint) {
        e = 5;
    }
    GridType::Options opt(0);  // not staggered
    GridType::Representation g_rep = grid.GetRepresentation(opt);
    IndexGlobal ind_glob = grid.GetOffsetCells();
    for (auto& e : ind_glob)
        e += 2;
    Index ind_loc = g_rep.MapGlobalToLocal(ind_glob);
    dare::Matrix::MatrixBlock<GridType, LO, SC, N> mblock_local(&g_rep,
                g_rep.MapIndexToOrdinalLocalInternal(g_rep.MapLocalToInternal(ind_loc)), size_hint);
    dare::Matrix::MatrixBlock<GridType, GO, SC, N> mblock_global(&g_rep,
                g_rep.MapIndexToOrdinalGlobalInternal(g_rep.MapGlobalToInternal(ind_loc)), size_hint);
}
