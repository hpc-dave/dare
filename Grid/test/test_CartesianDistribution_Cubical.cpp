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

#include <gtest/gtest.h>
#include "test_DistributionFunctions.h"
#include "Utilities/Vector.h"
#include "Grid/Cartesian.h"

TEST_P(ConsistencyTest, CubicalOneDim) {
    dare::mpi::ExecutionManager exman;
    if (!exman.AmIRoot()) {
        SUCCEED();
    }
    int num_proc = GetParam();
    const std::size_t Dim{1};
    using Grid = dare::Grid::Cartesian<Dim>;
    using VecLO = typename Grid::VecLO;
    using VecGO = typename Grid::VecGO;
    VecGO resolution_global(num_proc*5);
    std::vector<VecLO> vec_res_local;
    std::vector<VecGO> vec_offsets;
    dare::Grid::details::CartesianDistribution_Cubical(num_proc,
                                                      resolution_global,
                                                      &vec_res_local, &vec_offsets,
                                                      false);

    dare::Grid::test::details::TestSumCells(resolution_global, vec_res_local);
}

TEST_P(ConsistencyTest, CubicalTwoDim) {
    dare::mpi::ExecutionManager exman;
    if (!exman.AmIRoot()) {
        SUCCEED();
    }
    int num_proc = GetParam();
    const std::size_t Dim{2};
    using Grid = dare::Grid::Cartesian<Dim>;
    using VecLO = typename Grid::VecLO;
    using VecGO = typename Grid::VecGO;
    VecGO resolution_global(num_proc * 5, num_proc * 5);
    std::vector<VecLO> vec_res_local;
    std::vector<VecGO> vec_offsets;
    dare::Grid::details::CartesianDistribution_Cubical(num_proc,
                                                      resolution_global,
                                                      &vec_res_local, &vec_offsets,
                                                      false);

    dare::Grid::test::details::TestSumCells(resolution_global, vec_res_local);
}

TEST_P(ConsistencyTest, CubicalThreeDim) {
    dare::mpi::ExecutionManager exman;
    if (!exman.AmIRoot()) {
        SUCCEED();
    }
    int num_proc = GetParam();
    const std::size_t Dim{3};
    using Grid = dare::Grid::Cartesian<Dim>;
    using VecLO = typename Grid::VecLO;
    using VecGO = typename Grid::VecGO;
    VecGO resolution_global(num_proc * 5, num_proc * 5, num_proc * 5);
    std::vector<VecLO> vec_res_local;
    std::vector<VecGO> vec_offsets;
    dare::Grid::details::CartesianDistribution_Cubical(num_proc,
                                                      resolution_global,
                                                      &vec_res_local, &vec_offsets,
                                                      false);

    dare::Grid::test::details::TestSumCells(resolution_global, vec_res_local);
}

INSTANTIATE_TEST_SUITE_P(CartesianDistributionTest,
                         ConsistencyTest,
                         testing::Values(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 17,
                         20, 23, 28, 32, 50, 75, 100, 192, 1013));
