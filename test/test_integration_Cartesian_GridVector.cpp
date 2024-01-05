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

#include <memory>

#include "../Grid/Cartesian.h"
#include "../Data/GridVector.h"

namespace dare::test {

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
}  // namespace dare::test

template <std::size_t Dim>
class IntegrationTestCartesianGridVector : public testing::Test {
public:
    using GridType = dare::Grid::Cartesian<Dim>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using SC = double;
    static const std::size_t N{3};
    using GridVector = dare::Data::GridVector<GridType, SC, N>;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(&exec_man,
                                          dare::test::GetResolutionTestCartesianGrid<Dim, GO>(),
                                          dare::test::GetSizeTestCartesianGrid<Dim, SC>(),
                                          num_ghost);
    }

    std::unique_ptr<GridType> grid;
    dare::mpi::ExecutionManager exec_man;
};

using IntegrationTestCartesianGridVector1D = IntegrationTestCartesianGridVector<1>;
using IntegrationTestCartesianGridVector2D = IntegrationTestCartesianGridVector<2>;
using IntegrationTestCartesianGridVector3D = IntegrationTestCartesianGridVector<3>;

TEST_F(IntegrationTestCartesianGridVector1D, Initialize) {
    GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector2D, Initialize) {
    GridType::Options opt(0, 0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector3D, Initialize) {
    GridType::Options opt(0, 0, 0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector1D, InitializeStaggered) {
    GridType::Options opt(1);  // staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector2D, InitializeStaggeredX) {
    GridType::Options opt(1, 0);  // staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector2D, InitializeStaggeredY) {
    GridType::Options opt(0, 1);  // staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector3D, InitializeStaggeredX) {
    GridType::Options opt(1, 0, 0);  // staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector3D, InitializeStaggeredY) {
    GridType::Options opt(0, 1, 0);  // staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector3D, InitializeStaggeredZ) {
    GridType::Options opt(0, 0, 1);  // staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    LO num_data_entries = data.GetSize();
    LO num_expected_entries = g_rep.GetNumberLocalCells() * N;
    EXPECT_EQ(num_data_entries, num_expected_entries);
}

TEST_F(IntegrationTestCartesianGridVector3D, UseIndex) {
    GridType::Options opt(0, 0, 0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    for (std::size_t n{0}; n < data.GetSize(); n++) {
        data[n] = n;
    }
    GridType::Index res = g_rep.GetLocalResolution();
    for (LO i{0}; i < res.i(); i++) {
        for (LO j{0}; j < res.j(); j++) {
            for (LO k{0}; k < res.k(); k++) {
                GridType::Index ind(i, j, k);
                LO id = g_rep.MapIndexToOrdinalLocal(ind);
                for (std::size_t n{0}; n < N; n++) {
                    SC value_expect = id * N + n;
                    SC value_found = data.At(ind, n);
                    EXPECT_EQ(value_expect, value_found);
                }
            }
        }
    }
}
