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

#include "../MatrixSystem/Operators_Cartesian.h"

namespace dare::test {

template <std::size_t Dim, typename GO>
dare::utils::Vector<Dim, GO> GetResolutionTestCartesianOperators() {
    dare::utils::Vector<Dim, GO> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 10 + n;
    return res;
}
template <std::size_t Dim, typename SC>
dare::utils::Vector<Dim, SC> GetSizeTestCartesianOperators() {
    dare::utils::Vector<Dim, SC> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

}  // namespace dare::test

template <std::size_t Dim>
class IntegrationTestCartesianOperators : public testing::Test {
public:
    static const std::size_t N{3};
    using GridType = dare::Grid::Cartesian<Dim>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using Index = typename GridType::Index;
    using CenterMatrixStencil = dare::Data::CenterMatrixStencil<GridType, N>;
    using CenterValueStencil = dare::Data::CenterValueStencil<GridType, N>;
    using FaceMatrixStencil = dare::Data::FaceMatrixStencil<GridType, N>;
    using FaceValueStencil = dare::Data::FaceValueStencil<GridType, N>;
    using Gradient = dare::Matrix::Gradient<GridType>;
    using Divergence = dare::Matrix::Divergence<GridType>;
    using MatrixBlock = dare::Matrix::MatrixBlock<GridType, LO, SC, N>;
    using VecSC = dare::utils::Vector<Dim, SC>;
    using Positions = typename Gradient::Positions;
    using Field = dare::Data::GridVector<GridType, SC, N>;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(&exec_man,
                                          dare::test::GetResolutionTestCartesianOperators<Dim, GO>(),
                                          dare::test::GetSizeTestCartesianOperators<Dim, SC>(),
                                          num_ghost);
    }

    std::unique_ptr<GridType> grid;         //!< the grid
    dare::mpi::ExecutionManager exec_man;   //!< the execution manager
};

using IntegrationTestCartesianOperators1D = IntegrationTestCartesianOperators<1>;
using IntegrationTestCartesianOperators2D = IntegrationTestCartesianOperators<2>;
using IntegrationTestCartesianOperators3D = IntegrationTestCartesianOperators<3>;

TEST_F(IntegrationTestCartesianOperators1D, GradientMatrixValues) {
    GridType::Options opt{0};   // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    MatrixBlock mb(&grid_rep, ordinal_internal);
    Gradient grad(grid_rep, ordinal_internal);

    auto s = grad(mb);
    static_assert(std::is_same_v<decltype(s), FaceMatrixStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1.) / grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s.GetValueCenter(Positions::WEST, n), dn_r[0]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::WEST, n), -dn_r[0]);

        EXPECT_EQ(s.GetValueCenter(Positions::EAST, n), -dn_r[0]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::EAST, n), dn_r[0]);
    }
}

TEST_F(IntegrationTestCartesianOperators2D, GradientMatrixValues) {
    GridType::Options opt{0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    MatrixBlock mb(&grid_rep, ordinal_internal);
    Gradient grad(grid_rep, ordinal_internal);

    auto s = grad(mb);
    static_assert(std::is_same_v<decltype(s), FaceMatrixStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1., 1.) / grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s.GetValueCenter(Positions::WEST, n), dn_r[0]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::WEST, n), -dn_r[0]);

        EXPECT_EQ(s.GetValueCenter(Positions::EAST, n), -dn_r[0]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::EAST, n), dn_r[0]);

        EXPECT_EQ(s.GetValueCenter(Positions::SOUTH, n), dn_r[1]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::SOUTH, n), -dn_r[1]);

        EXPECT_EQ(s.GetValueCenter(Positions::NORTH, n), -dn_r[1]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::NORTH, n), dn_r[1]);
    }
}

TEST_F(IntegrationTestCartesianOperators3D, GradientMatrixValues) {
    GridType::Options opt{0, 0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    MatrixBlock mb(&grid_rep, ordinal_internal);
    Gradient grad(grid_rep, ordinal_internal);

    auto s = grad(mb);
    static_assert(std::is_same_v<decltype(s), FaceMatrixStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1., 1., 1.) / grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s.GetValueCenter(Positions::WEST, n), dn_r[0]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::WEST, n), -dn_r[0]);

        EXPECT_EQ(s.GetValueCenter(Positions::EAST, n), -dn_r[0]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::EAST, n), dn_r[0]);

        EXPECT_EQ(s.GetValueCenter(Positions::SOUTH, n), dn_r[1]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::SOUTH, n), -dn_r[1]);

        EXPECT_EQ(s.GetValueCenter(Positions::NORTH, n), -dn_r[1]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::NORTH, n), dn_r[1]);

        EXPECT_EQ(s.GetValueCenter(Positions::BOTTOM, n), dn_r[2]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::BOTTOM, n), -dn_r[2]);

        EXPECT_EQ(s.GetValueCenter(Positions::TOP, n), -dn_r[2]);
        EXPECT_EQ(s.GetValueNeighbor(Positions::TOP, n), dn_r[2]);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, GradientFaceValuesFromField) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Gradient grad(grid_rep, ordinal_internal);

    Field field("centers", grid_rep);

    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * n;
        }
    }

    auto s = grad(field);
    static_assert(std::is_same_v<decltype(s), FaceValueStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1.) / grid_rep.GetDistances();

    LO ordinal_local = grid_rep.MapInternalToLocal(ordinal_internal);
    Index ind = grid_rep.MapOrdinalToIndexLocal(ordinal_local);
    Index ind_w{ind}, ind_e{ind};
    ind_w.i()--;
    ind_e.i()++;
    for (std::size_t n{0}; n < N; n++) {
        double v_west = field.At(ind_w, n);
        double v_center = field.At(ind, n);
        double v_east = field.At(ind_e, n);
        EXPECT_EQ(s.GetValue(Positions::WEST, n), (v_center - v_west)   * dn_r[0]);
        EXPECT_EQ(s.GetValue(Positions::EAST, n), (v_east   - v_center) * dn_r[0]);
    }
}

TEST_F(IntegrationTestCartesianOperators2D, GradientFaceValuesFromField) {
    GridType::Options opt{0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Gradient grad(grid_rep, ordinal_internal);

    Field field("centers", grid_rep);

    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        for (LO j{0}; j < grid_rep.GetLocalResolution().j(); j++) {
            Index ind(i, j);
            for (std::size_t n{0}; n < N; n++) {
                field.At(ind, n) = j * n;
            }
        }
    }

    auto s = grad(field);
    static_assert(std::is_same_v<decltype(s), FaceValueStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1., 1.) / grid_rep.GetDistances();

    LO ordinal_local = grid_rep.MapInternalToLocal(ordinal_internal);
    Index ind = grid_rep.MapOrdinalToIndexLocal(ordinal_local);
    Index ind_s{ind}, ind_n{ind};
    ind_s.j()--;
    ind_n.j()++;
    for (std::size_t n{0}; n < N; n++) {
        double v_south = field.At(ind_s, n);
        double v_center = field.At(ind, n);
        double v_north = field.At(ind_n, n);
        EXPECT_EQ(s.GetValue(Positions::WEST, n), 0.);
        EXPECT_EQ(s.GetValue(Positions::EAST, n), 0.);
        EXPECT_EQ(s.GetValue(Positions::SOUTH, n), (v_center - v_south) * dn_r[1]);
        EXPECT_EQ(s.GetValue(Positions::NORTH, n), (v_north - v_center) * dn_r[1]);
    }
}

TEST_F(IntegrationTestCartesianOperators3D, GradientFaceValuesFromField) {
    GridType::Options opt{0, 0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Gradient grad(grid_rep, ordinal_internal);

    Field field("centers", grid_rep);

    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        for (LO j{0}; j < grid_rep.GetLocalResolution().j(); j++) {
            for (LO k{0}; k < grid_rep.GetLocalResolution().k(); k++) {
                Index ind(i, j, k);
                for (std::size_t n{0}; n < N; n++) {
                    field.At(ind, n) = k * n;
                }
            }
        }
    }

    auto s = grad(field);
    static_assert(std::is_same_v<decltype(s), FaceValueStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1., 1., 1.) / grid_rep.GetDistances();

    LO ordinal_local = grid_rep.MapInternalToLocal(ordinal_internal);
    Index ind = grid_rep.MapOrdinalToIndexLocal(ordinal_local);
    Index ind_bot{ind}, ind_top{ind};
    ind_bot.k()--;
    ind_top.k()++;
    for (std::size_t n{0}; n < N; n++) {
        double v_bottom = field.At(ind_bot, n);
        double v_center = field.At(ind, n);
        double v_top = field.At(ind_top, n);
        EXPECT_EQ(s.GetValue(Positions::WEST, n), 0.);
        EXPECT_EQ(s.GetValue(Positions::EAST, n), 0.);
        EXPECT_EQ(s.GetValue(Positions::SOUTH, n), 0.);
        EXPECT_EQ(s.GetValue(Positions::NORTH, n), 0.);
        EXPECT_EQ(s.GetValue(Positions::BOTTOM, n), (v_center - v_bottom) * dn_r[2]);
        EXPECT_EQ(s.GetValue(Positions::TOP, n), (v_top - v_center) * dn_r[2]);
    }
}
