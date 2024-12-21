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

#include "Grid/Cartesian/Operators_Cartesian.h"
#include "Equations/TimeDiscretizationSchemes.h"

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
    using TimeDisc = dare::Matrix::EULER_BACKWARD;
    using Index = typename GridType::Index;
    using CenterMatrixStencil = dare::Data::CenterMatrixStencil<GridType, SC, N>;
    using CenterMatrixStencil1C = dare::Data::CenterMatrixStencil<GridType, SC, 1>;
    using CenterValueStencil = dare::Data::CenterValueStencil<GridType, SC, N>;
    using FaceMatrixStencil = dare::Data::FaceMatrixStencil<GridType, SC, N>;
    using FaceValueStencil = dare::Data::FaceValueStencil<GridType, SC, N>;
    using Gradient = dare::Matrix::Gradient<GridType>;
    using Divergence = dare::Matrix::Divergence<GridType, TimeDisc>;
    template <typename FluxLimiter>
    using TVD = dare::Matrix::TVD<GridType, SC, FluxLimiter>;
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

TEST_F(IntegrationTestCartesianOperators1D, TypeTraitTest) {
    static_assert(dare::is_center_matrix_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_center_value_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_face_value_stencil_v<CenterMatrixStencil>);

    static_assert(dare::is_face_matrix_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_face_value_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_center_value_stencil_v<FaceMatrixStencil>);

    static_assert(dare::is_center_value_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_face_value_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<CenterValueStencil>);

    static_assert(dare::is_face_value_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_center_value_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<FaceValueStencil>);
}

TEST_F(IntegrationTestCartesianOperators2D, TypeTraitTest) {
    static_assert(dare::is_center_matrix_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_center_value_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_face_value_stencil_v<CenterMatrixStencil>);

    static_assert(dare::is_face_matrix_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_face_value_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_center_value_stencil_v<FaceMatrixStencil>);

    static_assert(dare::is_center_value_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_face_value_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<CenterValueStencil>);

    static_assert(dare::is_face_value_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_center_value_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<FaceValueStencil>);
}

TEST_F(IntegrationTestCartesianOperators3D, TypeTraitTest) {
    static_assert(dare::is_center_matrix_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_center_value_stencil_v<CenterMatrixStencil>);
    static_assert(!dare::is_face_value_stencil_v<CenterMatrixStencil>);

    static_assert(dare::is_face_matrix_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_face_value_stencil_v<FaceMatrixStencil>);
    static_assert(!dare::is_center_value_stencil_v<FaceMatrixStencil>);

    static_assert(dare::is_center_value_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_face_value_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<CenterValueStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<CenterValueStencil>);

    static_assert(dare::is_face_value_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_face_matrix_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_center_value_stencil_v<FaceValueStencil>);
    static_assert(!dare::is_center_matrix_stencil_v<FaceValueStencil>);
}

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

    auto s1 = grad(field, 1);
    static_assert(std::is_same_v<decltype(s1), dare::Data::FaceValueStencil<GridType, SC, 1>>, "Type is wrong!");
    double v_west = field.At(ind_w, 1);
    double v_center = field.At(ind, 1);
    double v_east = field.At(ind_e, 1);
    EXPECT_EQ(s1.GetValue(Positions::WEST, 0), (v_center - v_west) * dn_r[0]);
    EXPECT_EQ(s1.GetValue(Positions::EAST, 0), (v_east - v_center) * dn_r[0]);
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

    auto s1 = grad(field, 1);
    static_assert(std::is_same_v<decltype(s1), dare::Data::FaceValueStencil<GridType, SC, 1>>, "Type is wrong!");
    std::size_t n = 1;
    double v_south = field.At(ind_s, n);
    double v_center = field.At(ind, n);
    double v_north = field.At(ind_n, n);
    EXPECT_EQ(s.GetValue(Positions::WEST, n), 0.);
    EXPECT_EQ(s.GetValue(Positions::EAST, n), 0.);
    EXPECT_EQ(s.GetValue(Positions::SOUTH, n), (v_center - v_south) * dn_r[1]);
    EXPECT_EQ(s.GetValue(Positions::NORTH, n), (v_north - v_center) * dn_r[1]);
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

    auto s1 = grad(field, 1);
    static_assert(std::is_same_v<decltype(s1), dare::Data::FaceValueStencil<GridType, SC, 1>>, "Type is wrong!");
    std::size_t n = 1;

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

TEST_F(IntegrationTestCartesianOperators1D, GradientFaceValuesFromStencil) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Gradient grad(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_center = 1.;
    double v_east = 2.;
    CenterValueStencil s_c;
    for (std::size_t n{0}; n < N; n++) {
        s_c.SetValue(Positions::WEST, n, v_west * n);
        s_c.SetValue(Positions::CENTER, n, v_center * n);
        s_c.SetValue(Positions::EAST, n, v_east * n);
    }

    auto s = grad(s_c);
    static_assert(std::is_same_v<decltype(s), FaceValueStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1.) / grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s.GetValue(Positions::WEST, n), (v_center - v_west) * n * dn_r[0]);
        EXPECT_EQ(s.GetValue(Positions::EAST, n), (v_east - v_center) * n * dn_r[0]);
    }

    auto s1 = grad(s_c, 1);
    static_assert(std::is_same_v<decltype(s1), dare::Data::FaceValueStencil<GridType, SC, 1>>, "Type is wrong!");
    EXPECT_EQ(s1.GetValue(Positions::WEST, 0), (v_center - v_west) * dn_r[0]);
    EXPECT_EQ(s1.GetValue(Positions::EAST, 0), (v_east - v_center) * dn_r[0]);
}

TEST_F(IntegrationTestCartesianOperators2D, GradientFaceValuesFromStencil) {
    GridType::Options opt{0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Gradient grad(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_center = 1.;
    double v_east = 2.;
    double v_south = -1.;
    double v_north = -2.;
    CenterValueStencil s_c;
    for (std::size_t n{0}; n < N; n++) {
        s_c.SetValue(Positions::WEST, n, v_west * n);
        s_c.SetValue(Positions::SOUTH, n, v_south * n);
        s_c.SetValue(Positions::CENTER, n, v_center * n);
        s_c.SetValue(Positions::NORTH, n, v_north * n);
        s_c.SetValue(Positions::EAST, n, v_east * n);
    }

    auto s = grad(s_c);
    static_assert(std::is_same_v<decltype(s), FaceValueStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1., 1.) / grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s.GetValue(Positions::WEST, n), (v_center - v_west) * n * dn_r[0]);
        EXPECT_EQ(s.GetValue(Positions::EAST, n), (v_east - v_center) * n * dn_r[0]);
        EXPECT_EQ(s.GetValue(Positions::SOUTH, n), (v_center - v_south) * n * dn_r[1]);
        EXPECT_EQ(s.GetValue(Positions::NORTH, n), (v_north - v_center) * n * dn_r[1]);
    }

    auto s1 = grad(s_c, 1);
    static_assert(std::is_same_v<decltype(s1), dare::Data::FaceValueStencil<GridType, SC, 1>>, "Type is wrong!");
    EXPECT_EQ(s1.GetValue(Positions::WEST, 0), (v_center - v_west) * dn_r[0]);
    EXPECT_EQ(s1.GetValue(Positions::EAST, 0), (v_east - v_center) * dn_r[0]);
    EXPECT_EQ(s1.GetValue(Positions::SOUTH, 0), (v_center - v_south) * dn_r[1]);
    EXPECT_EQ(s1.GetValue(Positions::NORTH, 0), (v_north - v_center) * dn_r[1]);
}

TEST_F(IntegrationTestCartesianOperators3D, GradientFaceValuesFromStencil) {
    GridType::Options opt{0, 0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Gradient grad(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_center = 1.;
    double v_east = 2.;
    double v_south = -1.;
    double v_north = -2.;
    double v_bottom = 3.;
    double v_top = 4.;
    CenterValueStencil s_c;
    for (std::size_t n{0}; n < N; n++) {
        s_c.SetValue(Positions::WEST, n, v_west * n);
        s_c.SetValue(Positions::SOUTH, n, v_south * n);
        s_c.SetValue(Positions::CENTER, n, v_center * n);
        s_c.SetValue(Positions::NORTH, n, v_north * n);
        s_c.SetValue(Positions::EAST, n, v_east * n);
        s_c.SetValue(Positions::BOTTOM, n, v_bottom * n);
        s_c.SetValue(Positions::TOP, n, v_top * n);
    }

    auto s = grad(s_c);
    static_assert(std::is_same_v<decltype(s), FaceValueStencil>, "Type is wrong!");

    VecSC dn_r = VecSC(1., 1., 1.) / grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s.GetValue(Positions::WEST, n), (v_center - v_west) * n * dn_r[0]);
        EXPECT_EQ(s.GetValue(Positions::EAST, n), (v_east - v_center) * n * dn_r[0]);
        EXPECT_EQ(s.GetValue(Positions::SOUTH, n), (v_center - v_south) * n * dn_r[1]);
        EXPECT_EQ(s.GetValue(Positions::NORTH, n), (v_north - v_center) * n * dn_r[1]);
        EXPECT_EQ(s.GetValue(Positions::BOTTOM, n), (v_center - v_bottom) * n * dn_r[2]);
        EXPECT_EQ(s.GetValue(Positions::TOP, n), (v_top - v_center) * n * dn_r[2]);
    }

    auto s1 = grad(s_c, 1);
    static_assert(std::is_same_v<decltype(s1), dare::Data::FaceValueStencil<GridType, SC, 1>>, "Type is wrong!");
    EXPECT_EQ(s1.GetValue(Positions::WEST, 0), (v_center - v_west) * dn_r[0]);
    EXPECT_EQ(s1.GetValue(Positions::EAST, 0), (v_east - v_center) * dn_r[0]);
    EXPECT_EQ(s1.GetValue(Positions::SOUTH, 0), (v_center - v_south) * dn_r[1]);
    EXPECT_EQ(s1.GetValue(Positions::NORTH, 0), (v_north - v_center) * dn_r[1]);
    EXPECT_EQ(s1.GetValue(Positions::BOTTOM, 0), (v_center - v_bottom) * dn_r[2]);
    EXPECT_EQ(s1.GetValue(Positions::TOP, 0), (v_top - v_center) * dn_r[2]);
}

TEST_F(IntegrationTestCartesianOperators1D, DivergenceFaceValueStencil) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Divergence div(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_east = 2.;
    FaceValueStencil s_f;
    for (std::size_t n{0}; n < N; n++) {
        s_f.SetValue(Positions::WEST, n, v_west * n);
        s_f.SetValue(Positions::EAST, n, v_east * n);
    }

    auto v_res = div(s_f);
    static_assert(std::is_same_v<decltype(v_res), dare::utils::Vector<N, SC>>, "Type is wrong!");

    VecSC A = grid_rep.GetFaceArea();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v_res[n], (v_east - v_west) * n* A[0]);
    }
}

TEST_F(IntegrationTestCartesianOperators2D, DivergenceFaceValueStencil) {
    GridType::Options opt{0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Divergence div(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_east = 2.;
    double v_south = 3.;
    double v_north = 4.;
    FaceValueStencil s_f;
    for (std::size_t n{0}; n < N; n++) {
        s_f.SetValue(Positions::WEST, n, v_west * n);
        s_f.SetValue(Positions::EAST, n, v_east * n);
        s_f.SetValue(Positions::SOUTH, n, v_south * n);
        s_f.SetValue(Positions::NORTH, n, v_north * n);
    }

    auto v_res = div(s_f);
    static_assert(std::is_same_v<decltype(v_res), dare::utils::Vector<N, SC>>, "Type is wrong!");

    VecSC A = grid_rep.GetFaceArea();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v_res[n], (v_east - v_west) * n * A[0]
                          + (v_north - v_south) * n * A[1]);
    }
}

TEST_F(IntegrationTestCartesianOperators3D, DivergenceFaceValueStencil) {
    GridType::Options opt{0, 0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Divergence div(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_east = 2.;
    double v_south = 3.;
    double v_north = 4.;
    double v_bottom = 5.;
    double v_top = 6.;
    FaceValueStencil s_f;
    for (std::size_t n{0}; n < N; n++) {
        s_f.SetValue(Positions::WEST, n, v_west * n);
        s_f.SetValue(Positions::EAST, n, v_east * n);
        s_f.SetValue(Positions::SOUTH, n, v_south * n);
        s_f.SetValue(Positions::NORTH, n, v_north * n);
        s_f.SetValue(Positions::BOTTOM, n, v_bottom * n);
        s_f.SetValue(Positions::TOP, n, v_top * n);
    }

    auto v_res = div(s_f);
    static_assert(std::is_same_v<decltype(v_res), dare::utils::Vector<N, SC>>, "Type is wrong!");

    VecSC A = grid_rep.GetFaceArea();
    VecSC dn_r;
    for (auto& e : dn_r)
        e = 1.;
    dn_r /= grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v_res[n], (v_east - v_west) * n * A[0]
                          + (v_north - v_south) * n * A[1]
                          + (v_top - v_bottom) * n * A[2]);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, DivergenceFaceMatrixStencil) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Divergence div(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_east = 2.;
    FaceMatrixStencil s_f;
    for (std::size_t n{0}; n < N; n++) {
        s_f.SetValues(Positions::WEST, n, v_west * n, -v_west * n);
        s_f.SetValues(Positions::EAST, n, v_east * n, -v_east * n);
    }

    auto s_c = div(s_f);
    static_assert(std::is_same_v<decltype(s_c), CenterMatrixStencil>, "Type is wrong!");

    VecSC A = grid_rep.GetFaceArea();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s_c.GetValue(Positions::WEST, n), -v_west * n * A[0]);
        EXPECT_EQ(s_c.GetValue(Positions::EAST, n), v_east * n * A[0]);
        EXPECT_EQ(s_c.GetValue(Positions::CENTER, n), (v_west - v_east) * n * A[0]);
    }
}

TEST_F(IntegrationTestCartesianOperators2D, DivergenceFaceMatrixStencil) {
    GridType::Options opt{0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Divergence div(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_east = 2.;
    double v_south = 4.;
    double v_north = 6.;
    FaceMatrixStencil s_f;
    for (std::size_t n{0}; n < N; n++) {
        s_f.SetValues(Positions::WEST, n, v_west * n, -v_west * n);
        s_f.SetValues(Positions::EAST, n, v_east * n, -v_east * n);
        s_f.SetValues(Positions::SOUTH, n, v_south * n, -v_south * n);
        s_f.SetValues(Positions::NORTH, n, v_north * n, -v_north * n);
    }

    auto s_c = div(s_f);
    static_assert(std::is_same_v<decltype(s_c), CenterMatrixStencil>, "Type is wrong!");

    VecSC A = grid_rep.GetFaceArea();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_NEAR(s_c.GetValue(Positions::WEST, n), -v_west * n * A[0],
                    std::numeric_limits<SC>::epsilon() * std::abs(v_west * A[0] *n));
        EXPECT_NEAR(s_c.GetValue(Positions::EAST, n), v_east * n * A[0],
                    std::numeric_limits<SC>::epsilon() * std::abs(v_east * (A[0] * n)));
        EXPECT_NEAR(s_c.GetValue(Positions::SOUTH, n), -v_south * n * A[1],
                    std::numeric_limits<SC>::epsilon() * std::abs(v_south * A[0] * n));
        EXPECT_NEAR(s_c.GetValue(Positions::NORTH, n), v_north * n * A[1],
                    std::numeric_limits<SC>::epsilon() * std::abs(v_north * A[0] * n));
        EXPECT_NEAR(s_c.GetValue(Positions::CENTER, n),
                (v_west - v_east) * n * A[0] + (v_south - v_north) * n * A[1],
                v_north * std::numeric_limits<SC>::epsilon() * std::abs(A[0]));
    }
}

TEST_F(IntegrationTestCartesianOperators3D, DivergenceFaceMatrixStencil) {
    GridType::Options opt{0, 0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;

    Divergence div(grid_rep, ordinal_internal);

    double v_west = 0.;
    double v_east = 2.;
    double v_south = 4.;
    double v_north = 6.;
    double v_bottom = 8.;
    double v_top = 10.;
    FaceMatrixStencil s_f;
    for (std::size_t n{0}; n < N; n++) {
        s_f.SetValues(Positions::WEST, n, v_west * n, -v_west * n);
        s_f.SetValues(Positions::EAST, n, v_east * n, -v_east * n);
        s_f.SetValues(Positions::SOUTH, n, v_south * n, -v_south * n);
        s_f.SetValues(Positions::NORTH, n, v_north * n, -v_north * n);
        s_f.SetValues(Positions::BOTTOM, n, v_bottom * n, -v_bottom * n);
        s_f.SetValues(Positions::TOP, n, v_top * n, -v_top * n);
    }

    auto s_c = div(s_f);
    static_assert(std::is_same_v<decltype(s_c), CenterMatrixStencil>, "Type is wrong!");

    VecSC A = grid_rep.GetFaceArea();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(s_c.GetValue(Positions::WEST, n), -v_west * n * A[0]);
        EXPECT_EQ(s_c.GetValue(Positions::EAST, n), v_east * n * A[0]);
        EXPECT_EQ(s_c.GetValue(Positions::SOUTH, n), -v_south * n * A[1]);
        EXPECT_EQ(s_c.GetValue(Positions::NORTH, n), v_north * n * A[1]);
        EXPECT_EQ(s_c.GetValue(Positions::BOTTOM, n), -v_bottom * n * A[2]);
        EXPECT_EQ(s_c.GetValue(Positions::TOP, n), v_top * n * A[2]);
        EXPECT_NEAR(s_c.GetValue(Positions::CENTER, n), (v_west - v_east) * n * A[0]
                                                    + (v_south - v_north) * n * A[1]
                                                    + (v_bottom - v_top) * n * A[2], 1e-14);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, MatrixBlockIntegration) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(3, 3, 3);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    Divergence div(grid_rep, ordinal_internal);
    Gradient grad(grid_rep, ordinal_internal);

    mb = -1. * div(grad(mb));

    VecSC A = grid_rep.GetFaceArea();
    VecSC dn_r;
    for (auto& e : dn_r)
        e = 1.;
    dn_r /= grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), -A[0] * dn_r[0]);
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), -A[0] * dn_r[0]);
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), 2.* A[0] * dn_r[0]);
    }
}

TEST_F(IntegrationTestCartesianOperators2D, MatrixBlockIntegration) {
    GridType::Options opt{0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(5, 5, 5);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    Divergence div(grid_rep, ordinal_internal);
    Gradient grad(grid_rep, ordinal_internal);

    mb = -1. * div(grad(mb));

    VecSC A = grid_rep.GetFaceArea();
    VecSC dn_r;
    for (auto& e : dn_r)
        e = 1.;
    dn_r /= grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), -A[0] * dn_r[0]);
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), -A[0] * dn_r[0]);
        EXPECT_EQ(mb.Get(n, n, Positions::SOUTH), -A[1] * dn_r[1]);
        EXPECT_EQ(mb.Get(n, n, Positions::NORTH), -A[1] * dn_r[1]);
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), 2. * A[0] * dn_r[0] + 2. * A[1] * dn_r[1]);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, MatrixBlockSet) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(3, 3, 3);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    SC v_center = 2., v_west = -1., v_east = -2., v_rhs = 1.5;
    CenterMatrixStencil1C comp;
    for (std::size_t n{0}; n < N; n++) {
        comp.SetValue(Positions::CENTER, 0, v_center + n);
        comp.SetValue(Positions::WEST, 0, v_west + n);
        comp.SetValue(Positions::EAST, 0, v_east + n);
        comp.SetRHS(0, v_rhs + n);
        mb.Set(n, comp);
    }

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), v_west + n);
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), v_east + n);
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), v_center + n);
        EXPECT_EQ(mb.GetRhs(n), v_rhs + n);
    }
}

TEST_F(IntegrationTestCartesianOperators2D, MatrixBlockSet) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(3, 3, 3);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    SC v_center = 2., v_west = -1., v_east = -2., v_south = -3., v_north = -4., v_rhs = 1.5;
    CenterMatrixStencil1C comp;
    for (std::size_t n{0}; n < N; n++) {
        comp.SetValue(Positions::CENTER, 0, v_center + n);
        comp.SetValue(Positions::WEST, 0, v_west + n);
        comp.SetValue(Positions::EAST, 0, v_east + n);
        comp.SetValue(Positions::SOUTH, 0, v_south + n);
        comp.SetValue(Positions::NORTH, 0, v_north + n);
        comp.SetRHS(0, v_rhs + n);
        mb.Set(n, comp);
    }

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), v_west + n);
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), v_east + n);
        EXPECT_EQ(mb.Get(n, n, Positions::SOUTH), v_south + n);
        EXPECT_EQ(mb.Get(n, n, Positions::NORTH), v_north + n);
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), v_center + n);
        EXPECT_EQ(mb.GetRhs(n), v_rhs + n);
    }
}

TEST_F(IntegrationTestCartesianOperators3D, MatrixBlockSet) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(3, 3, 3);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    SC v_center = 2., v_west = -1., v_east = -2.,
       v_south = -3., v_north = -4.,
       v_bottom = -5., v_top = -6., v_rhs = 1.5;
    CenterMatrixStencil1C comp;
    for (std::size_t n{0}; n < N; n++) {
        comp.SetValue(Positions::CENTER, 0, v_center + n);
        comp.SetValue(Positions::WEST, 0, v_west + n);
        comp.SetValue(Positions::EAST, 0, v_east + n);
        comp.SetValue(Positions::SOUTH, 0, v_south + n);
        comp.SetValue(Positions::NORTH, 0, v_north + n);
        comp.SetValue(Positions::BOTTOM, 0, v_bottom + n);
        comp.SetValue(Positions::TOP, 0, v_top + n);
        comp.SetRHS(0, v_rhs + n);
        mb.Set(n, comp);
    }

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), v_west + n);
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), v_east + n);
        EXPECT_EQ(mb.Get(n, n, Positions::SOUTH), v_south + n);
        EXPECT_EQ(mb.Get(n, n, Positions::NORTH), v_north + n);
        EXPECT_EQ(mb.Get(n, n, Positions::BOTTOM), v_bottom + n);
        EXPECT_EQ(mb.Get(n, n, Positions::TOP), v_top + n);
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), v_center + n);
        EXPECT_EQ(mb.GetRhs(n), v_rhs + n);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, MatrixBlockAdd) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(3, 3, 3);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    SC v_center = 2., v_west = -1., v_east = -2., v_rhs = 1.5;
    CenterMatrixStencil1C comp;
    for (std::size_t n{0}; n < N; n++) {
        comp.SetValue(Positions::CENTER, 0, v_center + n);
        comp.SetValue(Positions::WEST, 0, v_west + n);
        comp.SetValue(Positions::EAST, 0, v_east + n);
        comp.SetRHS(0, v_rhs + n);
        mb.Set(n, comp);
        mb.Add(n, comp);
    }

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), 2. * (v_west + n));
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), 2. * (v_east + n));
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), 2. * (v_center + n));
        EXPECT_EQ(mb.GetRhs(n), 2.*(v_rhs + n));
    }
}

TEST_F(IntegrationTestCartesianOperators2D, MatrixBlockAdd) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(3, 3, 3);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    SC v_center = 2., v_west = -1., v_east = -2.,
       v_south = -3., v_north = -4., v_rhs = 1.5;
    CenterMatrixStencil1C comp;
    for (std::size_t n{0}; n < N; n++) {
        comp.SetValue(Positions::CENTER, 0, v_center + n);
        comp.SetValue(Positions::WEST, 0, v_west + n);
        comp.SetValue(Positions::EAST, 0, v_east + n);
        comp.SetValue(Positions::SOUTH, 0, v_south + n);
        comp.SetValue(Positions::NORTH, 0, v_north + n);
        comp.SetRHS(0, v_rhs + n);
        mb.Set(n, comp);
        mb.Add(n, comp);
    }

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), 2. * (v_west + n));
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), 2. * (v_east + n));
        EXPECT_EQ(mb.Get(n, n, Positions::SOUTH), 2. * (v_south + n));
        EXPECT_EQ(mb.Get(n, n, Positions::NORTH), 2. * (v_north + n));
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), 2. * (v_center + n));
        EXPECT_EQ(mb.GetRhs(n), 2. * (v_rhs + n));
    }
}

TEST_F(IntegrationTestCartesianOperators3D, MatrixBlockAdd) {
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(3, 3, 3);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    SC v_center = 2., v_west = -1., v_east = -2.,
       v_south = -3., v_north = -4.,
       v_bottom = -5., v_top = -6., v_rhs = 1.5;
    CenterMatrixStencil1C comp;
    for (std::size_t n{0}; n < N; n++) {
        comp.SetValue(Positions::CENTER, 0, v_center + n);
        comp.SetValue(Positions::WEST, 0, v_west + n);
        comp.SetValue(Positions::EAST, 0, v_east + n);
        comp.SetValue(Positions::SOUTH, 0, v_south + n);
        comp.SetValue(Positions::NORTH, 0, v_north + n);
        comp.SetValue(Positions::BOTTOM, 0, v_bottom + n);
        comp.SetValue(Positions::TOP, 0, v_top + n);
        comp.SetRHS(0, v_rhs + n);
        mb.Set(n, comp);
        mb.Add(n, comp);
    }

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), 2. * (v_west + n));
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), 2. * (v_east + n));
        EXPECT_EQ(mb.Get(n, n, Positions::SOUTH), 2. * (v_south + n));
        EXPECT_EQ(mb.Get(n, n, Positions::NORTH), 2. * (v_north + n));
        EXPECT_EQ(mb.Get(n, n, Positions::BOTTOM), 2. * (v_bottom + n));
        EXPECT_EQ(mb.Get(n, n, Positions::TOP), 2. * (v_top + n));
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), 2. * (v_center + n));
        EXPECT_EQ(mb.GetRhs(n), 2. * (v_rhs + n));
    }
}

TEST_F(IntegrationTestCartesianOperators3D, MatrixBlockIntegration) {
    GridType::Options opt{0, 0, 0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    LO ordinal_internal = 0;
    dare::utils::Vector<N, std::size_t> size_hint(7, 7, 7);
    MatrixBlock mb(&grid_rep, ordinal_internal, size_hint);
    Divergence div(grid_rep, ordinal_internal);
    Gradient grad(grid_rep, ordinal_internal);

    mb = -1. * div(grad(mb));

    VecSC A = grid_rep.GetFaceArea();
    VecSC dn_r;
    for (auto& e : dn_r)
        e = 1.;
    dn_r /= grid_rep.GetDistances();

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mb.Get(n, n, Positions::WEST), -A[0] * dn_r[0]);
        EXPECT_EQ(mb.Get(n, n, Positions::EAST), -A[0] * dn_r[0]);
        EXPECT_EQ(mb.Get(n, n, Positions::SOUTH), -A[1] * dn_r[1]);
        EXPECT_EQ(mb.Get(n, n, Positions::NORTH), -A[1] * dn_r[1]);
        EXPECT_EQ(mb.Get(n, n, Positions::BOTTOM), -A[2] * dn_r[2]);
        EXPECT_EQ(mb.Get(n, n, Positions::TOP), -A[2] * dn_r[2]);
        EXPECT_EQ(mb.Get(n, n, Positions::CENTER), 2. * A[0] * dn_r[0] + 2. * A[1] * dn_r[1] + 2. * A[2] * dn_r[2]);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDCDSValueTest) {
    using Scheme = dare::Matrix::CDS;
    using CNB = dare::Grid::CartesianNeighbor;
    const SC eps_tol{1e3};
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * n;
        }
    }

    auto v = tvd.Interpolate(field);
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceValueStencil<GridType, SC, N>>);

    Index ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = 0.5 * (field.At(ind, n) + field.At(ind_l, n));
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        ind_l.i() += 2;
        v_ex = 0.5 * (field.At(ind, n) + field.At(ind_l, n));
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }

    // Downwind velocity
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);
    v = tvd_down.Interpolate(field);

    ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = 0.5 * (field.At(ind, n) + field.At(ind_l, n));
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        ind_l.i() += 2;
        v_ex = 0.5 * (field.At(ind, n) + field.At(ind_l, n));
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDCDSMatrixTest) {
    using Scheme = dare::Matrix::CDS;
    using CNB = dare::Grid::CartesianNeighbor;
    const SC eps_tol{1e3};
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * n;
        }
    }

    auto v = tvd * field;
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceMatrixStencil<GridType, SC, N>>);

    Index ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), velocity[0]);

        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = -0.5 * (field.At(ind, n) - field.At(ind_l, n)) * velocity[0];
        EXPECT_NEAR(v.GetRHS(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        ind_l.i() += 2;
        v_ex = -0.5 * (field.At(ind_l, n) - field.At(ind, n)) * velocity[0];
        EXPECT_NEAR(v.GetRHS(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // Downwind velocity
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);
    v = tvd_down * field;

    ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), 0.);

        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = -0.5 * (field.At(ind_l, n) - field.At(ind, n)) * velocity[0];
        EXPECT_NEAR(v.GetRHS(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        ind_l.i() += 2;
        v_ex = -0.5 * (field.At(ind, n) - field.At(ind_l, n)) * velocity[0];
        EXPECT_NEAR(v.GetRHS(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDUPWINDValueTest) {
    using Scheme = dare::Matrix::UPWIND;
    using CNB = dare::Grid::CartesianNeighbor;
    const SC eps_tol{1e3};
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * n;
        }
    }

    auto v = tvd.Interpolate(field);
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceValueStencil<GridType, SC, N>>);

    Index ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = field.At(ind_l, n);
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        ind_l.i() += 1;
        v_ex = field.At(ind_l, n);
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }

    // down wind
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);

    v = tvd_down.Interpolate(field);

    for (std::size_t n{0}; n < N; n++) {
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        ind_l.i() += 1;
        v_ex = field.At(ind_l, n);
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDUPWINDMatrixTest) {
    using Scheme = dare::Matrix::UPWIND;
    using CNB = dare::Grid::CartesianNeighbor;
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * n;
        }
    }

    auto v = tvd * field;
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceMatrixStencil<GridType, SC, N>>);

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), velocity[0]);

        SC v_ex = 0.;
        EXPECT_EQ(v.GetRHS(CNB::WEST, n), v_ex);
        EXPECT_EQ(v.GetRHS(CNB::EAST, n), v_ex);
    }

    // Downwind velocity
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);
    v = tvd_down * field;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), 0.);

        SC v_ex = 0.;
        EXPECT_EQ(v.GetRHS(CNB::WEST, n), v_ex);
        EXPECT_EQ(v.GetRHS(CNB::EAST, n), v_ex);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDMINMODValueTest) {
    using Scheme = dare::Matrix::MINMOD;
    using CNB = dare::Grid::CartesianNeighbor;
    const SC eps_tol{1e3};
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * i * n;
        }
    }

    auto v = tvd.Interpolate(field);
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceValueStencil<GridType, SC, N>>);

    Index ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        Index ind_w{ind};
        Index ind_ww{ind};
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_ww = field.At(ind_ww, n);
        SC r_w = (phi_w - phi_ww) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = std::max(0., std::min(r_w, 1.));
        SC v_ex = phi_w + 0.5 * r_w * (phi_c - phi_w);
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        Index ind_e{ind};
        ind_e.i() += 1;
        SC phi_e = field.At(ind_e, n);
        SC r_e = (phi_c - phi_w) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = std::max(0., std::min(r_e, 1.));
        v_ex = phi_c + 0.5 * r_e * (phi_e - phi_c);
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }

    // down wind
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);

    v = tvd_down.Interpolate(field);

    for (std::size_t n{0}; n < N; n++) {
        Index ind_w{ind};
        Index ind_e{ind};
        ind_w.i() -= 1;
        ind_e.i() += 1;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_e = field.At(ind_e, n);
        SC r_w = (phi_e - phi_c) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = std::max(0., std::min(r_w, 1.));
        SC v_ex = phi_c + 0.5 * r_w * (phi_w - phi_c);
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        Index ind_ee{ind};
        ind_ee.i() += 2;
        SC phi_ee = field.At(ind_ee, n);
        SC r_e = (phi_ee - phi_e) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = std::max(0., std::min(r_e, 1.));
        v_ex = phi_e + 0.5 * r_e * (phi_c - phi_e);
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDMINMODMatrixTest) {
    using Scheme = dare::Matrix::MINMOD;
    using CNB = dare::Grid::CartesianNeighbor;
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * i * n;
        }
    }

    auto v = tvd * field;
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceMatrixStencil<GridType, SC, N>>);

    Index ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), velocity[0]);

        Index ind_w{ind};
        Index ind_ww{ind};
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_ww = field.At(ind_ww, n);
        SC r_w = (phi_w - phi_ww) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = std::max(0., std::min(r_w, 1.));
        SC v_ex_w = -0.5 * r_w * (phi_c - phi_w) * velocity[0];
        Index ind_e{ind};
        ind_e.i() += 1;
        SC phi_e = field.At(ind_e, n);
        SC r_e = (phi_c - phi_w) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = std::max(0., std::min(r_e, 1.));
        SC v_ex_e = -0.5 * r_e * (phi_e - phi_c) * velocity[0];
        EXPECT_EQ(v.GetRHS(CNB::WEST, n), v_ex_w);
        EXPECT_EQ(v.GetRHS(CNB::EAST, n), v_ex_e);
    }

    // Downwind velocity
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);
    v = tvd_down * field;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), 0.);

        Index ind_w{ind};
        Index ind_e{ind};
        ind_w.i() -= 1;
        ind_e.i() += 1;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_e = field.At(ind_e, n);
        SC r_w = (phi_e - phi_c) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = std::max(0., std::min(r_w, 1.));
        SC v_ex_w = -0.5 * r_w * (phi_w - phi_c) * velocity[0];
        Index ind_ee{ind};
        ind_ee.i() += 2;
        SC phi_ee = field.At(ind_ee, n);
        SC r_e = (phi_ee - phi_e) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = std::max(0., std::min(r_e, 1.));
        SC v_ex_e = -0.5 * r_e * (phi_c - phi_e) * velocity[0];
        EXPECT_EQ(v.GetRHS(CNB::WEST, n), v_ex_w);
        EXPECT_EQ(v.GetRHS(CNB::EAST, n), v_ex_e);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDVANALBADAValueTest) {
    using Scheme = dare::Matrix::VANALBADA;
    using CNB = dare::Grid::CartesianNeighbor;
    const SC eps_tol{1e3};
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * n;
        }
    }

    auto v = tvd.Interpolate(field);
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceValueStencil<GridType, SC, N>>);

    Index ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        Index ind_w{ind};
        Index ind_ww{ind};
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_ww = field.At(ind_ww, n);
        SC r_w = (phi_w - phi_ww) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = (r_w + r_w * r_w) / (1. + r_w * r_w);
        SC v_ex = phi_w + 0.5 * r_w * (phi_c - phi_w);
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        Index ind_e{ind};
        ind_e.i() += 1;
        SC phi_e = field.At(ind_e, n);
        SC r_e = (phi_c - phi_w) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = (r_e + r_e * r_e) / (1. + r_e * r_e);
        v_ex = phi_c + 0.5 * r_e * (phi_e - phi_c);
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }

    // down wind
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);

    v = tvd_down.Interpolate(field);

    for (std::size_t n{0}; n < N; n++) {
        Index ind_w{ind};
        Index ind_e{ind};
        ind_w.i() -= 1;
        ind_e.i() += 1;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_e = field.At(ind_e, n);
        SC r_w = (phi_e - phi_c) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = (r_w + r_w * r_w) / (1. + r_w * r_w);
        SC v_ex = phi_c + 0.5 * r_w * (phi_w - phi_c);
        EXPECT_NEAR(v.GetValue(CNB::WEST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
        Index ind_ee{ind};
        ind_ee.i() += 2;
        SC phi_ee = field.At(ind_ee, n);
        SC r_e = (phi_ee - phi_e) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = (r_e + r_e * r_e) / (1. + r_e * r_e);
        v_ex = phi_e + 0.5 * r_e * (phi_c - phi_e);
        EXPECT_NEAR(v.GetValue(CNB::EAST, n), v_ex, eps_tol * std::numeric_limits<SC>::epsilon() * v_ex);
    }
}

TEST_F(IntegrationTestCartesianOperators1D, TVDVANALBADAMatrixTest) {
    using Scheme = dare::Matrix::VANALBADA;
    using CNB = dare::Grid::CartesianNeighbor;
    GridType::Options opt{0};  // not staggered
    auto grid_rep = grid->GetRepresentation(opt);
    VecSC velocity{1.};  // upwind velocity
    LO internal_ordinal{0};
    TVD<Scheme> tvd(grid_rep, internal_ordinal, velocity);

    Field field("values", grid_rep);
    for (LO i{0}; i < grid_rep.GetLocalResolution().i(); i++) {
        Index ind(i);
        for (std::size_t n{0}; n < N; n++) {
            field.At(ind, n) = i * i * n;
        }
    }

    auto v = tvd * field;
    static_assert(std::is_same_v<decltype(v), dare::Data::FaceMatrixStencil<GridType, SC, N>>);

    Index ind = grid_rep.MapOrdinalToIndexLocal(grid_rep.MapInternalToLocal(internal_ordinal));
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), velocity[0]);

        Index ind_w{ind};
        Index ind_ww{ind};
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_ww = field.At(ind_ww, n);
        SC r_w = (phi_w - phi_ww) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = (r_w + r_w * r_w) / (1. + r_w * r_w);
        SC v_ex_w = -0.5 * r_w * (phi_c - phi_w) * velocity[0];
        Index ind_e{ind};
        ind_e.i() += 1;
        SC phi_e = field.At(ind_e, n);
        SC r_e = (phi_c - phi_w) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = (r_e + r_e * r_e) / (1. + r_e * r_e);
        SC v_ex_e = -0.5 * r_e * (phi_e - phi_c) * velocity[0];
        EXPECT_EQ(v.GetRHS(CNB::WEST, n), v_ex_w);
        EXPECT_EQ(v.GetRHS(CNB::EAST, n), v_ex_e);
    }

    // Downwind velocity
    velocity[0] = -1.;
    TVD<Scheme> tvd_down(grid_rep, internal_ordinal, velocity);
    v = tvd_down * field;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(v.GetValueNeighbor(CNB::WEST, n), 0.);
        EXPECT_EQ(v.GetValueCenter(CNB::WEST, n), velocity[0]);
        EXPECT_EQ(v.GetValueNeighbor(CNB::EAST, n), velocity[0]);
        EXPECT_EQ(v.GetValueCenter(CNB::EAST, n), 0.);

        Index ind_w{ind};
        Index ind_e{ind};
        ind_w.i() -= 1;
        ind_e.i() += 1;
        SC phi_c = field.At(ind, n);
        SC phi_w = field.At(ind_w, n);
        SC phi_e = field.At(ind_e, n);
        SC r_w = (phi_e - phi_c) / (phi_c - phi_w);
        if (std::isnan(r_w))
            r_w = 0.;
        r_w = (r_w + r_w * r_w) / (1. + r_w * r_w);
        SC v_ex_w = -0.5 * r_w * (phi_w - phi_c) * velocity[0];
        Index ind_ee{ind};
        ind_ee.i() += 2;
        SC phi_ee = field.At(ind_ee, n);
        SC r_e = (phi_ee - phi_e) / (phi_e - phi_c);
        if (std::isnan(r_e))
            r_e = 0.;
        r_e = (r_e + r_e * r_e) / (1. + r_e * r_e);
        SC v_ex_e = -0.5 * r_e * (phi_c - phi_e) * velocity[0];
        EXPECT_EQ(v.GetRHS(CNB::WEST, n), v_ex_w);
        EXPECT_EQ(v.GetRHS(CNB::EAST, n), v_ex_e);
    }
}
