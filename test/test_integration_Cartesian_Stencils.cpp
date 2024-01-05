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

#include "../MatrixSystem/Stencils_Cartesian.h"

template <std::size_t Dim>
class IntegrationTestCartesianStencils : public testing::Test {
public:
    using GridType = dare::Grid::Cartesian<Dim>;
    using CenterMatrixStencil = dare::Data::CenterMatrixStencil<GridType>;

    void SetUp() {
        // const LO num_ghost{2};
        // grid = std::make_unique<GridType>(&exec_man,
        //                                   dare::test::GetResolutionTestCartesianGrid<Dim, GO>(),
        //                                   dare::test::GetSizeTestCartesianGrid<Dim, SC>(),
        //                                   num_ghost);
    }

    // std::unique_ptr<GridType> grid;
    // dare::mpi::ExecutionManager exec_man;
};

using IntegrationTestCartesianStencils1D = IntegrationTestCartesianStencils<1>;
using IntegrationTestCartesianStencils2D = IntegrationTestCartesianStencils<2>;
using IntegrationTestCartesianStencils3D = IntegrationTestCartesianStencils<3>;

TEST_F(IntegrationTestCartesianStencils1D, CenterMatrixBasicOperations) {
    using Positions = CenterMatrixStencil::Positions;
    CenterMatrixStencil stencil, stencil_compare;
    double c_center = 1.;
    double c_west = -1.;
    double c_east = -2.;
    auto Reset = [&]() {
        stencil.SetValue(Positions::CENTER, c_center);
        stencil.SetValue(Positions::WEST, c_west);
        stencil.SetValue(Positions::EAST, c_east);
    };
    Reset();

    /*
     *  Operations with scalars
     */
    // Multiplication with scalar
    stencil_compare = 2. * stencil;
    EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER), c_center * 2.);
    EXPECT_EQ(stencil_compare.GetValue(Positions::WEST), c_west * 2.);
    EXPECT_EQ(stencil_compare.GetValue(Positions::EAST), c_east * 2.);
    Reset();

    // Division by scalar
    stencil_compare = stencil / 2.;
    EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER), c_center / 2.);
    EXPECT_EQ(stencil_compare.GetValue(Positions::WEST), c_west / 2.);
    EXPECT_EQ(stencil_compare.GetValue(Positions::EAST), c_east / 2.);
    Reset();

    /*
     *  Operations with stencils
     */
    CenterMatrixStencil stencil_op;
    double op_center = 0.1;
    double op_west = 0.2;
    double op_east = 0.3;
    stencil_op.SetValue(Positions::CENTER, op_center);
    stencil_op.SetValue(Positions::WEST, op_west);
    stencil_op.SetValue(Positions::EAST, op_east);

    // Addition of stencil
    stencil_compare = stencil + stencil_op;
    EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER), c_center + op_center);
    EXPECT_EQ(stencil_compare.GetValue(Positions::WEST), c_west + op_west);
    EXPECT_EQ(stencil_compare.GetValue(Positions::EAST), c_east + op_east);
    Reset();

    // Subtraction of stencil
    stencil_compare = stencil - stencil_op;
    EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER), c_center - op_center);
    EXPECT_EQ(stencil_compare.GetValue(Positions::WEST), c_west - op_west);
    EXPECT_EQ(stencil_compare.GetValue(Positions::EAST), c_east - op_east);
    Reset();

    // Multiplication with stencil
    stencil_compare = stencil * stencil_op;
    EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER), c_center * op_center);
    EXPECT_EQ(stencil_compare.GetValue(Positions::WEST), c_west * op_west);
    EXPECT_EQ(stencil_compare.GetValue(Positions::EAST), c_east * op_east);
    Reset();

    // Division by stencil
    stencil_compare = stencil / stencil_op;
    EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER), c_center / op_center);
    EXPECT_EQ(stencil_compare.GetValue(Positions::WEST), c_west / op_west);
    EXPECT_EQ(stencil_compare.GetValue(Positions::EAST), c_east / op_east);
    Reset();
}

TEST_F(IntegrationTestCartesianStencils2D, CenterMatrixBasicOperations) {
    CenterMatrixStencil stencil;
}

TEST_F(IntegrationTestCartesianStencils3D, CenterMatrixBasicOperations) {
    CenterMatrixStencil stencil;
}
