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
    static const std::size_t N{3};
    using GridType = dare::Grid::Cartesian<Dim>;
    using SC = double;
    using CenterMatrixStencil = dare::Data::CenterMatrixStencil<GridType, SC, N>;
    using CenterValueStencil = dare::Data::CenterValueStencil<GridType, SC, N>;
    using FaceMatrixStencil = dare::Data::FaceMatrixStencil<GridType, SC, N>;
    using FaceValueStencil = dare::Data::FaceValueStencil<GridType, SC, N>;

    void SetUp() {
    }
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
        for (std::size_t n{0}; n < N; n++){
            stencil.SetValue(Positions::CENTER, n, c_center + n);
            stencil.SetValue(Positions::WEST, n, c_west + n);
            stencil.SetValue(Positions::EAST, n, c_east + n);
        }
    };
    Reset();

    /*
     *  Operations with scalars
     */
    // Multiplication with scalar
    stencil_compare = 2. * stencil;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * 2.);
    }
    Reset();

        // Division by scalar
    stencil_compare = stencil / 2.;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / 2.);
    }
    Reset();

    /*
     *  Operations with stencils
     */
    CenterMatrixStencil stencil_op;
    double op_center = 0.1;
    double op_west = 0.2;
    double op_east = 0.3;
    for (std::size_t n{0}; n < N; n++) {
        stencil_op.SetValue(Positions::CENTER, n, op_center + n);
        stencil_op.SetValue(Positions::WEST, n, op_west + n);
        stencil_op.SetValue(Positions::EAST, n, op_east + n);
    }

    // Addition of stencil
    stencil_compare = stencil + stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) + (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) + (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) + (op_east + n));
    }
    Reset();

    // Subtraction of stencil
    stencil_compare = stencil - stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) - (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) - (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) - (op_east + n));
    }
    Reset();

    // Multiplication with stencil
    stencil_compare = stencil * stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) * (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * (op_east + n));
    }
    Reset();

    // Division by stencil
    stencil_compare = stencil / stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) / (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / (op_east + n));
    }
    Reset();

    // Set all values
    stencil.SetAll(1.);
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil.GetValue(Positions::CENTER, n), 1.);
        EXPECT_EQ(stencil.GetValue(Positions::WEST, n), 1.);
        EXPECT_EQ(stencil.GetValue(Positions::EAST, n), 1.);
    }
    Reset();
}

TEST_F(IntegrationTestCartesianStencils2D, CenterMatrixBasicOperations) {
    using Positions = CenterMatrixStencil::Positions;
    CenterMatrixStencil stencil, stencil_compare;
    double c_center = 1.;
    double c_west = -1.;
    double c_east = -2.;
    double c_south = -3.;
    double c_north = -4.;

    auto Reset = [&]() {
        for (std::size_t n{0}; n < N; n++) {
            stencil.SetValue(Positions::CENTER, n, c_center + n);
            stencil.SetValue(Positions::WEST, n, c_west + n);
            stencil.SetValue(Positions::EAST, n, c_east + n);
            stencil.SetValue(Positions::SOUTH, n, c_south + n);
            stencil.SetValue(Positions::NORTH, n, c_north + n);
        }
    };
    Reset();

    /*
     *  Operations with scalars
     */
    // Multiplication with scalar
    stencil_compare = 2. * stencil;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) * 2.);
    }
    Reset();

    // Division by scalar
    stencil_compare = stencil / 2.;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) / 2.);
    }
    Reset();

    /*
     *  Operations with stencils
     */
    CenterMatrixStencil stencil_op;
    double op_center = 0.1;
    double op_west = 0.2;
    double op_east = 0.3;
    double op_south = 0.4;
    double op_north = 0.5;

    for (std::size_t n{0}; n < N; n++) {
        stencil_op.SetValue(Positions::CENTER, n, op_center + n);
        stencil_op.SetValue(Positions::WEST, n, op_west + n);
        stencil_op.SetValue(Positions::EAST, n, op_east + n);
        stencil_op.SetValue(Positions::SOUTH, n, op_south + n);
        stencil_op.SetValue(Positions::NORTH, n, op_north + n);
    }

    // Addition of stencil
    stencil_compare = stencil + stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) + (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) + (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) + (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) + (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) + (op_north + n));
    }
    Reset();

    // Subtraction of stencil
    stencil_compare = stencil - stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) - (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) - (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) - (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) - (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) - (op_north + n));
    }
    Reset();

    // Multiplication with stencil
    stencil_compare = stencil * stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) * (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) * (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) * (op_north + n));
    }
    Reset();

    // Division by stencil
    stencil_compare = stencil / stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) / (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) / (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) / (op_north + n));
    }
    Reset();
}

TEST_F(IntegrationTestCartesianStencils3D, CenterMatrixBasicOperations) {
    using Positions = CenterMatrixStencil::Positions;
    CenterMatrixStencil stencil, stencil_compare;
    double c_center = 1.;
    double c_west = -1.;
    double c_east = -2.;
    double c_south = -3.;
    double c_north = -4.;
    double c_bottom = -5.;
    double c_top = -6.;
    auto Reset = [&]() {
        for (std::size_t n{0}; n < N; n++) {
            stencil.SetValue(Positions::CENTER, n, c_center + n);
            stencil.SetValue(Positions::WEST, n, c_west + n);
            stencil.SetValue(Positions::EAST, n, c_east + n);
            stencil.SetValue(Positions::SOUTH, n, c_south + n);
            stencil.SetValue(Positions::NORTH, n, c_north + n);
            stencil.SetValue(Positions::BOTTOM, n, c_bottom + n);
            stencil.SetValue(Positions::TOP, n, c_top + n);
        }
    };
    Reset();

    /*
     *  Operations with scalars
     */
    // Multiplication with scalar
    stencil_compare = 2. * stencil;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) * 2.);
        }
    Reset();

    // Division by scalar
    stencil_compare = stencil / 2.;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) / 2.);
    }
    Reset();

    /*
     *  Operations with stencils
     */
    CenterMatrixStencil stencil_op;
    double op_center = 0.1;
    double op_west = 0.2;
    double op_east = 0.3;
    double op_south = 0.4;
    double op_north = 0.5;
    double op_bottom = 0.6;
    double op_top = 0.7;
    for (std::size_t n{0}; n < N; n++) {
        stencil_op.SetValue(Positions::CENTER, n, op_center + n);
        stencil_op.SetValue(Positions::WEST, n, op_west + n);
        stencil_op.SetValue(Positions::EAST, n, op_east + n);
        stencil_op.SetValue(Positions::SOUTH, n, op_south + n);
        stencil_op.SetValue(Positions::NORTH, n, op_north + n);
        stencil_op.SetValue(Positions::BOTTOM, n, op_bottom + n);
        stencil_op.SetValue(Positions::TOP, n, op_top + n);
    }

    // Addition of stencil
    stencil_compare = stencil + stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) + (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) + (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) + (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) + (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) + (op_north + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) + (op_bottom + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) + (op_top + n));
    }
    Reset();

    // Subtraction of stencil
    stencil_compare = stencil - stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) - (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) - (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) - (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) - (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) - (op_north + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) - (op_bottom + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) - (op_top + n));
    }
    Reset();

    // Multiplication with stencil
    stencil_compare = stencil * stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) * (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) * (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) * (op_north + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) * (op_bottom + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) * (op_top + n));
    }
    Reset();

    // Division by stencil
    stencil_compare = stencil / stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) / (op_center + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / (op_east + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) / (op_south + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) / (op_north + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) / (op_bottom + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) / (op_top + n));
    }
    Reset();
}

TEST_F(IntegrationTestCartesianStencils1D, CenterValueInteroperability) {
    using Positions = CenterValueStencil::Positions;
    CenterValueStencil stencil, stencil_compare;
    double c_center = 1.;
    double c_west = -1.;
    double c_east = -2.;
    auto Reset = [&]() {
        for (std::size_t n{0}; n < N; n++) {
            stencil.SetValue(Positions::CENTER, n, c_center + n);
            stencil.SetValue(Positions::WEST, n, c_west + n);
            stencil.SetValue(Positions::EAST, n, c_east + n);
        }
    };
    Reset();

    CenterMatrixStencil matrix_stencil;
    for (std::size_t n{0}; n < N; n++) {
        matrix_stencil.SetValue(Positions::CENTER, n, 1.);
        matrix_stencil.SetValue(Positions::WEST, n, 1.);
        matrix_stencil.SetValue(Positions::EAST, n, 1.);
    }

    stencil_compare = matrix_stencil + stencil;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) + 1.);
    }
    Reset();

    stencil_compare = stencil - matrix_stencil;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) - 1.);
    }
}

TEST_F(IntegrationTestCartesianStencils2D, CenterValueInteroperability) {
    using Positions = CenterValueStencil::Positions;
    CenterValueStencil stencil, stencil_compare;
    double c_center = 1.;
    double c_west = -1.;
    double c_east = -2.;
    double c_south = -3.;
    double c_north = -4.;
    auto Reset = [&]() {
        for (std::size_t n{0}; n < N; n++) {
            stencil.SetValue(Positions::CENTER, n, c_center + n);
            stencil.SetValue(Positions::WEST, n, c_west + n);
            stencil.SetValue(Positions::EAST, n, c_east + n);
            stencil.SetValue(Positions::SOUTH, n, c_south + n);
            stencil.SetValue(Positions::NORTH, n, c_north + n);
        }
    };
    Reset();

    CenterMatrixStencil matrix_stencil;
    for (std::size_t n{0}; n < N; n++) {
        matrix_stencil.SetValue(Positions::CENTER, n, 1.);
        matrix_stencil.SetValue(Positions::WEST, n, 1.);
        matrix_stencil.SetValue(Positions::EAST, n, 1.);
        matrix_stencil.SetValue(Positions::SOUTH, n, 1.);
        matrix_stencil.SetValue(Positions::NORTH, n, 1.);
    }

    stencil_compare = matrix_stencil + stencil;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) + 1.);
    }
    Reset();

    stencil_compare = stencil - matrix_stencil;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) - 1.);
    }
}

TEST_F(IntegrationTestCartesianStencils3D, CenterValueInteroperability) {
    using Positions = CenterValueStencil::Positions;
    CenterValueStencil stencil, stencil_compare;
    double c_center = 1.;
    double c_west = -1.;
    double c_east = -2.;
    double c_south = -3.;
    double c_north = -4.;
    double c_bottom = -5.;
    double c_top = -6.;

    auto Reset = [&]() {
        for (std::size_t n{0}; n < N; n++) {
            stencil.SetValue(Positions::CENTER, n, c_center + n);
            stencil.SetValue(Positions::WEST, n, c_west + n);
            stencil.SetValue(Positions::EAST, n, c_east + n);
            stencil.SetValue(Positions::SOUTH, n, c_south + n);
            stencil.SetValue(Positions::NORTH, n, c_north + n);
            stencil.SetValue(Positions::BOTTOM, n, c_bottom + n);
            stencil.SetValue(Positions::TOP, n, c_top + n);
        }
    };
    Reset();

    CenterMatrixStencil matrix_stencil;
    for (std::size_t n{0}; n < N; n++) {
        matrix_stencil.SetValue(Positions::CENTER, n, 1.);
        matrix_stencil.SetValue(Positions::WEST, n, 1.);
        matrix_stencil.SetValue(Positions::EAST, n, 1.);
        matrix_stencil.SetValue(Positions::SOUTH, n, 1.);
        matrix_stencil.SetValue(Positions::NORTH, n, 1.);
        matrix_stencil.SetValue(Positions::BOTTOM, n, 1.);
        matrix_stencil.SetValue(Positions::TOP, n, 1.);
    }

    stencil_compare = matrix_stencil + stencil;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) + 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) + 1.);
    }
    Reset();

    stencil_compare = stencil - matrix_stencil;

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::CENTER, n), (c_center + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::SOUTH, n), (c_south + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::NORTH, n), (c_north + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::BOTTOM, n), (c_bottom + n) - 1.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::TOP, n), (c_top + n) - 1.);
    }
}

TEST_F(IntegrationTestCartesianStencils1D, FaceMatrixBasicOperations) {
    using Positions = FaceMatrixStencil::Positions;
    FaceMatrixStencil stencil, stencil_compare;
    double c_center = 1.;
    double c_west = -1.;
    double c_east = -2.;
    auto Reset = [&]() {
        for (std::size_t n{0}; n < N; n++) {
            stencil.SetValues(Positions::WEST, n, c_west + n, c_center + c_west + n);
            stencil.SetValues(Positions::EAST, n, c_east + n, c_center + c_east + n);
        }
    };
    Reset();

    /*
     *  Operations with scalars
     */
    // Multiplication with scalar
    stencil_compare = 2. * stencil;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::WEST, n), (c_west + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::WEST, n), (c_center + c_west + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::EAST, n), (c_east + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::EAST, n), (c_center + c_east + n) * 2.);
    }
    Reset();

    // Division by scalar
    stencil_compare = stencil / 2.;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::WEST, n), (c_west + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::WEST, n), (c_center + c_west + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::EAST, n), (c_east + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::EAST, n), (c_center + c_east + n) / 2.);
    }
    Reset();

    /*
     *  Operations with stencils
     */
    FaceMatrixStencil stencil_op;
    double op_center = 0.1;
    double op_west = 0.2;
    double op_east = 0.3;
    for (std::size_t n{0}; n < N; n++) {
        stencil_op.SetValues(Positions::WEST, n, op_west + n, op_center + op_west + n);
        stencil_op.SetValues(Positions::EAST, n, op_east + n, op_center + op_east + n);
    }

    // Addition of stencil
    stencil_compare = stencil + stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::WEST, n), (c_west + n) + (op_west + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::WEST, n),
                    (c_center + c_west + n) + (op_center + op_west + n));
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::EAST, n), (c_east + n) + (op_east + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::EAST, n),
                    (c_center + c_east + n) + (op_center + op_east + n));
    }
    Reset();

    // Subtraction of stencil
    stencil_compare = stencil - stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::WEST, n), (c_west + n) - (op_west + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::WEST, n),
                  (c_center + c_west + n) - (op_center + op_west + n));
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::EAST, n), (c_east + n) - (op_east + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::EAST, n),
                  (c_center + c_east + n) - (op_center + op_east + n));
    }
    Reset();

    // Multiplication with stencil
    stencil_compare = stencil * stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::WEST, n), (c_west + n) * (op_west + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::WEST, n),
                  (c_center + c_west + n) * (op_center + op_west + n));
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::EAST, n), (c_east + n) * (op_east + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::EAST, n),
                  (c_center + c_east + n) * (op_center + op_east + n));
    }
    Reset();

    // Division by stencil
    stencil_compare = stencil / stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::WEST, n), (c_west + n) / (op_west + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::WEST, n),
                  (c_center + c_west + n) / (op_center + op_west + n));
        EXPECT_EQ(stencil_compare.GetValueNeighbor(Positions::EAST, n), (c_east + n) / (op_east + n));
        EXPECT_EQ(stencil_compare.GetValueCenter(Positions::EAST, n),
                  (c_center + c_east + n) / (op_center + op_east + n));
    }
    Reset();

    // Set all values
    stencil.SetAll(1.);
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil.GetValueNeighbor(Positions::WEST, n), 1.);
        EXPECT_EQ(stencil.GetValueCenter(Positions::WEST, n), 1.);
        EXPECT_EQ(stencil.GetValueNeighbor(Positions::EAST, n), 1.);
        EXPECT_EQ(stencil.GetValueCenter(Positions::EAST, n), 1.);
    }
    Reset();
}

TEST_F(IntegrationTestCartesianStencils1D, FaceValueBasicOperations) {
    using Positions = FaceValueStencil::Positions;
    FaceValueStencil stencil, stencil_compare;
    double c_west = -1.;
    double c_east = -2.;
    auto Reset = [&]() {
        for (std::size_t n{0}; n < N; n++) {
            stencil.SetValue(Positions::WEST, n, c_west + n);
            stencil.SetValue(Positions::EAST, n, c_east + n);
        }
    };
    Reset();

    /*
     *  Operations with scalars
     */
    // Multiplication with scalar
    stencil_compare = 2. * stencil;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * 2.);
    }
    Reset();

    // Division by scalar
    stencil_compare = stencil / 2.;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / 2.);
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / 2.);
    }
    Reset();

    /*
     *  Operations with stencils
     */
    FaceValueStencil stencil_op;
    double op_west = 0.2;
    double op_east = 0.3;
    for (std::size_t n{0}; n < N; n++) {
        stencil_op.SetValue(Positions::WEST, n, op_west + n);
        stencil_op.SetValue(Positions::EAST, n, op_east + n);
    }

    // Addition of stencil
    stencil_compare = stencil + stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) + (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) + (op_east + n));
    }
    Reset();

    // Subtraction of stencil
    stencil_compare = stencil - stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) - (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) - (op_east + n));
    }
    Reset();

    // Multiplication with stencil
    stencil_compare = stencil * stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) * (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) * (op_east + n));
    }
    Reset();

    // Division by stencil
    stencil_compare = stencil / stencil_op;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_compare.GetValue(Positions::WEST, n), (c_west + n) / (op_west + n));
        EXPECT_EQ(stencil_compare.GetValue(Positions::EAST, n), (c_east + n) / (op_east + n));
    }
    Reset();

    // Set all values
    stencil.SetAll(1.);
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil.GetValue(Positions::WEST, n), 1.);
        EXPECT_EQ(stencil.GetValue(Positions::EAST, n), 1.);
    }
    Reset();

    // Test operations with with FaceMatrixStencil
    FaceMatrixStencil stencil_f, stencil_f_compare;
    double cf_center = 1.;
    double cf_west = -1.;
    double cf_east = -2.;
    for (std::size_t n{0}; n < N; n++) {
        stencil_f.SetValues(Positions::WEST, n, cf_west + n, cf_center + cf_west + n);
        stencil_f.SetValues(Positions::EAST, n, cf_east + n, cf_center + cf_east + n);
    }

    // Multiplication with FaceMatrixStencil
    stencil_f_compare = stencil * stencil_f;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(stencil_f_compare.GetValueNeighbor(Positions::WEST, n), (cf_west + n) * (c_west + n));
        EXPECT_EQ(stencil_f_compare.GetValueCenter(Positions::WEST, n),
                  (cf_center + cf_west + n) * (c_west + n));
        EXPECT_EQ(stencil_f_compare.GetValueNeighbor(Positions::EAST, n), (cf_east + n) * (c_east + n));
        EXPECT_EQ(stencil_f_compare.GetValueCenter(Positions::EAST, n),
                  (cf_center + cf_east + n) * (c_east + n));
    }
    Reset();
}

// TODO(@Dave): implement remaining cases for 2D and 3D
