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

#include "Equations/DDT.h"
#include "Grid/Cartesian.h"

namespace dare::test {

template <std::size_t Dim, typename GO>
dare::utils::Vector<Dim, GO> GetResolutionTestDDT() {
    dare::utils::Vector<Dim, GO> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 10 + n;
    return res;
}
template <std::size_t Dim, typename SC>
dare::utils::Vector<Dim, SC> GetSizeTestDDT() {
    dare::utils::Vector<Dim, SC> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

}  // namespace dare::test

class DDTTest : public testing::Test {
public:
    static const std::size_t N{3};
    using GridType = dare::Grid::Cartesian<3>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using Index = typename GridType::Index;
    using CenterMatrixStencil = dare::Data::CenterMatrixStencil<GridType, SC, N>;
    using CenterValueStencil = dare::Data::CenterValueStencil<GridType, SC, N>;
    using FieldType = dare::Data::Field<GridType, SC, N>;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(&exec_man,
                                          dare::test::GetResolutionTestDDT<3, GO>(),
                                          dare::test::GetSizeTestDDT<3, SC>(),
                                          num_ghost);
    }

    std::unique_ptr<GridType> grid;        //!< the grid
    dare::mpi::ExecutionManager exec_man;  //!< the execution manager
};

TEST_F(DDTTest, EulerBackward) {
    using TimeDiscretization = dare::Matrix::EULER_BACKWARD;
    const std::size_t NUM_TFIELD{TimeDiscretization::NUM_TIMESTEPS + 1};
    GridType::Options opt(0, 0, 0);
    GridType::Representation grep = grid->GetRepresentation(opt);
    SC phi_base{0.3};
    SC p1_base{-0.2};
    SC p2_base{0.16};
    FieldType field_phi("phi", grep, NUM_TFIELD);
    FieldType field_p1("property_1", grep, NUM_TFIELD);
    FieldType field_p2("property_2", grep, NUM_TFIELD);

    for (std::size_t t{0}; t < NUM_TFIELD; t++) {
        SC f{static_cast<SC>(t) + 1.};
        field_phi.SetValues(f * phi_base, t);
        field_p1.SetValues(f * p1_base, t);
        field_p2.SetValues(f * p2_base, t);
    }
    LO ordinal = 0;
    SC dt = 0.1;
    dare::Matrix::DDT<GridType, TimeDiscretization> ddt(grep, ordinal, dt);
    SC dV_dt = grep.GetCellVolume(ordinal) / dt;
    auto s_1 = ddt(field_phi);
    auto s_2 = ddt(field_p1, field_phi);
    auto s_3 = ddt(field_p1, field_p2, field_phi);

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_NEAR(s_1.Center(n), dV_dt , std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(s_1.GetRHS(n), 0., std::numeric_limits<SC>::epsilon());
    }
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_NEAR(s_2.Center(n), p1_base  * dV_dt, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(s_2.GetRHS(n), 0., std::numeric_limits<SC>::epsilon());
    }
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_NEAR(s_3.Center(n), p2_base * p1_base * dV_dt, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(s_3.GetRHS(n), 0., std::numeric_limits<SC>::epsilon());
    }
}
