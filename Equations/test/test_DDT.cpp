/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder

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

#include <limits>
#include <memory>
#include <random>

#include "Grid/Cartesian.h"  // for whatever reason this one needs to be loaded in first

#include "Equations/DDT.h"
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

/*!
 * @brief fixture for testing time partial derivative with Cartesian grid
 */
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
    const std::size_t NUM_TFIELD{2};
    GridType::Options opt(0, 0, 0);
    GridType::Representation grep = grid->GetRepresentation(opt);
    FieldType field_phi("phi", grep, NUM_TFIELD);
    FieldType field_p1("property_1", grep, NUM_TFIELD);
    FieldType field_p2("property_2", grep, NUM_TFIELD);
    SC phi_base{0.3};
    SC p1_base{-0.2};
    SC p2_base{0.16};

    for (std::size_t t{0}; t < NUM_TFIELD; t++) {
        SC f{static_cast<SC>(t) + 1.};
        field_phi.SetValues(f * phi_base, t);
        field_p1.SetValues(f * p1_base, t);
        field_p2.SetValues(f * p2_base, t);
    }

    LO ordinal = 0;
    SC dt = 0.1;
    dare::Matrix::DDT<GridType> ddt(grep, ordinal, dt);
    SC dV_dt = grep.GetCellVolume(ordinal) / dt;
    auto s_1 = ddt(field_phi);
    auto s_2 = ddt(field_p1, field_phi);
    auto s_3 = ddt(field_p1, field_p2, field_phi);

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_NEAR(s_1.Center(n), dV_dt, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(s_1.GetRHS(n), dV_dt * 2. * phi_base, std::numeric_limits<SC>::epsilon());
    }

    for (std::size_t n{0}; n < N; n++) {
        EXPECT_NEAR(s_2.Center(n), p1_base  * dV_dt, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(s_2.GetRHS(n), dV_dt * p1_base * 4. * phi_base, std::numeric_limits<SC>::epsilon());
    }
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_NEAR(s_3.Center(n), p2_base * p1_base * dV_dt, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(s_3.GetRHS(n), dV_dt * p2_base * p1_base * 8. * phi_base, std::numeric_limits<SC>::epsilon());
    }
}

TEST_F(DDTTest, EulerBackwardStaggered) {
    const std::size_t NUM_TFIELD{2};
    GridType::Index ind(2, 3, 4);
    GridType::Index ind_nb(1, 3, 4);
    GridType::Options opt_s(0, 0, 0);
    GridType::Options opt_x(1, 0, 0);
    GridType::Representation grep_s = grid->GetRepresentation(opt_s);
    GridType::Representation grep_x = grid->GetRepresentation(opt_x);

    FieldType field_phi("phi", grep_x, NUM_TFIELD);
    FieldType field_p1("property_1", grep_s, NUM_TFIELD);
    FieldType field_p2("property_2", grep_s, NUM_TFIELD);

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * distribution(generator); };

    for (std::size_t t{0}; t < NUM_TFIELD; t++) {
        for (std::size_t n{0}; n < field_phi.GetDataVector(t).GetSize(); n++)
            field_phi.GetDataVector(t).At(n) = GetRandValue();
        for (std::size_t n{0}; n < field_p1.GetDataVector(t).GetSize(); n++)
            field_p1.GetDataVector(t).At(n) = GetRandValue();
        for (std::size_t n{0}; n < field_p2.GetDataVector(t).GetSize(); n++)
            field_p2.GetDataVector(t).At(n) = GetRandValue();
    }

    LO ordinal = grep_x.MapIndexToOrdinalLocal(ind);
    LO ordinal_internal = grep_x.MapIndexToOrdinalLocalInternal(grep_x.MapLocalToInternal(ind));
    SC dt = 0.1;
    dare::Matrix::DDT<GridType> ddt(grep_x, ordinal_internal, dt);
    SC dV_dt = grep_x.GetCellVolume(ordinal) / dt;
    auto s_1 = ddt(field_phi);
    auto s_2 = ddt(field_p1, field_phi);
    auto s_3 = ddt(field_p1, field_p2, field_phi);

    for (std::size_t n{0}; n < N; n++) {
        SC v_0 = s_1.Center(n);
        SC v_1 = s_1.GetRHS(n);
        SC v_0_ex = dV_dt;
        SC v_1_ex = dV_dt * field_phi.GetDataVector(1).At(ind, n);
        EXPECT_NEAR(v_0, v_0_ex, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(v_1, v_1_ex, std::numeric_limits<SC>::epsilon());
    }

    for (std::size_t n{0}; n < N; n++) {
        SC p1_0 = 0.5 * (field_p1.GetDataVector().At(ind, n) + field_p1.GetDataVector().At(ind_nb, n));
        SC p1_1 = 0.5 * (field_p1.GetDataVector(1).At(ind, n) + field_p1.GetDataVector(1).At(ind_nb, n));
        SC v_0 = s_2.Center(n);
        SC v_1 = s_2.GetRHS(n);
        SC v_0_ex = p1_0 * dV_dt;
        SC v_1_ex = p1_1 * dV_dt * field_phi.GetDataVector(1).At(ind, n);
        EXPECT_NEAR(v_0, v_0_ex, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(v_1, v_1_ex, std::numeric_limits<SC>::epsilon());
    }

    for (std::size_t n{0}; n < N; n++) {
        SC p1_0 = 0.5 * (field_p1.GetDataVector().At(ind, n) + field_p1.GetDataVector().At(ind_nb, n));
        SC p1_1 = 0.5 * (field_p1.GetDataVector(1).At(ind, n) + field_p1.GetDataVector(1).At(ind_nb, n));
        SC p2_0 = 0.5 * (field_p2.GetDataVector().At(ind, n) + field_p2.GetDataVector().At(ind_nb, n));
        SC p2_1 = 0.5 * (field_p2.GetDataVector(1).At(ind, n) + field_p2.GetDataVector(1).At(ind_nb, n));
        SC v_0 = s_3.Center(n);
        SC v_1 = s_3.GetRHS(n);
        SC v_0_ex = p2_0 * p1_0 * dV_dt;
        SC v_1_ex = p2_1 * p1_1 * dV_dt * field_phi.GetDataVector(1).At(ind, n);
        EXPECT_NEAR(v_0, v_0_ex, std::numeric_limits<SC>::epsilon());
        EXPECT_NEAR(v_1, v_1_ex, std::numeric_limits<SC>::epsilon());
    }
}
