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

#include <random>

#include "Grid/DefaultTypes.h"
#include "Math/Interpolation_Cartesian.h"

namespace dare::test {

template <std::size_t Dim>
dare::utils::Vector<Dim, defaults::GlobalOrdinalType> GetResolutionTestCartesianInterpolation() {
    dare::utils::Vector < Dim, defaults::GlobalOrdinalType> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 10 + n;
    return res;
}
template <std::size_t Dim>
dare::utils::Vector<Dim, defaults::ScalarType> GetSizeTestCartesianInterpolation() {
    dare::utils::Vector<Dim, defaults::ScalarType> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

}  // namespace dare::test

class IntegrationTestCartesianInterpolation : public testing::Test {
public:
    static const std::size_t Dim{3};
    static const std::size_t N{Dim};
    using GridType = dare::Grid::Cartesian<Dim>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using Index = typename GridType::Index;
    using Field = dare::Data::GridVector<GridType, SC, N>;
    using Options = typename GridType::Options;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(&exec_man,
                                          dare::test::GetResolutionTestCartesianInterpolation<Dim>(),
                                          dare::test::GetSizeTestCartesianInterpolation<Dim>(),
                                          num_ghost);
    }

    std::unique_ptr<GridType> grid;        //!< the grid
    dare::mpi::ExecutionManager exec_man;  //!< the execution manager
};

TEST_F(IntegrationTestCartesianInterpolation, GetInterpolationIndices0DimTest) {
    Index ind{2, 3, 4};
    Options off_rel{0, 0, 0};

    // vectors with affected dimensions
    dare::utils::Vector<0, std::size_t> dim_aff;

    // for testing indices
    dare::utils::Vector<1, bool> found_ind;

    // zero affected dimension (value at center)
    auto vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                         off_rel,
                                                                         dim_aff);
    static_assert(std::is_same_v< decltype(vlist), dare::utils::Vector<1, Index> >);
    for (std::size_t n{0}; n < Dim; n++)
        EXPECT_EQ(vlist[0][n], ind[n]);
}

TEST_F(IntegrationTestCartesianInterpolation, GetInterpolationIndices1DimTest) {
    Index ind{2, 3, 4};
    Options off_rel{0, 0, 0};

    // vectors with affected dimensions
    dare::utils::Vector<1, std::size_t> dim_aff;

    // for testing indices
    dare::utils::Vector<2, bool> found_ind;

    // WEST
    off_rel.SetAllValues(0);
    off_rel[0] = -1;
    dim_aff[0] = 0;
    auto vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                         off_rel,
                                                                         dim_aff);
    static_assert(std::is_same_v<decltype(vlist), dare::utils::Vector<2, Index> >);
    dare::utils::Vector<2, Index> vlist_ex;  // expected values
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() -= 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 2; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 2; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // EAST
    off_rel.SetAllValues(0);
    off_rel[0] = 1;
    dim_aff[0] = 0;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() += 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 2; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 2; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // SOUTH
    off_rel.SetAllValues(0);
    off_rel[1] = -1;
    dim_aff[0] = 1;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].j() -= 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 2; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 2; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // NORTH
    off_rel.SetAllValues(0);
    off_rel[1] = 1;
    dim_aff[0] = 1;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].j() += 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 2; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 2; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // BOTTOM
    off_rel.SetAllValues(0);
    off_rel[2] = -1;
    dim_aff[0] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].k() -= 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 2; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 2; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // TOP
    off_rel.SetAllValues(0);
    off_rel[2] = 1;
    dim_aff[0] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].k() += 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 2; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 2; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }
}

TEST_F(IntegrationTestCartesianInterpolation, GetInterpolationIndices2DimTest) {
    Index ind{2, 3, 4};
    Options off_rel{0, 0, 0};

    // vectors with affected dimensions
    dare::utils::Vector<2, std::size_t> dim_aff;

    // for testing indices
    dare::utils::Vector<4, bool> found_ind;

    // SOUTH-WEST
    off_rel.SetAllValues(0);
    off_rel[0] = -1;
    off_rel[1] = -1;
    dim_aff[0] = 0;
    dim_aff[1] = 1;
    auto vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                         off_rel,
                                                                         dim_aff);
    static_assert(std::is_same_v<decltype(vlist), dare::utils::Vector<4, Index> >);
    dare::utils::Vector<4, Index> vlist_ex;  // expected values
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() -= 1;
    vlist_ex[1].i() -= 1;
    vlist_ex[1].j() -= 1;
    vlist_ex[2].j() -= 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // NORTH-WEST
    off_rel.SetAllValues(0);
    off_rel[0] = -1;
    off_rel[1] = 1;
    dim_aff[0] = 0;
    dim_aff[1] = 1;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);

    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() -= 1;
    vlist_ex[1].i() -= 1;
    vlist_ex[1].j() += 1;
    vlist_ex[2].j() += 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // NORTH-EAST
    off_rel.SetAllValues(0);
    off_rel[0] = 1;
    off_rel[1] = 1;
    dim_aff[0] = 1;  // just some wild mixing
    dim_aff[1] = 0;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);

    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() += 1;
    vlist_ex[1].i() += 1;
    vlist_ex[1].j() += 1;
    vlist_ex[2].j() += 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // SOUTH-EAST
    off_rel.SetAllValues(0);
    off_rel[0] = 1;
    off_rel[1] = -1;
    dim_aff[0] = 1;  // just some wild mixing
    dim_aff[1] = 0;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);

    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() += 1;
    vlist_ex[1].i() += 1;
    vlist_ex[1].j() -= 1;
    vlist_ex[2].j() -= 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // BOTTOM-WEST
    off_rel.SetAllValues(0);
    off_rel[0] = -1;
    off_rel[2] = -1;
    dim_aff[0] = 0;
    dim_aff[1] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);

    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() -= 1;
    vlist_ex[1].i() -= 1;
    vlist_ex[1].k() -= 1;
    vlist_ex[2].k() -= 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // TOP-EAST
    off_rel.SetAllValues(0);
    off_rel[0] = 1;
    off_rel[2] = 1;
    dim_aff[0] = 0;
    dim_aff[1] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);

    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() += 1;
    vlist_ex[1].i() += 1;
    vlist_ex[1].k() += 1;
    vlist_ex[2].k() += 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // BOTTOM-SOUTH
    off_rel.SetAllValues(0);
    off_rel[1] = -1;
    off_rel[2] = -1;
    dim_aff[0] = 1;
    dim_aff[1] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);

    vlist_ex.SetAllValues(ind);
    vlist_ex[0].j() -= 1;
    vlist_ex[1].j() -= 1;
    vlist_ex[1].k() -= 1;
    vlist_ex[2].k() -= 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // TOP-NORTH
    off_rel.SetAllValues(0);
    off_rel[1] = 1;
    off_rel[2] = 1;
    dim_aff[0] = 1;
    dim_aff[1] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);

    vlist_ex.SetAllValues(ind);
    vlist_ex[0].j() += 1;
    vlist_ex[1].j() += 1;
    vlist_ex[1].k() += 1;
    vlist_ex[2].k() += 1;
    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 4; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 4; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }
}

TEST_F(IntegrationTestCartesianInterpolation, GetInterpolationIndices3DimTest) {
    Index ind{2, 3, 4};
    Options off_rel{0, 0, 0};

    // vectors with affected dimensions
    dare::utils::Vector<3, std::size_t> dim_aff;

    // for testing indices
    dare::utils::Vector<8, bool> found_ind;

    // BOTTOM-SOUTH-WEST
    off_rel.SetAllValues(0);
    off_rel[0] = -1;
    off_rel[1] = -1;
    off_rel[2] = -1;
    dim_aff[0] = 0;
    dim_aff[1] = 1;
    dim_aff[2] = 2;
    auto vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                         off_rel,
                                                                         dim_aff);
    static_assert(std::is_same_v<decltype(vlist), dare::utils::Vector<8, Index> >);
    dare::utils::Vector<8, Index> vlist_ex;  // expected values
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() -= 1;
    vlist_ex[1].i() -= 1;
    vlist_ex[1].j() -= 1;
    vlist_ex[2].j() -= 1;
    vlist_ex[3].i() -= 1;
    vlist_ex[3].k() -= 1;
    vlist_ex[4].i() -= 1;
    vlist_ex[4].j() -= 1;
    vlist_ex[4].k() -= 1;
    vlist_ex[5].j() -= 1;
    vlist_ex[5].k() -= 1;
    vlist_ex[6].k() -= 1;

    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 8; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 8; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // TOP-NORTH-EAST
    off_rel.SetAllValues(0);
    off_rel[0] = 1;
    off_rel[1] = 1;
    off_rel[2] = 1;
    dim_aff[0] = 0;
    dim_aff[1] = 1;
    dim_aff[2] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() += 1;
    vlist_ex[1].i() += 1;
    vlist_ex[1].j() += 1;
    vlist_ex[2].j() += 1;
    vlist_ex[3].i() += 1;
    vlist_ex[3].k() += 1;
    vlist_ex[4].i() += 1;
    vlist_ex[4].j() += 1;
    vlist_ex[4].k() += 1;
    vlist_ex[5].j() += 1;
    vlist_ex[5].k() += 1;
    vlist_ex[6].k() += 1;

    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 8; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 8; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }

    // BOTTOM-NORTH-EAST
    off_rel.SetAllValues(0);
    off_rel[0] = 1;
    off_rel[1] = 1;
    off_rel[2] = -1;
    dim_aff[0] = 0;
    dim_aff[1] = 1;
    dim_aff[2] = 2;
    vlist = dare::math::details::Cartesian::GetInterpolationIndices(ind,
                                                                    off_rel,
                                                                    dim_aff);
    vlist_ex.SetAllValues(ind);
    vlist_ex[0].i() += 1;
    vlist_ex[1].i() += 1;
    vlist_ex[1].j() += 1;
    vlist_ex[2].j() += 1;
    vlist_ex[3].i() += 1;
    vlist_ex[3].k() -= 1;
    vlist_ex[4].i() += 1;
    vlist_ex[4].j() += 1;
    vlist_ex[4].k() -= 1;
    vlist_ex[5].j() += 1;
    vlist_ex[5].k() -= 1;
    vlist_ex[6].k() -= 1;

    found_ind.SetAllValues(false);
    for (std::size_t n{0}; n < 8; n++) {
        for (const auto v : vlist)
            found_ind[n] |= (v == vlist_ex[n]);
    }
    for (std::size_t n{0}; n < 8; n++) {
        EXPECT_TRUE(found_ind[n]) << "Did not find the expected Index: " << vlist_ex[n];
    }
}

TEST_F(IntegrationTestCartesianInterpolation, Interpolation1DimLinearTest) {
    Options off_rel{0, 0, 0};
    Options opt{0, 0, 0};           // not staggered, in this test it's either way irrelevant
    dare::utils::Vector<1, std::size_t> dim_aff{0};
    auto grep = grid->GetRepresentation(opt);
    Field field("grid", grep);
    Index ind(2, 3, 4);
    dare::utils::Vector<Dim, SC> grad;
    grad[0] = 1.;
    grad[1] = 1.5;
    grad[2] = -0.3;
    grad /= grep.GetDistances();

    for (LO i{0}; i < grep.GetLocalResolution().i(); i++) {
        for (LO j{0}; j < grep.GetLocalResolution().j(); j++) {
            for (LO k{0}; k < grep.GetLocalResolution().k(); k++) {
                field.At(Index(i, j, k), 0) = i * grad[0];
                field.At(Index(i, j, k), 1) = j * grad[1];
                field.At(Index(i, j, k), 2) = k * grad[2];
            }
        }
    }

    // WEST
    std::size_t n_comp = 0;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = -1;
    dim_aff[0] = n_comp;

    SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n_comp);
    SC v_ex = (ind[n_comp] * grad[n_comp] + (ind[n_comp] + off_rel[n_comp]) * grad[n_comp]) * 0.5;
    EXPECT_EQ(v, v_ex);

    // EAST
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = 1;
    dim_aff[0] = n_comp;

    v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n_comp);
    v_ex = (ind[n_comp] * grad[n_comp] + (ind[n_comp] + off_rel[n_comp]) * grad[n_comp]) * 0.5;
    EXPECT_EQ(v, v_ex);

    // SOUTH
    n_comp = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = -1;
    dim_aff[0] = n_comp;

    v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n_comp);
    v_ex = (ind[n_comp] * grad[n_comp] + (ind[n_comp] + off_rel[n_comp]) * grad[n_comp]) * 0.5;
    EXPECT_EQ(v, v_ex);

    // NORTH
    n_comp = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = 1;
    dim_aff[0] = n_comp;

    v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n_comp);
    v_ex = (ind[n_comp] * grad[n_comp] + (ind[n_comp] + off_rel[n_comp]) * grad[n_comp]) * 0.5;
    EXPECT_EQ(v, v_ex);

    // BOTTOM
    n_comp = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = -1;
    dim_aff[0] = n_comp;

    v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n_comp);
    v_ex = (ind[n_comp] * grad[n_comp] + (ind[n_comp] + off_rel[n_comp]) * grad[n_comp]) * 0.5;
    EXPECT_EQ(v, v_ex);

    // TOP
    n_comp = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = 1;
    dim_aff[0] = n_comp;

    v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n_comp);
    v_ex = (ind[n_comp] * grad[n_comp] + (ind[n_comp] + off_rel[n_comp]) * grad[n_comp]) * 0.5;
    EXPECT_EQ(v, v_ex);

    // Multiple components
    // WEST
    n_comp = 0;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = -1;
    dim_aff[n_comp] = 0;

    dare::utils::Vector<N, SC> v_m
        = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    dare::utils::Vector<N, SC> v_m_ex;
    for (std::size_t n{0}; n < N; n++) {
        v_m_ex[n] = (ind[n] * grad[n] + (ind[n] + off_rel[n]) * grad[n]) * 0.5;
    }
    EXPECT_EQ(v_m, v_m_ex);

    // Multiple components
    // EAST
    n_comp = 0;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = 1;
    dim_aff[0] = n_comp;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        v_m_ex[n] = (ind[n] * grad[n] + (ind[n] + off_rel[n]) * grad[n]) * 0.5;
    }
    EXPECT_EQ(v_m, v_m_ex);

    // Multiple components
    // SOUTH
    n_comp = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = -1;
    dim_aff[0] = n_comp;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    
    for (std::size_t n{0}; n < N; n++) {
        v_m_ex[n] = (ind[n] * grad[n] + (ind[n] + off_rel[n]) * grad[n]) * 0.5;
    }
    EXPECT_EQ(v_m, v_m_ex);

    // NORTH
    n_comp = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = 1;
    dim_aff[0] = n_comp;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    
    for (std::size_t n{0}; n < N; n++) {
        v_m_ex[n] = (ind[n] * grad[n] + (ind[n] + off_rel[n]) * grad[n]) * 0.5;
    }
    EXPECT_EQ(v_m, v_m_ex);

    // BOTTOM
    n_comp = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = -1;
    dim_aff[0] = n_comp;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    
    for (std::size_t n{0}; n < N; n++) {
        v_m_ex[n] = (ind[n] * grad[n] + (ind[n] + off_rel[n]) * grad[n]) * 0.5;
    }
    EXPECT_EQ(v_m, v_m_ex);

    // TOP
    n_comp = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[n_comp] = 1;
    dim_aff[0] = n_comp;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    
    for (std::size_t n{0}; n < N; n++) {
        v_m_ex[n] = (ind[n] * grad[n] + (ind[n] + off_rel[n]) * grad[n]) * 0.5;
    }
    EXPECT_EQ(v_m, v_m_ex);
}

TEST_F(IntegrationTestCartesianInterpolation, Interpolation2DimLinearTest) {
    Options off_rel{0, 0, 0};
    Options opt{0, 0, 0};  // not staggered, in this test it's either way irrelevant
    dare::utils::Vector<2, std::size_t> dim_aff{0};
    auto grep = grid->GetRepresentation(opt);
    Field field("grid", grep);
    Index ind(2, 3, 4);
    const double tol_eps{10.};

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * distribution(generator); };

    for (LO i{0}; i < grep.GetLocalResolution().i(); i++) {
        for (LO j{0}; j < grep.GetLocalResolution().j(); j++) {
            for (LO k{0}; k < grep.GetLocalResolution().k(); k++) {
                field.At(Index(i, j, k), 0) = GetRandValue();
                field.At(Index(i, j, k), 1) = GetRandValue();
                field.At(Index(i, j, k), 2) = GetRandValue();
            }
        }
    }

    // SOUTH - WEST
    std::size_t dir_1 = 0;
    std::size_t dir_2 = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = -1;
    off_rel[dir_2] = -1;
    dim_aff[0] = dir_1;
    dim_aff[1] = dir_2;

    auto v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    static_assert(std::is_same_v<decltype(v_m), dare::utils::Vector<N, SC>>);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // NORTH - WEST
    dir_1 = 0;
    dir_2 = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = -1;
    off_rel[dir_2] = 1;
    dim_aff[0] = dir_1;
    dim_aff[1] = dir_2;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // NORTH - EAST
    dir_1 = 0;
    dir_2 = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = 1;
    off_rel[dir_2] = 1;
    dim_aff[0] = dir_1;
    dim_aff[1] = dir_2;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // SOUTH - EAST
    dir_1 = 0;
    dir_2 = 1;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = 1;
    off_rel[dir_2] = -1;
    dim_aff[0] = dir_1;
    dim_aff[1] = dir_2;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // BOTTOM - WEST
    dir_1 = 0;
    dir_2 = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = -1;
    off_rel[dir_2] = -1;
    dim_aff[0] = dir_1;
    dim_aff[1] = dir_2;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP - NORTH
    dir_1 = 0;
    dir_2 = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = 1;
    off_rel[dir_2] = 1;
    dim_aff[0] = dir_1;
    dim_aff[1] = dir_2;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }
}

TEST_F(IntegrationTestCartesianInterpolation, Interpolation3DimLinearTest) {
    Options off_rel{0, 0, 0};
    Options opt{0, 0, 0};  // not staggered, in this test it's either way irrelevant
    dare::utils::Vector<3, std::size_t> dim_aff{0, 0, 0};
    auto grep = grid->GetRepresentation(opt);
    Field field("grid", grep);
    Index ind(2, 3, 4);
    const double tol_eps{10.};

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * distribution(generator); };

    for (LO i{0}; i < grep.GetLocalResolution().i(); i++) {
        for (LO j{0}; j < grep.GetLocalResolution().j(); j++) {
            for (LO k{0}; k < grep.GetLocalResolution().k(); k++) {
                for (std::size_t n{0}; n < N; n++) {
                    field.At(Index(i, j, k), n) = GetRandValue();
                }
            }
        }
    }

    // BOTTOM - SOUTH - WEST
    std::size_t dir_1 = 0;
    std::size_t dir_2 = 1;
    std::size_t dir_3 = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = -1;
    off_rel[dir_2] = -1;
    off_rel[dir_3] = -1;
    dim_aff[0] = dir_2;  // some wild mixing
    dim_aff[1] = dir_3;
    dim_aff[2] = dir_1;

    auto v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    static_assert(std::is_same_v<decltype(v_m), dare::utils::Vector<N, SC>>);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l = ind;
        ind_l[dir_3] += off_rel[dir_3];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.125;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // BOTTOM - NORTH - EAST
    dir_1 = 0;
    dir_2 = 1;
    dir_3 = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = 1;
    off_rel[dir_2] = 1;
    off_rel[dir_3] = -1;
    dim_aff[0] = dir_2;  // some wild mixing
    dim_aff[1] = dir_3;
    dim_aff[2] = dir_1;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l = ind;
        ind_l[dir_3] += off_rel[dir_3];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.125;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP - NORTH - WEST
    dir_1 = 0;
    dir_2 = 1;
    dir_3 = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = -1;
    off_rel[dir_2] = 1;
    off_rel[dir_3] = 1;
    dim_aff[0] = dir_2;  // some wild mixing
    dim_aff[1] = dir_3;
    dim_aff[2] = dir_1;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l = ind;
        ind_l[dir_3] += off_rel[dir_3];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.125;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP - SOUTH - EAST
    dir_1 = 0;
    dir_2 = 1;
    dir_3 = 2;
    off_rel.SetAllValues(0);
    dim_aff.SetAllValues(0);
    off_rel[dir_1] = 1;
    off_rel[dir_2] = -1;
    off_rel[dir_3] = 1;
    dim_aff[0] = dir_2;  // some wild mixing
    dim_aff[1] = dir_3;
    dim_aff[2] = dir_1;

    v_m = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::details::Cartesian::InterpolateCartesianLinear(ind, field, off_rel, dim_aff, n);
        Index ind_l{ind};
        SC v_ex = field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l = ind;
        ind_l[dir_3] += off_rel[dir_3];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] += off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        ind_l[dir_2] += off_rel[dir_2];
        v_ex += field.At(ind_l, n);
        ind_l[dir_1] -= off_rel[dir_1];
        v_ex += field.At(ind_l, n);
        v_ex *= 0.125;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }
}

TEST_F(IntegrationTestCartesianInterpolation, InterpolateToScalarFaceTest) {
    Options opt_source_scalar{0, 0, 0};
    Options opt_source_staggerX{1, 0, 0};
    Options opt_source_staggerY{0, 1, 0};
    Options opt_source_staggerZ{0, 0, 1};
    Options opt_source_staggerXYZ{1, 1, 1};
    Options opt_target{0, 0, 0};

    using Representation = typename dare::Grid::Cartesian<Dim>::Representation;
    Representation grep_source_scalar = grid->GetRepresentation(opt_source_scalar);
    Representation grep_source_staggerX = grid->GetRepresentation(opt_source_staggerX);
    Representation grep_source_staggerY = grid->GetRepresentation(opt_source_staggerY);
    Representation grep_source_staggerZ = grid->GetRepresentation(opt_source_staggerZ);
    Representation grep_source_staggerXYZ = grid->GetRepresentation(opt_source_staggerXYZ);
    Representation grep_target = grid->GetRepresentation(opt_target);

    Field field_scalar("scalar_field", grep_source_scalar);
    Field field_x("x_field", grep_source_staggerX);
    Field field_y("y_field", grep_source_staggerY);
    Field field_z("z_field", grep_source_staggerZ);
    Field field_xyz("xyz_field", grep_source_staggerXYZ);

    Index ind(2, 3, 4);
    const double tol_eps{10.};

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * distribution(generator); };

    for (std::size_t n{0}; n < field_scalar.GetSize(); n++)
        field_scalar[n] = GetRandValue();
    for (std::size_t n{0}; n < field_x.GetSize(); n++)
        field_x[n] = GetRandValue();
    for (std::size_t n{0}; n < field_y.GetSize(); n++)
        field_y[n] = GetRandValue();
    for (std::size_t n{0}; n < field_z.GetSize(); n++)
        field_z[n] = GetRandValue();
    for (std::size_t n{0}; n < field_xyz.GetSize(); n++)
        field_xyz[n] = GetRandValue();

    // WEST (Scalar -> Scalar)
    dare::utils::Vector<N, SC> v_m =
        dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_scalar);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_scalar, n);
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = 0.5 * (field_scalar.At(ind, n) + field_scalar.At(ind_l, n));
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (Xstaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_x);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_x, n);
        Index ind_l{ind};
        SC v_ex = field_x.At(ind_l, n);
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (Ystaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_y);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_y, n);
        Index ind_l{ind};
        SC v_ex = field_y.At(ind_l, n);
        ind_l.i() -= 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.i() += 1;
        v_ex += field_y.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (Zstaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_z);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_z, n);
        Index ind_l{ind};
        SC v_ex = field_z.At(ind_l, n);
        ind_l.i() -= 1;
        v_ex += field_z.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_z.At(ind_l, n);
        ind_l.i() += 1;
        v_ex += field_z.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (XYZstaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_xyz);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_xyz, n);
        Index ind_l{ind};
        SC v_ex = field_xyz.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.j() -= 1;
        v_ex += field_xyz.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Scalar -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_scalar);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_scalar, n);
        Index ind_l{ind};
        ind_l.k() += 1;
        SC v_ex = 0.5 * (field_scalar.At(ind, n) + field_scalar.At(ind_l, n));
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Xstaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_x);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_x, n);
        Index ind_l{ind};
        SC v_ex = field_x.At(ind_l, n);
        ind_l.i() += 1;
        v_ex += field_x.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_x.At(ind_l, n);
        ind_l.i() -= 1;
        v_ex += field_x.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Ystaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_y);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_y, n);
        Index ind_l{ind};
        SC v_ex = field_y.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.j() -= 1;
        v_ex += field_y.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Zstaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_z);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_z, n);
        Index ind_l{ind};
        ind_l.k() += 1;
        SC v_ex = field_z.At(ind_l, n);
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (XYZstaggered -> Scalar)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_xyz);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_xyz, n);
        Index ind_l{ind};
        ind_l.k() += 1;
        SC v_ex = field_xyz.At(ind_l, n);
        ind_l.i() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.i() -= 1;
        v_ex += field_xyz.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }
}

TEST_F(IntegrationTestCartesianInterpolation, InterpolateToXStaggeredFaceTest) {
    Options opt_source_scalar{0, 0, 0};
    Options opt_source_staggerX{1, 0, 0};
    Options opt_source_staggerY{0, 1, 0};
    Options opt_source_staggerZ{0, 0, 1};
    Options opt_source_staggerXYZ{1, 1, 1};
    Options opt_target{1, 0, 0};

    using Representation = typename dare::Grid::Cartesian<Dim>::Representation;
    Representation grep_source_scalar = grid->GetRepresentation(opt_source_scalar);
    Representation grep_source_staggerX = grid->GetRepresentation(opt_source_staggerX);
    Representation grep_source_staggerY = grid->GetRepresentation(opt_source_staggerY);
    Representation grep_source_staggerZ = grid->GetRepresentation(opt_source_staggerZ);
    Representation grep_source_staggerXYZ = grid->GetRepresentation(opt_source_staggerXYZ);
    Representation grep_target = grid->GetRepresentation(opt_target);

    Field field_scalar("scalar_field", grep_source_scalar);
    Field field_x("x_field", grep_source_staggerX);
    Field field_y("y_field", grep_source_staggerY);
    Field field_z("z_field", grep_source_staggerZ);
    Field field_xyz("xyz_field", grep_source_staggerXYZ);

    Index ind(2, 3, 4);
    const double tol_eps{20.};

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * distribution(generator); };

    for (std::size_t n{0}; n < field_scalar.GetSize(); n++)
        field_scalar[n] = GetRandValue();
    for (std::size_t n{0}; n < field_x.GetSize(); n++)
        field_x[n] = GetRandValue();
    for (std::size_t n{0}; n < field_y.GetSize(); n++)
        field_y[n] = GetRandValue();
    for (std::size_t n{0}; n < field_z.GetSize(); n++)
        field_z[n] = GetRandValue();
    for (std::size_t n{0}; n < field_xyz.GetSize(); n++)
        field_xyz[n] = GetRandValue();

    // WEST (Scalar -> XStagger)
    dare::utils::Vector<N, SC> v_m =
        dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_scalar);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_scalar, n);
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = field_scalar.At(ind_l, n);
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (Xstaggered -> XStaggered)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_x);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_x, n);
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = 0.5 * (field_x.At(ind_l, n) + field_x.At(ind, n));
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (Ystaggered -> XStaggered)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_y);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_y, n);
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = field_y.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_y.At(ind_l, n);
        v_ex *= 0.5;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (Zstaggered -> XStaggered)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_z);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_z, n);
        Index ind_l{ind};
        ind_l.i() -= 1;
        SC v_ex = field_z.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_z.At(ind_l, n);
        v_ex *= 0.5;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // WEST (XYZstaggered -> XStaggered)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_xyz);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::WEST, field_xyz, n);
        Index ind_l{ind};
        SC v_ex = field_xyz.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.j() -= 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l = ind;
        ind_l.i() -= 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_xyz.At(ind_l, n);
        ind_l.j() -= 1;
        v_ex += field_xyz.At(ind_l, n);
        v_ex *= 0.125;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Scalar -> XStagger)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_scalar);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_scalar, n);
        Index ind_l{ind};
        SC v_ex = field_scalar.At(ind_l, n);
        ind_l.i() -= 1;
        v_ex += field_scalar.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_scalar.At(ind_l, n);
        ind_l.i() += 1;
        v_ex += field_scalar.At(ind_l, n);
        v_ex *= 0.25;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Xstaggered -> XStagger)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_x);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_x, n);
        Index ind_l{ind};
        SC v_ex = field_x.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_x.At(ind_l, n);
        v_ex *= 0.5;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Ystaggered -> XStaggered)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_y);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_y, n);
        Index ind_l{ind};
        SC v_ex = field_y.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.j() -= 1;
        v_ex += field_y.At(ind_l, n);
        ind_l = ind;
        ind_l.i() -= 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.j() += 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.k() += 1;
        v_ex += field_y.At(ind_l, n);
        ind_l.j() -= 1;
        v_ex += field_y.At(ind_l, n);
        v_ex *= 0.125;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }

    // TOP (Zstaggered -> XStagger)
    v_m = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_z);
    for (std::size_t n{0}; n < N; n++) {
        SC v = dare::math::InterpolateToFace(grep_target, ind, GridType::NeighborID::TOP, field_z, n);
        Index ind_l{ind};
        ind_l.k() += 1;
        SC v_ex = field_z.At(ind_l, n);
        ind_l.i() -= 1;
        v_ex += field_z.At(ind_l, n);
        v_ex *= 0.5;
        EXPECT_NEAR(v, v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
        EXPECT_NEAR(v_m[n], v_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_ex));
    }
}
