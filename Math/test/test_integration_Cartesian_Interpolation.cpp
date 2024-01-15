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

#include "Math/Interpolation_Cartesian.h"
#include "Grid/DefaultTypes.h"

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
    static const std::size_t N{3};
    static const std::size_t Dim{3};
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
