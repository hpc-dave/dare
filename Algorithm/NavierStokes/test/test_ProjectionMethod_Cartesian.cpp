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

#include <type_traits>

#include "Algorithm/NavierStokes/ProjectionMethod.h"
#include "Grid/Cartesian.h"
#include "Algorithm/NavierStokes/ProjectionMethod_Cartesian.h"


namespace dare::test {

template <std::size_t Dim>
dare::Vector<Dim, defaults::GlobalOrdinalType> GetResolutionTestPMCartesian() {
    dare::Vector<Dim, defaults::GlobalOrdinalType> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 10 + n;
    return res;
}
template <std::size_t Dim>
dare::Vector<Dim, defaults::ScalarType> GetSizeTestPMCartesian() {
    dare::Vector<Dim, defaults::ScalarType> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

// a dummy for the boundary strategy
struct BStrat {
    template <typename T>
    void operator()(T t) {}
};

}  // namespace dare::test

/*!
 * @brief Fixture for testing the interpolation functions with the Cartesian Grid
 */
template<std::size_t D>
class ProjectionMethodCartesianTestFixture : public testing::Test {
public:
    static const std::size_t Dim{D};
    static const std::size_t N{Dim};
    using GridType = dare::Cartesian<Dim>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using Index = typename GridType::Index;
    using Field = dare::GridVector<GridType, SC, N>;
    using Options = typename GridType::Options;

    using PDefault = dare::PMPropertyInfoDefault<GridType>;
    using NDefault = dare::PMNumericalInfoDefault;
    using BStrat = dare::test::BStrat;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(&exec_man,
                                          dare::test::GetResolutionTestPMCartesian<Dim>(),
                                          dare::test::GetSizeTestPMCartesian<Dim>(),
                                          num_ghost);
    }

    std::unique_ptr<GridType> grid;   //!< the grid
    dare::ExecutionManager exec_man;  //!< the execution manager
};

using ProjectionMethodCartesian1DFixture = ProjectionMethodCartesianTestFixture<1>;
using ProjectionMethodCartesian2DFixture = ProjectionMethodCartesianTestFixture<2>;
using ProjectionMethodCartesian3DFixture = ProjectionMethodCartesianTestFixture<3>;

TEST_F(ProjectionMethodCartesian1DFixture, Initialization) {
    dare::ProjectionMethod<GridType, BStrat, PDefault, NDefault> pm;
    dare::test::BStrat bstrat;
    pm.Initialize(grid, bstrat);
}
