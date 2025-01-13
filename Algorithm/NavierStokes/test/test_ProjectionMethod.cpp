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

namespace dare::algorithm::test {
template <typename T>
using DensityInfo = dare::algorithm::PMDensityInfo<T>;

template <typename GridType, typename PDict, typename SDict>
using PM = dare::algorithm::ProjectionMethod<GridType, PDict, SDict>;

}  // namespace dare::algorithm::test

TEST(ProjectionMethodTest, PropertyInfo) {
    using GridType = dare::Grid::Cartesian<1>;
    using PDefault = dare::algorithm::PMPropertyInfoDefault<GridType>;
    using SDefault = dare::algorithm::PMSolverInfoDefault<GridType>;

    struct PDict_empty {
    };

    static_assert(std::is_same_v<
                  dare::algorithm::ProjectionMethod<GridType, PDict_empty, SDefault>::DensityInfo,
                  PDefault::density>);

    static_assert(std::is_same_v<
                  dare::algorithm::ProjectionMethod<GridType, PDict_empty, SDefault>::ViscosityInfo,
                  PDefault::viscosity
                  >);

    static_assert(std::is_same_v<
                  dare::algorithm::ProjectionMethod<GridType, PDict_empty, SDefault>::PorosityInfo,
                  PDefault::porosity>);

    static_assert(std::is_same_v<
                  dare::algorithm::ProjectionMethod<GridType, PDict_empty, SDefault>::ExplicitForceInfo,
                  PDefault::explicit_force>);

    static_assert(std::is_same_v<
                  dare::algorithm::ProjectionMethod<GridType, PDict_empty, SDefault>::ImplicitForceInfo,
                  PDefault::implicit_force>);

    static_assert(std::is_same_v<
                  dare::algorithm::ProjectionMethod<GridType, PDict_empty, SDefault>::CompressibilityInfo,
                  PDefault::compressible>);

}
