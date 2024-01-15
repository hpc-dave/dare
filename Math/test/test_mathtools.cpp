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

#include "Grid/DefaultTypes.h"
#include "Math/Pow.h"
#include "Math/Divisors.h"

TEST(MathTools, PowTest) {
    static_assert(dare::math::Pow<0, 2>() == 0);
    static_assert(dare::math::Pow<1, 2>() == 1);
    static_assert(dare::math::Pow<2, 0>() == 1);
    static_assert(dare::math::Pow<2, 1>() == 2);
    static_assert(dare::math::Pow<2, 2>() == 4);
    static_assert(dare::math::Pow<2, 3>() == 8);
    static_assert(dare::math::Pow<2, 4>() == 16);
    static_assert(dare::math::Pow<3, 3>() == 27);
    static_assert(dare::math::Pow<-1, 0>() == 1);
    static_assert(dare::math::Pow<-1, 1>() == -1);
    static_assert(dare::math::Pow<-1, 2>() == 1);
    static_assert(dare::math::Pow<-1, 3>() == -1);

    const double v{1.5};
    const double ONE{1.};
    EXPECT_EQ(dare::math::Pow<0>(v), ONE);
    EXPECT_EQ(dare::math::Pow<1>(v), v);
    EXPECT_EQ(dare::math::Pow<2>(v), v * v);
    EXPECT_EQ(dare::math::Pow<3>(v), v * v * v);
    EXPECT_EQ(dare::math::Pow<-1>(v), ONE / v);
    EXPECT_EQ(dare::math::Pow<-2>(v), (ONE / v / v));
    EXPECT_EQ(dare::math::Pow(v, 0), ONE);
    EXPECT_EQ(dare::math::Pow(v, 1), v);
    EXPECT_EQ(dare::math::Pow(v, 2), v * v);
    EXPECT_EQ(dare::math::Pow(v, -1), ONE / v);
    EXPECT_EQ(dare::math::Pow(v, -2), ONE / v / v);
}

TEST(MathTools, IsRootOfTest) {
    static_assert(!dare::math::IsRootOf<0, 1>());
    static_assert(dare::math::IsRootOf<1, 1>());
    static_assert(dare::math::IsRootOf<2, 2>());
    static_assert(dare::math::IsRootOf<2, 1>());
    static_assert(!dare::math::IsRootOf<3, 2>());
    static_assert(dare::math::IsRootOf<4, 2>());
    static_assert(!dare::math::IsRootOf<5, 2>());
    static_assert(!dare::math::IsRootOf<6, 2>());
    static_assert(dare::math::IsRootOf<8, 2>());
}

TEST(MathTools, DivisionTest) {
    constexpr double ONE{1.}, HALF{0.5};
    static_assert(dare::math::Divisor<double, 1>() == ONE);
    static_assert(dare::math::Divisor<double, 2>() == HALF);
    static_assert(dare::math::Divisor<double, 4>() == HALF * HALF);
    static_assert(dare::math::Divisor<double, 8>() == HALF * HALF * HALF);
}

TEST(MathTools, DivideTest) {
    const double v{1.5}, HALF{0.5};
    EXPECT_EQ(dare::math::Divide<1>(v), v);
    EXPECT_EQ(dare::math::Divide<2>(v), v * HALF);
    EXPECT_EQ(dare::math::Divide<4>(v), v * HALF * HALF);
    EXPECT_EQ(dare::math::Divide<8>(v), v * HALF * HALF * HALF);
}
