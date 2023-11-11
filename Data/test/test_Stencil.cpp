/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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
#include "../Stencil.h"

TEST(StencilTest, Initialize) {
    using SC = double;
    int size{7};
    dare::Data::Stencil<SC> stencil_empty, stencil_test("test", size);

    EXPECT_EQ(stencil_empty.GetSize(), 0);
    EXPECT_EQ(stencil_test.GetSize(), size);
}

TEST(StencilTest, GetSet) {
    using SC = double;
    int size{7};
    dare::Data::Stencil<SC> stencil_test("test", size);
    for (int n{0}; n < size; n++) {
        stencil_test.InsertValue(n, n);
        EXPECT_EQ(stencil_test.GetValue(n), static_cast<SC>(n));
    }

    for (int n{0}; n < size; n++) {
        stencil_test.GetValue(n) = static_cast<SC>(n + size);
        EXPECT_EQ(stencil_test.GetValue(n), static_cast<SC>(n + size));
    }
}

TEST(StencilTest, Copy) {
    using SC = double;
    int size{7};
    dare::Data::Stencil<SC> stencil_test("test", size);

    for (int n{0}; n < size; n++) {
        stencil_test.GetValue(n) = static_cast<SC>(n + size);
        ASSERT_EQ(stencil_test.GetValue(n), static_cast<SC>(n + size)) << "GetSet already broken";
    }

    dare::Data::Stencil<SC> stencil_construct(stencil_test), stencil_copy;
    stencil_copy = stencil_test;

    EXPECT_EQ(stencil_construct.GetSize(), size);
    EXPECT_EQ(stencil_copy.GetSize(), size);
    for (int n{0}; n < size; n++) {
        EXPECT_EQ(stencil_construct.GetValue(n), static_cast<SC>(n + size));
        EXPECT_EQ(stencil_copy.GetValue(n), static_cast<SC>(n + size));
    }
}
