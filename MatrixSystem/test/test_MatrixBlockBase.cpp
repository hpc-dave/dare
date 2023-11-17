/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder
 *
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
#include "../MatrixBlockBase.h"

/*!
 * @brief Fixture for testing MatrixBlockBase
 */
class MatrixBlockBaseTest : public testing::Test {
public:
    using O = int;                      //!< ordinal type
    using SC = double;                  //!< scalar type
    static const std::size_t N = 5;     //!< number of components
    using SizeHintVector = dare::utils::Vector<N, std::size_t>;  // type of vector for size hint
};

TEST_F(MatrixBlockBaseTest, Initialization) {
    O node = 13;
    SizeHintVector size_hint;
    for (std::size_t n{0}; n < N; n++)
        size_hint[n] = 7;
    dare::Matrix::MatrixBlockBase<O, SC, N> mblock_default_construct;
    dare::Matrix::MatrixBlockBase<O, SC, N> mblock_init_construct(node, size_hint);
    dare::Matrix::MatrixBlockBase<O, SC, N> mblock_copy_construct(mblock_init_construct);

    mblock_default_construct.Initialize(node, size_hint);

    // check if row is consistent
    O row_expect = node * N;
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mblock_default_construct.GetRow(n), row_expect + n);
        EXPECT_EQ(mblock_init_construct.GetRow(n), row_expect + n);
        EXPECT_EQ(mblock_copy_construct.GetRow(n), row_expect + n);
    }

    // check if allocated sizes are consistent
    for (std::size_t n{0}; n < N; n++) {
        EXPECT_EQ(mblock_default_construct.GetNumEntries(n), size_hint[n]);
        EXPECT_EQ(mblock_init_construct.GetNumEntries(n), size_hint[n]);
        EXPECT_EQ(mblock_copy_construct.GetNumEntries(n), size_hint[n]);
    }
}

TEST_F(MatrixBlockBaseTest, SizeHint) {
    O node = 13;
    SizeHintVector size_hint;
    for (std::size_t n{0}; n < N; n++)
        size_hint[n] = n;
    dare::Matrix::MatrixBlockBase<O, SC, N> mblock(node, size_hint);

    for (std::size_t n{0}; n < N; n++) {
        ASSERT_EQ(mblock.GetNumEntries(n), n);
    }
}

TEST_F(MatrixBlockBaseTest, Copy) {
    O node = 13;
    SizeHintVector size_hint;
    for (std::size_t n{0}; n < N; n++)
        size_hint[n] = n;
    dare::Matrix::MatrixBlockBase<O, SC, N> mblock_src(node, size_hint);
    dare::Matrix::MatrixBlockBase<O, SC, N> mblock_dst;

    O col_count{0};
    SC val_count{0.};
    SC rhs_count{0.};
    SC init_count{0.};
    O d_col{2};
    SC d_val{1.5}, d_rhs{0.5}, d_init{0.7};
    for (std::size_t n{0}; n < N; n++) {
        for (std::size_t col_id{0}; col_id < mblock_src.GetNumEntries(n); col_id++) {
            mblock_src.SetCoefficient(n, col_count, val_count);
            col_count += d_col;
            val_count += d_val;
        }
        mblock_src.SetRhs(n, rhs_count);
        mblock_src.SetInitialGuess(n, init_count);
        rhs_count += d_rhs;
        init_count += d_init;
    }

    mblock_dst = mblock_src;

    col_count = 0;
    val_count = 0.;
    rhs_count = 0.;
    init_count = 0.;

    for (std::size_t n{0}; n < N; n++) {
        for (std::size_t col_id{0}; col_id < mblock_src.GetNumEntries(n); col_id++) {
            EXPECT_EQ(mblock_dst.GetCoefficientByOrdinal(n, col_count), val_count);
            col_count += d_col;
            val_count += d_val;
        }
        EXPECT_EQ(mblock_dst.GetRhs(n), rhs_count);
        EXPECT_EQ(mblock_dst.GetInitialGuess(n), init_count);
        EXPECT_EQ(mblock_dst.GetRow(n), node * N + n);
        rhs_count += d_rhs;
        init_count += d_init;
    }
}

TEST_F(MatrixBlockBaseTest, GettersAndSetter) {
    O node = 17;
    SizeHintVector size_hint;
    for (std::size_t n{0}; n < N; n++)
        size_hint[n] = n;
    dare::Matrix::MatrixBlockBase<O, SC, N> mblock(node, size_hint);

    O col_count{0};
    SC val_count{0.};
    SC rhs_count{0.};
    SC init_count{0.};
    O d_col{2};
    SC d_val{1.5}, d_rhs{0.5}, d_init{0.7};
    std::vector<O> cmp_col[N];
    std::vector<SC> cmp_coef[N];
    std::vector<SC> cmp_rhs, cmp_init;
    for (std::size_t n{0}; n < N; n++) {
        for (std::size_t col_id{0}; col_id < mblock.GetNumEntries(n); col_id++) {
            mblock.SetCoefficient(n, col_count, val_count);
            cmp_col[n].push_back(col_count);
            cmp_coef[n].push_back(val_count);
            col_count += d_col;
            val_count += d_val;
        }
        mblock.SetRhs(n, rhs_count);
        mblock.SetInitialGuess(n, init_count);
        cmp_rhs.push_back(rhs_count);
        cmp_init.push_back(init_count);
        rhs_count += d_rhs;
        init_count += d_init;
    }

    // Regular setters
    for (std::size_t n{0}; n < N; n++) {
        for (std::size_t col_id{0}; col_id < mblock.GetNumEntries(n); col_id++) {
            EXPECT_EQ(mblock.GetOrdinalByPosition(n, col_id), cmp_col[n][col_id]);
            EXPECT_EQ(mblock.GetCoefficientByPosition(n, col_id), cmp_coef[n][col_id]);
            EXPECT_EQ(mblock.GetCoefficientByOrdinal(n, cmp_col[n][col_id]), cmp_coef[n][col_id]);
            cmp_col[n].push_back(col_count);
            cmp_coef[n].push_back(val_count);
        }
        EXPECT_EQ(mblock.GetRhs(n), cmp_rhs[n]);
        EXPECT_EQ(mblock.GetInitialGuess(n), cmp_init[n]);
    }

    // Array based setters
    for (std::size_t n{0}; n < N; n++) {
        const auto& ordinals = mblock.GetColumnOrdinals(n);
        const auto& coefs = mblock.GetColumnValues(n);
        for (std::size_t i{0}; i < mblock.GetNumEntries(n); i++) {
            EXPECT_EQ(ordinals[i], cmp_col[n][i]);
            EXPECT_EQ(coefs[i], cmp_coef[n][i]);
        }
    }
}
