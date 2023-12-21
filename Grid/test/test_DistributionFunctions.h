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

#ifndef GRID_TEST_TEST_DISTRIBUTIONFUNCTIONS_H_
#define GRID_TEST_TEST_DISTRIBUTIONFUNCTIONS_H_

#include <vector>
#include "../../Utilities/Vector.h"

namespace dare::Grid::test {
template <std::size_t Dim, typename T>
using Vector = dare::utils::Vector<Dim, T>;

namespace details {

template <std::size_t Dim, typename LO, typename GO>
void TestSumCells(const Vector<Dim, GO>& resolution_global, const std::vector<Vector<Dim, LO>>& vec_res_local) {
    GO num_cells_total{1};
    for (auto e : resolution_global)
        num_cells_total *= e;
    GO sum_cells{0};
    for (const auto& r_sub : vec_res_local) {
        LO n_sub{1};
        for (LO dim : r_sub)
            n_sub *= dim;
        sum_cells += n_sub;
    }
    EXPECT_EQ(sum_cells, num_cells_total)
        << "Total number of cells deviates in " << Dim << "D: "
        << num_cells_total << " != " << sum_cells << "!";
}

template <std::size_t Dim, typename LO, typename GO>
void TestForMissingCells(const Vector<Dim, GO>& resolution_global,
                         const std::vector<Vector<Dim, LO>>& vec_res_local,
                         const std::vector<Vector<Dim, GO>>& vec_offsets) {
    int num_proc = vec_res_local.size();
    GO num_cells_total{1};
    for (auto e : resolution_global)
        num_cells_total *= e;

    Vector<Dim, GO> glob_hsum;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        glob_hsum[dim] = 1;
        for (std::size_t n{dim + 1}; n < Dim; n++)
            glob_hsum[dim] *= resolution_global[n];
    }

    for (GO id{0}; id < static_cast<GO>(num_cells_total); id++) {
        GO found_id{0};
        // loop through all subdomains and find the ID
#pragma omp parallel for reduction(+ : found_id)
        for (int n_sub = 0; n_sub < num_proc; n_sub++) {
            Vector<Dim, GO> offset = vec_offsets[n_sub];
            Vector<Dim, GO> res_loc = vec_res_local[n_sub];
            GO num_loc_cells{1};
            for (std::size_t dim{0}; dim < Dim; dim++) {
                num_loc_cells *= res_loc[dim];
            }

            Vector<Dim, LO> loc_hsum;
            for (std::size_t dim{0}; dim < Dim; dim++) {
                loc_hsum[dim] = 1;
                for (std::size_t n{dim + 1}; n < Dim; n++)
                    loc_hsum[dim] *= res_loc[n];
            }

            for (LO loc_id{0}; loc_id < num_loc_cells; loc_id++) {
                LO n_loc{loc_id};
                Vector<Dim, LO> ind_loc;
                for (std::size_t dim{0}; dim < Dim; dim++) {
                    ind_loc[dim] = n_loc / loc_hsum[dim];
                    n_loc -= ind_loc[dim] * loc_hsum[dim];
                }
                Vector<Dim, GO> ind_glob = offset + ind_loc;
                GO id_glob{0};
                for (std::size_t dim{0}; dim < Dim; dim++)
                    id_glob += glob_hsum[dim] * ind_glob[dim];
                found_id += id_glob == id;
            }
        }
        EXPECT_FALSE(found_id == 0) << "Cell with ID = " << id
                                    << " was not found after distribution in " << Dim << "D!";
        EXPECT_FALSE(found_id > 1) << "Cell with ID = " << id
                                   << " was not found " << found_id << " after distribution in " << Dim << "D!";
        EXPECT_TRUE(found_id == 1) << "Error during distribution, see above for details";
    }
}
}  // namespace details
}  // namespace dare::Grid::test

class ConsistencyTest : public testing::TestWithParam<int> {
    // You can implement all the usual fixture class members here.
    // To access the test parameter, call GetParam() from class
    // TestWithParam<T>.
};

#endif  // GRID_TEST_TEST_DISTRIBUTIONFUNCTIONS_H_
