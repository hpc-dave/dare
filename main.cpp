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

#include <iostream>
#include "Utilities/Vector.h"
#include "MPI/ExecutionManager.h"
#include "Grid/Cartesian.h"
#include "ScopeGuard/ScopeGuard.h"

int main(int argc, char* argv[]) {
    const std::size_t Dim = 3;
    using LO = int32_t;
    using GO = int64_t;
    using SC = double;
    int rank_stop = 1;

    dare::ScopeGuard scope_guard(argc, argv);

    dare::mpi::ExecutionManager exman;

    // int dims[] = {0};
    // MPI_Dims_create(exman.GetNumberProcesses(), 1, dims);

    // if (exman.AmIRoot()) {
    //     for (int n = 0; n < 1; n++) {
    //         std::cout << dims[n] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    GO i_size{10}, j_size{10};
    dare::utils::Vector<Dim, GO> res(30, 20, 10);
    dare::utils::Vector<Dim, LO> res_i = res;
    dare::utils::Vector<Dim, SC> size(1., 1., 1.);
    res_i = size;
    LO num_ghost = 2;
    dare::Grid::Cartesian<Dim> grid(&exman, res, size, num_ghost, dare::utils::Vector<Dim, LO>(1, 0, 0));
    int stop = 1;
    if (exman.GetRank() == rank_stop)
        while (stop == 1) {
        }
    auto rep = grid.GetRepresentation(dare::utils::Vector<Dim, LO>());
    rep.PrintDistribution("distribution.csv");

    return 0;
}
