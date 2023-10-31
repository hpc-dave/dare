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

int main(int argc, char* argv[]) {
    const std::size_t Dim = 3;
    using LO = int32_t;
    using GO = int64_t;
    using SC = double;

    MPI_Init(&argc, &argv);

    MPI_Request request[2];
    // MPI_Request_free(&request);

    // std::cout << "trying MPI_wait" << std::endl;
    // MPI_Wait(&request, &status);

    // std::cout << "success" << std::endl;
    dare::mpi::ExecutionManager exman;
    std::vector<double> send_buffer(10), recv_buffer(10);

    for (int n = 0; n < 10; n++) {
        send_buffer[n] = 0. + n + exman.GetRank();
    }

    std::cout << "Rank " << exman.GetRank() << ": starting exchange" << std::endl;
    exman.Iexchange(send_buffer.data(), send_buffer.size(), recv_buffer.data(), recv_buffer.size(),
                    exman.GetRank(), 1000, &request[0], &request[1]);
    std::cout << "Rank " << exman.GetRank() << ": starting to wait" << std::endl;
    MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
    std::cout << "Rank " << exman.GetRank() << ": finished waiting" << std::endl;

    bool success{true};
    for (int n = 0; n < 10; n++)
        success &= send_buffer[n] == recv_buffer[n];

    if (success)
        std::cout << "Rank " << exman.GetRank() << ": success" << std::endl;
    else
        std::cout << "Rank " << exman.GetRank() << ": failure" << std::endl;
    GO i_size{10}, j_size{10};
    dare::utils::Vector<Dim, GO> res(30, 20, 10);
    dare::utils::Vector<Dim, LO> res_i = res;
    dare::utils::Vector<Dim, SC> size(1., 1., 1.);
    res_i = size;
    LO num_ghost = 2;
    dare::Grid::Cartesian<Dim> grid(&exman, res, size, num_ghost);
    auto rep = grid.GetRepresentation(dare::utils::Vector<Dim, LO>());
    rep.PrintDistribution("distribution.csv");
    MPI_Finalize();

    return 0;
}
