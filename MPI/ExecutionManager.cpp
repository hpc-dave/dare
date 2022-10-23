/*
 * MIT License
 *
 * Copyright (c) 2022 David Rieder

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

#include "ExecutionManager.h"

namespace dare::mpi {

ExecutionManager::ExecutionManager(MPI_Comm _communicator, Verbosity _output_level)
    : communicator(_communicator), output_level(_output_level) {
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &num_proc);

    is_root = rank == 0;

    is_serial = num_proc == 1;
};

ExecutionManager::~ExecutionManager(){};

std::ostream& ExecutionManager::operator()(Verbosity level) {
    return Print(level);
}

std::ostream& ExecutionManager::Print(Verbosity level) {
    if (level > output_level || !is_root)
        return black_hole_osteam;
    else
        return std::cout;
}

std::ostream& ExecutionManager::PrintAll(Verbosity level) {
    if (level > output_level)
        return black_hole_osteam;
    else
        return std::cout;
}

int ExecutionManager::GetNumberThreadsLocal() {
    return omp_get_num_threads();
}

int ExecutionManager::GetNumberThreadsAverage() {
    int num_threads = omp_get_max_threads();
    num_threads = Allsum(num_threads);

    num_threads /= GetNumberProcesses();

    return num_threads;
}

void ExecutionManager::Terminate(std::string function, std::string message) const {
    std::cerr << "Terminating execution due to error in " << function
              << " on rank " << GetRank() << "! Following error description was received: "
              << message << std::endl;
    MPI_Abort(communicator, EXIT_FAILURE);
}

void ExecutionManager::Barrier() {
    if (IsSerial())
        return;

    MPI_Barrier(MPI_COMM_WORLD);
}
}  // namespace dare::mpi
