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

#include <Kokkos_Core.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <iostream>

#include "Grid/Cartesian.h"
#include "MPI/ExecutionManager.h"
#include "ScopeGuard/ScopeGuard.h"
#include "Utilities/Vector.h"
#include "mpi.h"
// Tpetra  -- Vectors and Matrices
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Version.hpp>
// Xpetra  -- Wrapper for dual use of Tpetra and Epetra (required by MueLu)
#include <Xpetra_CrsMatrix.hpp>
// Belos   -- Iterative solvers
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Parameters.hpp>

#include <MatrixSystem/Trilinos.h>
#include "Data/Stencil.h"

int main(int argc, char* argv[]) {
    const std::size_t Dim = 3;
    using LO = int32_t;
    using GO = int64_t;
    using SC = double;
    int rank_stop = 1;

    dare::ScopeGuard scope_guard(&argc, &argv);

    // Kokkos::initialize();
{
    std::size_t num_entries = 100000;
    Kokkos::View<dare::Data::Stencil<double>*> k("test", 0);

    Kokkos::resize(k, num_entries);

    Kokkos::parallel_for(
        num_entries, KOKKOS_LAMBDA(std::size_t node) {
            k[node].Resize(7);
        });
    dare::mpi::ExecutionManager exman;

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
}
    // Kokkos::finalize();
    return 0;
}
