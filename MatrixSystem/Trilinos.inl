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

namespace dare::Matrix {

template <typename SC, typename LO, typename GO>
Trilinos<SC, LO, GO>::Trilinos(dare::mpi::ExecutionManager* exman)
    : exec_man(exman), g_stencil("global stencils", 0), l_stencil("local stencils", 0) {
    dare::utils::InitializationTracker::Initialize();
    comm = Teuchos::rcp(new Communicator(exec_man->GetCommunicator()));
}

template <typename SC, typename LO, typename GO>
Trilinos<SC, LO, GO>::~Trilinos() {}

template <typename SC, typename LO, typename GO>
template <typename GridRepresentation, std::size_t N, typename Lambda>
void Trilinos<SC, LO, GO>::Build(const typename Grid::Representation& grid,
                                 const dare::Data::GridVector<Grid, SC, N>& field,
                                 Lambda functor,
                                 bool rebuild) {
    exec_man->Terminate(__func__, "not implemented");

    if (A.is_null() || B.is_null() || x.is_null())
        BuildNew(grid, field, functor);
    else if (rebuild)
        BuildReplace(grid, field, functor);
    else
        BuildUpdate(grid, field, functor);

    const LO num_cells{grid.GetNumberLocalCellsInternal()};

    auto InsertMatrixEntries = KOKKOS_LAMBDA(const LO node){};

    Kokkos::parallel_for(num_cells, InsertMatrixEntries);
}

template <typename SC, typename LO, typename GO>
template <typename GridRepresentation, std::size_t N, typename Lambda>
void Trilinos<SC, LO, GO>::SetB(const typename Grid::Representation& grid,
                                const dare::Data::GridVector<Grid, SC, N>& field,
                                Lambda functor) {
    Kokkos::parallel_for(
        num_cells,
        KOKKOS_LAMBDA(const LO node) {
            if (node < l_stencil.size()) {
                LO local_internal = l_stencil[node];
                MatrixBlock<Grid, LO, SC, N> matrix_block(grid, local_internal);
                dare::utils::Vector<N, std::size_t> size_hint;
                for (std::size_t n{0}; n < N; n++)
                    size_hint[n] = A->getNumEntriesInLocalRow(matrix_block.GetRow(n));
                matrix_block.ProvideSizeHint(size_hint);
                // Set initial guess
                const LO local_full = grid.MapInternalToLocal(local_internal);
                for (std::size_t n{0}; n < N; n++)
                    matrix_block.SetInitialGuess(n, field.At(local_full * N + n));
                // call functor
                functor(&matrix_block);
                for (std::size_t n{0}; n < N; n++) {
                    B->replaceLocalValue(matrix_block.GetLocalRow(n),
                                         matrix_block.GetRhs(n));
                }
            } else {
                // initialize matrix block
                const LO local_internal = g_stencil[node - l_stencil.size()];
                const GO global_internal = grid.MapLocalToGlobalInternal(local_internal);
                MatrixBlock<Grid, GO, SC, N> matrix_block(grid, global_node);
                dare::utils::Vector<N, std::size_t> size_hint;
                for (std::size_t n{0}; n < N; n++)
                    size_hint[n] = A->getNumEntriesInGlobalRow(matrix_block.GetRow(n));
                matrix_block.ProvideSizeHint(size_hint);
                const LO local_full = grid.MapInternalToLocal(local_internal);
                for (std::size_t n{0}; n < N; n++)
                    matrix_block.SetInitialGuess(n, field.At(local_full * N + n));
                // call functor
                functor(&matrix_block);
                for (std::size_t n{0}; n < N; n++) {
                    B->replaceGlobalValue(matrix_block.GetLocalRow(n),
                                         matrix_block.GetRhs(n));
                }
            }
        });
}

template <typename SC, typename LO, typename GO>
const Teuchos::RCP<Trilinos<SC, LO, GO>::MatrixType>& Trilinos<SC, LO, GO>::GetA() const {
    return A;
}

template <typename SC, typename LO, typename GO>
const Teuchos::RCP<Trilinos<SC, LO, GO>::VectorType>& Trilinos<SC, LO, GO>::GetB() const {
    return B;
}

template <typename SC, typename LO, typename GO>
const Teuchos::RCP<Trilinos<SC, LO, GO>::VectorType>& Trilinos<SC, LO, GO>::GetX() const {
    return x;
}

template <typename SC, typename LO, typename GO>
Teuchos::RCP<Trilinos<SC, LO, GO>::VectorType>& Trilinos<SC, LO, GO>::GetX() {
    return x;
}

template <typename SC, typename LO, typename GO>
const Teuchos::RCP<Trilinos<SC, LO, GO>::OpType>& Trilinos<SC, LO, GO>::GetM() const {
    return M;
}

template <typename SC, typename LO, typename GO>
void Trilinos<SC, LO, GO>::SetInitialGuess(const SC value) {
    x->putScalar(value);
}

template <typename SC, typename LO, typename GO>
template <typename Grid, std::size_t N, typename Lambda>
void Trilinos<SC, LO, GO>::BuildNew(const typename Grid::Representation& grid,
                                    const dare::Data::GridVector<Grid, SC, N>& field,
                                    Lambda functor) {
    using MatrixBlockType = dare::Matrix::MatrixBlock<Grid, GO, SC, N>;
    if (map.is_null())
        exec_man->Terminate(__func__, "Cannot build new matrix system if map is not initialized");

    x = Teuchos::rcp(new VecType(map));
    B = Teuchos::rcp(new VecType(map));

    VecType* ptr_B = B.ptr();
    VecType* ptr_x = x.ptr();

    const LO num_cells{grid.GetNumberLocalCellsInternal()};
    const std::size_t num_local_rows{num_cells * N};

    Allocate memory for stencils
    Kokkos::View<MatrixBlockType*> list_matrix_blocks("listOfMatrixBlocks", num_cells);

    Kokkos::parallel_for(num_cells,
                         KOKKOS_LAMBDA(const LO node) {
                              //  initialize matrix block
                              MatrixBlock* matrix_block = &list_matrix_block[node];
                              matrix_block->Initialize(grid, node);
                              // set initial guess
                              for (std::size_t n{0}; n < N; n++)
                                  matrix_block.SetInitialGuess(n, field.At(node * N + n));
                              //  call functor
                              functor(matrix_block);
                              //  replace initial guess and rhs
                              for (std::size_t n{0}; n < N; n++) {
                                 B->replaceLocalValue((matrix_block.GetLocalRow(n),
                                                       matrix_block.GetRhs(n));
                                 x->replaceLocalValue((matrix_block.GetLocalRow(n),
                                                       matrix_block.GetInitialGuess(n));
                             }
                         });

    /*
     * Determine the numbers of entries to properly allocate the
     * arrays required to initialize the matrix
     */
    std::size_t num_entries{0};
    std::size_t num_local_stencils{0};
    Kokkos::View<int8_t*> is_stencil_local("local stencils", num_cells);
    Kokkos::View<std::size_t*> num_entries_per_row("number entries per row", num_cells * N);
    Kokkos::parallel_reduce(
        num_cells,
        KOKKOS_LAMBDA(const std::size_t node, std::size_t& num_entries_loc, std::size_t num_global_stencils) {
            for (std::size_t n{0}; n < N; n++) {
                std::size_t n_entry = list_matrix_block[n].GetNumEntries(n);
                num_entries_loc += n_entry;
                num_entries_per_row[node * N + n] = n_entry;
            }
            bool is_local = list_matrix_block[node].IsStencilLocal();
            num_local_stencils += is_local;
            is_stencil_local[node] = is_local;
        },
        num_entries, num_local_stencils);
    Kokkos::resize(l_stencil, num_local_stencils);
    Kokkos::resize(g_stencil, num_cells - num_local_stencils);

    {
        // I'm too stupid to do that in parallel, so I do it by copying everything to the host and doing it
        // there in serial. Here, the arrays with identifiers, which stencils are local and which are global are
        // populated
        Kokkos::View<int8_t*>::HostMirror is_stencil_local_host = Kokkos::create_mirror_view_and_copy(is_stencil_local);
        Kokkos::View<LO*>::HostMirror l_stencil_host = Kokkos::create_mirror_view_and_copy(l_stencil);
        Kokkos::View<LO*>::HostMirror g_stencil_host = Kokkos::create_mirror_view_and_copy(g_stencil);
        std::size_t count_l{0}, count_g{0};
        for (LO n{0}; n < num_cells; n++) {
            if (is_stencil_local_host[n]) {
                l_stencil_host[count_l] = n;
                ++count_l;
            } else {
                g_stencil_host[count_g] = n;
                ++count_g;
            }
        }
        if (count_l != num_local_stencils) {
            this->Terminate(__func__, "counted local stencils does not coincide with previously determined number!");
        }
        Kokkos::deep_copy(l_stencil, l_stencil_host);
        Kokkos::deep_copy(g_stencil, g_stencil_host);
    }

    /*
     * We need some offsets to work with, when substituting the values in the vectors
     */
    Kokkos::View<std::size_t*> row_offsets("row offsets", num_rows + 1);
    // for now serial, reconsider if really performance critical
    {
        Kokkos::View<std::size_t*>::HostMirror roff_h = Kokkos::create_mirror_view_and_copy(row_offsets);
        Kokkos::View<std::size_t*>::HostMirror npr_host = Kokkos::create_mirror_view_and_copy(num_entries_per_row);
        roff_h[0] = 0;
        for (std::size_t n{0}; n < num_cells*N; n++) {
            row_offsets[n + 1] = roff_h[n] + npr_host(num_row + n);
        }
        Kokkos::deep_copy(row_offsets, roff_h);  // copy back to device view
    }
    Kokkos::View<LO*> column_indices("column indices", num_entries_loc);
    Kokkos::View<SC*> matrix_values("matrix values", num_entries_loc);

    Kokkos::parallel_for(
        num_cells,
        KOKKOS_LAMBDA(LO node) {
            LO row = node * N;
            for (std::size_t n{0}; n < N; n++) {
                const auto& columns = list_matrix_block[n].GetColumnOrdinals(n);
                const auto& values = list_matrix_block[n].GetColumnValues(n);
                std::size_t offset = row_offsets[row];
                for (std::size_t num_entry{0}; num_entry < list_matrix_block[n].GetNumEntries(n); num_entry++) {
                    column_indices[offset] = columns[num_entry];
                    matrix_values[offset] = values[num_entry];
                    ++offset;
                }
                ++row;
            }
        });
    A = Teuchos::rcp(new MatrixType(map, map, row_offsets, column_indices, matrix_values));
    A->fillComplete();
}

template <typename SC, typename LO, typename GO>
template <typename Grid, std::size_t N, typename Lambda>
void Trilinos<SC, LO, GO>::BuildReplace(const typename Grid::Representation& grid,
                                        const dare::Data::GridVector<Grid, SC, N>& field,
                                        Lambda functor) {
    MatrixType* ptr_A = A.ptr();
    VectorType* ptr_x = x.ptr();
    VectorType* ptr_B = B.ptr();
    const LO num_cells{grid.GetNumberLocalCellsInternal()};

    if (A->isFillComplete())
        A->resumeFill();

    Kokkos::parallel_for(
        num_cells,
        KOKKOS_LAMBDA(const LO node) {
            if (node < l_stencil.size()) {
                // initialize matrix block
                LO local_internal = l_stencil[node];
                MatrixBlock<Grid, LO, SC, N> matrix_block(grid, local_internal);
                dare::utils::Vector<N, std::size_t> size_hint;
                for (std::size_t n{0}; n < N; n++)
                    size_hint[n] = A->getNumEntriesInLocalRow(matrix_block.GetRow(n));
                matrix_block.ProvideSizeHint(size_hint);
                // Set initial guess
                const LO local_full = grid.MapInternalToLocal(local_internal);
                for (std::size_t n{0}; n < N; n++)
                    matrix_block.SetInitialGuess(n, field.At(local_full * N + n));
                // call functor
                functor(&matrix_block);
                // replace matrix values
                for (std::size_t n{0}; n < N; n++) {
                    A->replaceLocalValues(matrix_block.GetRow(n),
                                          matrix_block.GetColumnOrdinals(n),
                                          matrix_block.GetColumnValues(n));
                B->replaceLocalValue((matrix_block.GetRow(n),
                                      matrix_block.GetRhs(n));
                x->replaceLocalValue((matrix_block.GetRow(n),
                                      matrix_block.GetInitialGuess(n));
                }
            } else {
                // initialize matrix block
                const LO local_internal = g_stencil[node - l_stencil.size()];
                const GO global_internal = grid.MapLocalToGlobalInternal(local_internal);
                MatrixBlock<Grid, GO, SC, N> matrix_block(grid, global_node);
                dare::utils::Vector<N, std::size_t> size_hint;
                for (std::size_t n{0}; n < N; n++)
                    size_hint[n] = A->getNumEntriesInGlobalRow(matrix_block.GetRow(n));
                matrix_block.ProvideSizeHint(size_hint);
                const LO local_full = grid.MapInternalToLocal(local_internal);
                for (std::size_t n{0}; n < N; n++)
                    matrix_block.SetInitialGuess(n, field.At(local_full * N + n));
                // call functor
                functor(&matrix_block);
                // replace matrix values
                for (std::size_t n{0}; n < N; n++) {
                    A->replaceGlobalValues(matrix_block.GetRow(n),
                                           matrix_block.GetColumnOrdinals(n),
                                           matrix_block.GetColumnValues(n));
                    B->replaceGlobalValue(matrix_block.GetRow(n),
                                          matrix_block.GetRhs(n));
                    x->replaceGlobalValue(matrix_block.GetRow(n),
                                          matrix_block.GetInitialGuess(n));
                }
            }
        });

    A->fillComplete();
}

template <typename SC, typename LO, typename GO>
template <typename GridRepresentation, std::size_t N, typename Lambda>
void Trilinos<SC, LO, GO>::BuildUpdate(const typename Grid::Representation& grid,
                                       const dare::Data::GridVector<Grid, SC, N>& field,
                                       Lambda functor) {
    BuildNew(grid, field, functor);
}

template <typename SC, typename LO, typename GO>
template <typename Grid, std::size_t N>
void Trilinos<SC, LO, GO>::AllocateMap(const typename Grid::Representation& grid) {
    using Teuchos::rcp;

    if (!dare::utils::InitializationTracker::IsInitialized()) {
        int rank{-1};
        MPI_Comm_Rank(MPI_COMM_WORLD, &rank);
        if (rank == 0)
            std::cerr << "Trilinos matrix system was not initialized before allocating the map, aborting!" << std::endl;
        MPI_Abort();
    }

    /*
     * Here we define the global distribution of the rows via a map
     * Additionally we define, from where we start counting, for us 0 is the standard.
     * Only madmen and noobs start counting from 1!!
     * Finally we can create the map!
     */

    const GO num_global_elements = grid.GetNumberGlobalCellsInternal() * N;
    const LO num_local_elements = grid.GetNumberLocalCellsInternal() * N;
    const GO index_base = 0;  // C++ starts counting from 0, so we stay consistent

    // Compute the locally stored rows by looping through the grid
    Kokkos::View<const GO*> elements_on_proc("local_glob_IDs", num_local_elements);
    Kokkos::parallel_for(
        grid.GetNumberLocalCellsInternal(),
        KOKKOS_LAMBDA(const LO node) {
            GO row = grid.MapLocalToGlobalInternal(node) * N;
            for (std::size_t n{0}; n < N; n++)
                elements_on_proc[node * N + n] = row + n;
        });
    // allocate the mach
    map = rcp(new MapType(num_global_elements, elements_on_proc, index_base, comm));
}

}  // namespace dare::Matrix
