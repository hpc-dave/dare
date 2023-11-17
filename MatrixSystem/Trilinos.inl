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
Trilinos<SC, LO, GO>::Trilinos(dare::mpi::ExecutionManager* exman) : exec_man(exman) {
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
        KOKKOS_LAMBDA(const LO node){
            // initialize matrix block
            // MatrixBlock matrix_block(grid, field, node);
            // dare::utils::Vector<N, std::size_hint> size_hint;
            // for (std::size_t n{0}; n < N; n++)
            //     size_hint[n] = A->getNumEntriesInLocalRow(matrix_block.GetLocalRow(n));
            // matrix_block.ProvideSizeHint(size_hint);
            // call functor
            // functor(&matrix_block);
            // replace matrix values
            // for (std::size_t n{0}; n < N; n++) {
            //     B->replaceLocalValue((matrix_block.GetLocalRow(n),
            //                           matrix_block.GetRhs(n));
            // }
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
template <typename GridRepresentation, std::size_t N, typename Lambda>
void Trilinos<SC, LO, GO>::BuildNew(const typename Grid::Representation& grid,
                                    const dare::Data::GridVector<Grid, SC, N>& field,
                                    Lambda functor) {
    if (map.is_null())
        exec_man->Terminate(__func__, "Cannot build new matrix system if map is not initialized");

    x = Teuchos::rcp(new VecType(map));
    B = Teuchos::rcp(new VecType(map));

    VecType* ptr_B = B.ptr();
    VecType* ptr_x = x.ptr();

    const LO num_cells{grid.GetNumberLocalCellsInternal()};
    const std::size_t num_local_rows{num_cells * N};

    // Allocate memory for stencils
    // Kokkos::View<MatrixBlock*> list_matrix_blocks("listOfMatrixBlocks", num_cells);

    Kokkos::parallel_for(num_cells,
                         KOKKOS_LAMBDA(const LO node){
                             // initialize matrix block
                             //  MatrixBlock* matrix_block = &list_matrix_block[node];
                             //  matrix_block->Initialize(grid, field, node);
                             // call functor
                             //  functor(matrix_block);
                             // replace initial guess and rhs
                             // for (std::size_t n{0}; n < N; n++) {
                             //     B->replaceLocalValue((matrix_block.GetLocalRow(n),
                             //                           matrix_block.GetRhs(n));
                             //     x->replaceLocalValue((matrix_block.GetLocalRow(n),
                             //                           matrix_block.GetInitialGuess(n));
                             // }
                         });

    /*
     * Determine the numbers of entries to properly allocate the
     * arrays required to initialize the matrix
     */
    std::size_t num_entries{0};
    Kokkos::parallel_reduce(
        num_cells,
        KOKKOS_LAMBDA(const std::size_t node, std::size_t& num_entries_loc) {
            // for (std::size_t n{0}; n < N; n++){
            //     num_entries_loc += list_matrix_block[n].GetNumEntries(n);
            // }
        });
    /*
     * We need some offsets to work with, when substituting the values in the vectors
     */
    Kokkos::View<std::size_t*> row_offsets("row offsets", num_rows + 1);
    row_offsets[0] = 0;
    // for now serial, reconsider if really performance critical
    for (std::size_t node{0}; node < num_cells; node++) {
        // LO num_row = node * N;
        // for (std::size_t n{0}; n < N; n++)
        //     row_offsets[num_row + n + 1] = row_offset[num_row + n] + list_matrix_block.GetNumEntries(n);
    }
    Kokkos::View<LO*> column_indices("column indices", num_entries_loc);
    Kokkos::View<SC*> matrix_values("matrix values", num_entries_loc);

    Kokkos::parallel_for(
        num_cells,
        KOKKOS_LAMBDA(LO node) {
            // LO row = node * N;
            // for (std::size_t n{0}; n < N; n++) {
            //     const auto& columns = list_matrix_block[n].GetColumnOrdinals(n);
            //     const auto& values = list_matrix_block[n].GetColumnValues(n);
            //     std::size_t offset = row_offsets[row];
            //     for (std::size_t num_entry{0}; num_entry < list_matrix_block[n].GetNumEntries(n); num_entry++) {
            //         column_indices[offset] = columns[num_entry];
            //         matrix_values[offset] = values[num_entry];
            //         ++offset;
            //     }
            //     ++row;
            // }
        });
    A = Teuchos::rcp(new MatrixType(map, map, row_offsets, column_indices, matrix_values);
    A->fillComplete();
}

template <typename SC, typename LO, typename GO>
template <typename GridRepresentation, std::size_t N, typename Lambda>
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
            // initialize matrix block
            // MatrixBlock matrix_block(grid, field, node);
            // dare::utils::Vector<N, std::size_t> size_hint;
            // for (std::size_t n{0}; n < N; n++)
            //     size_hint[n] = A->getNumEntriesInLocalRow(matrix_block.GetLocalRow(n));
            // matrix_block.ProvideSizeHint(size_hint);
            // call functor
            // functor(&matrix_block);
            // replace matrix values
            // for (std::size_t n{0}; n < N; n++) {
            //     A->replaceLocalValues(matrix_block.GetLocalRow(n),
            //                           matrix_block.GetColumnOrdinals(n),
            //                           matrix_block.GetColumnValues(n));
            //     B->replaceLocalValue((matrix_block.GetLocalRow(n),
            //                           matrix_block.GetRhs(n));
            //     x->replaceLocalValue((matrix_block.GetLocalRow(n),
            //                           matrix_block.GetInitialGuess(n));
            // }
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
