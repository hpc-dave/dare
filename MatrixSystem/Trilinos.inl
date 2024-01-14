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

namespace dare::Matrix {
template <typename SC>
Trilinos<SC>::Trilinos()
    : exec_man(nullptr), g_stencil(0), l_stencil(0) {
}

template <typename SC>
Trilinos<SC>::Trilinos(dare::mpi::ExecutionManager* exman)
    : exec_man(exman), g_stencil(0), l_stencil(0) {
    dare::utils::InitializationTracker::Initialize();
    comm = Teuchos::rcp(new Teuchos::MpiComm<int>(exec_man->GetCommunicator()));
}

template <typename SC>
Trilinos<SC>::~Trilinos() {}

template <typename SC>
void Trilinos<SC>::Initialize(dare::mpi::ExecutionManager* exman) {
    dare::utils::InitializationTracker::Initialize();
    exec_man = exman;
    comm = Teuchos::rcp(new Teuchos::MpiComm<int>(exec_man->GetCommunicator()));
}

template <typename SC>
template <typename Grid, std::size_t N, typename Lambda>
void Trilinos<SC>::Build(const typename Grid::Representation& grid,
                                 const dare::Data::GridVector<Grid, SC, N>& field,
                                 Lambda functor,
                                 bool rebuild) {
    if (A.is_null() || B.is_null() || x.is_null())
        BuildNew(grid, field, functor);
    else if (rebuild)
        BuildReplace(grid, field, functor);
    else
        BuildUpdate(grid, field, functor);
}

template <typename SC>
template <typename Grid, std::size_t N, typename Lambda>
void Trilinos<SC>::SetB(const typename Grid::Representation& grid,
                                const dare::Data::GridVector<Grid, SC, N>& field,
                                Lambda functor) {
    const std::size_t num_cells = grid.GetNumberLocalCellsInternal();
    Teuchos::Ptr<VecType> ptr_B = B.ptr();
    Teuchos::Ptr<MatrixType> ptr_A = A.ptr();
#pragma omp parallel for
    for (LO node = 0; node < num_cells; node++) {
        if (node < l_stencil.size()) {
            LO local_internal = l_stencil[node];
            MatrixBlock<Grid, LO, SC, N> matrix_block(grid, local_internal);
            dare::utils::Vector<N, std::size_t> size_hint;
            for (std::size_t n{0}; n < N; n++)
                size_hint[n] = ptr_A->getNumEntriesInLocalRow(matrix_block.GetRow(n));
            matrix_block.ProvideSizeHint(size_hint);
            // Set initial guess
            const LO local_full = grid.MapInternalToLocal(local_internal);
            for (std::size_t n{0}; n < N; n++)
                matrix_block.SetInitialGuess(n, field.At(local_full, n));
            // call functor
            functor(&matrix_block);
            for (std::size_t n{0}; n < N; n++) {
                ptr_B->replaceLocalValue(matrix_block.GetLocalRow(n),
                                         matrix_block.GetRhs(n));
            }
        } else {
            // initialize matrix block
            const LO local_internal = g_stencil[node - l_stencil.size()];
            const GO global_internal = grid.MapLocalToGlobalInternal(local_internal);
            MatrixBlock<Grid, GO, SC, N> matrix_block(grid, global_internal);
            dare::utils::Vector<N, std::size_t> size_hint;
            for (std::size_t n{0}; n < N; n++)
                size_hint[n] = ptr_A->getNumEntriesInGlobalRow(matrix_block.GetRow(n));
            matrix_block.ProvideSizeHint(size_hint);
            const LO local_full = grid.MapInternalToLocal(local_internal);
            for (std::size_t n{0}; n < N; n++)
                matrix_block.SetInitialGuess(n, field.At(local_full, n));
            // call functor
            functor(&matrix_block);
            for (std::size_t n{0}; n < N; n++) {
                ptr_B->replaceGlobalValue(matrix_block.GetLocalRow(n),
                                          matrix_block.GetRhs(n));
            }
        }
    }
}

template <typename SC>
Teuchos::RCP<typename Trilinos<SC>::MatrixType>& Trilinos<SC>::GetA() {
    return A;
}

template <typename SC>
Teuchos::RCP<typename Trilinos<SC>::VecType>& Trilinos<SC>::GetB() {
    return B;
}

template <typename SC>
Teuchos::RCP<typename Trilinos<SC>::VecType>& Trilinos<SC>::GetX() {
    return x;
}

template <typename SC>
Teuchos::RCP<typename Trilinos<SC>::OpType>& Trilinos<SC>::GetM() {
    return M;
}

template <typename SC>
Teuchos::RCP<const typename Trilinos<SC>::MapType>& Trilinos<SC>::GetMap() {
    return map;
}

template <typename SC>
Teuchos::RCP<const typename Trilinos<SC>::MatrixType> Trilinos<SC>::GetA() const {
    return A;
}

template <typename SC>
Teuchos::RCP<const typename Trilinos<SC>::VecType> Trilinos<SC>::GetB() const {
    return B;
}

template <typename SC>
Teuchos::RCP<const typename Trilinos<SC>::VecType> Trilinos<SC>::GetX() const {
    return x;
}

template <typename SC>
Teuchos::RCP<const typename Trilinos<SC>::OpType> Trilinos<SC>::GetM() const {
    return M;
}

template <typename SC>
Teuchos::RCP<const typename Trilinos<SC>::MapType> Trilinos<SC>::GetMap() const {
    return map;
}

template <typename SC>
void Trilinos<SC>::SetInitialGuess(const SC value) {
    x->putScalar(value);
}

template <typename SC>
template <typename Grid, std::size_t N, typename Lambda>
void Trilinos<SC>::BuildNew(const typename Grid::Representation& grid,
                                    const dare::Data::GridVector<Grid, SC, N>& field,
                                    Lambda functor) {
    using MatrixBlockType = dare::Matrix::MatrixBlock<Grid, GO, SC, N>;
    using DViewLO = Kokkos::DualView<LO*>;
    using DViewGO = Kokkos::DualView<GO*>;
    using DViewSC = Kokkos::DualView<SC*>;
    using DViewSizet = Kokkos::DualView<std::size_t*>;
    using HostSpace = typename DViewGO::host_mirror_space;
    using ExecutionSpace = typename DViewGO::execution_space;

    AllocateMap<Grid, N>(grid);

    if (map.is_null())
        exec_man->Terminate(__func__, "Cannot build new matrix system if map is not initialized");

    x = Teuchos::rcp(new VecType(map));
    B = Teuchos::rcp(new VecType(map));

    Teuchos::Ptr<VecType> ptr_B = B.ptr();
    Teuchos::Ptr<VecType> ptr_x = x.ptr();

    const LO num_cells{grid.GetNumberLocalCellsInternal()};
    const std::size_t num_local_rows{num_cells * N};

    // Allocate memory for stencils
    std::vector<MatrixBlockType> list_matrix_blocks(num_cells);

#pragma omp parallel for
    for (LO node = 0; node < num_cells; node++) {
        //  initialize matrix block
        MatrixBlockType* matrix_block = &list_matrix_blocks[node];
        const GO node_global = grid.MapLocalToGlobalInternal(node);
        matrix_block->Initialize(&grid, node_global);
        // set initial guess
        LO node_full = grid.MapInternalToLocal(node);
        for (std::size_t n{0}; n < N; n++)
            matrix_block->SetInitialGuess(n, field.At(node_full, n));
        //  call functor
        functor(matrix_block);
        //  replace initial guess and rhs
        for (std::size_t n{0}; n < N; n++) {
            ptr_B->replaceGlobalValue(matrix_block->GetRow(n),
                                      matrix_block->GetRhs(n));
            ptr_x->replaceGlobalValue(matrix_block->GetRow(n),
                                      matrix_block->GetInitialGuess(n));
        }
    }

    /*
     * Determine the numbers of entries to properly allocate the
     * arrays required to initialize the matrix
     */
    std::size_t num_entries{0};
    std::size_t num_local_stencils{0};
    Teuchos::Array<int8_t> is_stencil_local(num_cells);
    Teuchos::Array<std::size_t> num_entries_per_row(num_cells * N);
#pragma omp parallel for reduction(+:num_entries, num_local_stencils)
    for (LO node = 0; node < num_cells; node++) {
        for (std::size_t n{0}; n < N; n++) {
            std::size_t n_entry = list_matrix_blocks[node].GetNumEntries(n);
            num_entries += n_entry;
            num_entries_per_row[node * N + n] = n_entry;
        }
        bool is_local = list_matrix_blocks[node].IsStencilLocal();
        num_local_stencils += is_local;
        is_stencil_local[node] = is_local;
    }

    l_stencil.resize(num_local_stencils);
    g_stencil.resize(num_cells - num_local_stencils);

    // TODO(@Dave): test if overhead
    std::size_t count_l{0}, count_g{0};
    for (LO n{0}; n < num_cells; n++) {
        if (is_stencil_local[n]) {
            l_stencil[count_l] = n;
            ++count_l;
        } else {
            g_stencil[count_g] = n;
            ++count_g;
        }
    }
    if (count_l != num_local_stencils) {
        exec_man->Terminate(__func__,
                            "counted local stencils does not coincide with previously determined number!");
    }

    A = Teuchos::rcp(new MatrixType(map, num_entries_per_row()));
    Teuchos::Ptr<MatrixType> ptr_A = A.ptr();
// #pragma omp parallel for
    for (LO node = 0; node < num_cells; node++) {
        GO row = grid.MapLocalToGlobalInternal(node) * N;
        for (std::size_t n{0}; n < N; n++) {
            const auto& columns = list_matrix_blocks[node].GetColumnOrdinals(n);
            const auto& values = list_matrix_blocks[node].GetColumnValues(n);
            Teuchos::ArrayView<GO> cview(columns.data(), columns.size());
            Teuchos::ArrayView<SC> vview(values.data(), values.size());
// #pragma omp critical
            ptr_A->insertGlobalValues(row, cview, vview);
            ++row;
        }
    }

    A->fillComplete();
}

template <typename SC>
template <typename Grid, std::size_t N, typename Lambda>
void Trilinos<SC>::BuildReplace(const typename Grid::Representation& grid,
                                        const dare::Data::GridVector<Grid, SC, N>& field,
                                        Lambda functor) {
    Teuchos::Ptr<MatrixType> ptr_A = A.ptr();
    Teuchos::Ptr<VecType> ptr_x = x.ptr();
    Teuchos::Ptr<VecType> ptr_B = B.ptr();
    const LO num_cells{grid.GetNumberLocalCellsInternal()};

    if (A->isFillComplete())
        A->resumeFill();

#pragma omp parallel for
    for (LO node = 0; node < num_cells; node++) {
        if (node < l_stencil.size()) {
            // initialize matrix block
            const LO local_internal = l_stencil[node];
            const LO local_full = grid.MapInternalToLocal(local_internal);
            MatrixBlock<Grid, LO, SC, N> matrix_block(&grid, local_full);
            dare::utils::Vector<N, std::size_t> size_hint;
            for (std::size_t n{0}; n < N; n++)
                size_hint[n] = A->getNumEntriesInLocalRow(matrix_block.GetRow(n));
            matrix_block.ProvideSizeHint(size_hint);
            // Set initial guess
            for (std::size_t n{0}; n < N; n++)
                matrix_block.SetInitialGuess(n, field.At(local_full, n));
            // call functor
            functor(&matrix_block);
            // replace matrix values
            for (std::size_t n{0}; n < N; n++) {
                LO col_ptr[2];
                SC data_ptr[2];
                LO row = 1;
                LO n_e = 1;
                ptr_A->replaceLocalValues(matrix_block.GetRow(n),
                                          matrix_block.GetColumnOrdinals(n),
                                          matrix_block.GetColumnValues(n));
                ptr_B->replaceLocalValue(matrix_block.GetRow(n),
                                         matrix_block.GetRhs(n));
                ptr_x->replaceLocalValue(matrix_block.GetRow(n),
                                         matrix_block.GetInitialGuess(n));
            }
        } else {
            // initialize matrix block
            const LO local_internal = g_stencil[node - l_stencil.size()];
            const LO local_full = grid.MapInternalToLocal(local_internal);
            const GO global_internal = grid.MapLocalToGlobalInternal(local_internal);
            const GO global_full = grid.MapInternalToLocal(global_internal);
            MatrixBlock<Grid, GO, SC, N> matrix_block(&grid, global_internal);
            dare::utils::Vector<N, std::size_t> size_hint;
            for (std::size_t n{0}; n < N; n++)
                size_hint[n] = ptr_A->getNumEntriesInGlobalRow(matrix_block.GetRow(n));
            matrix_block.ProvideSizeHint(size_hint);

            for (std::size_t n{0}; n < N; n++)
                matrix_block.SetInitialGuess(n, field.At(local_full, n));
            // call functor
            functor(&matrix_block);
            // replace matrix values
            for (std::size_t n{0}; n < N; n++) {
                ptr_A->replaceGlobalValues(matrix_block.GetRow(n),
                                           matrix_block.GetColumnOrdinals(n),
                                           matrix_block.GetColumnValues(n));
                ptr_B->replaceGlobalValue(matrix_block.GetRow(n),
                                          matrix_block.GetRhs(n));
                ptr_x->replaceGlobalValue(matrix_block.GetRow(n),
                                          matrix_block.GetInitialGuess(n));
            }
        }
    }

    A->fillComplete();
}

template <typename SC>
template <typename Grid, std::size_t N, typename Lambda>
void Trilinos<SC>::BuildUpdate(const typename Grid::Representation& grid,
                                       const dare::Data::GridVector<Grid, SC, N>& field,
                                       Lambda functor) {
    BuildNew(grid, field, functor);
}

template <typename SC>
template <typename Grid, std::size_t N>
void Trilinos<SC>::AllocateMap(const typename Grid::Representation& grid) {
    using Teuchos::rcp;

    if (!dare::utils::InitializationTracker::IsInitialized()) {
        int rank{-1};
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0)
            std::cerr << "Trilinos matrix system was not initialized before allocating the map, aborting!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
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
    Kokkos::DualView<GO*> elements_on_proc("local_glob_IDs", num_local_elements);
#pragma omp parallel for
    for (LO node = 0; node < grid.GetNumberLocalCellsInternal(); node++) {
        GO row = grid.MapLocalToGlobalInternal(node) * N;
        for (std::size_t n{0}; n < N; n++) {
            elements_on_proc.h_view[node * N + n] = row + n;
        }
    }
    elements_on_proc.template modify<typename Kokkos::DualView<GO*>::host_mirror_space>();
    elements_on_proc.template sync<typename Kokkos::DualView<GO*>::execution_space>();
    // allocate the map
    map = rcp(new MapType(num_global_elements, elements_on_proc.d_view, index_base, comm));
}

}  // namespace dare::Matrix
