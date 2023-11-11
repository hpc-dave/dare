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

    const LO num_cells{grid.GetNumberLocalCellsInternal()};

    auto InsertMatrixEntries = KOKKOS_LAMBDA(const LO node){
    };

    Kokkos::parallel_for(num_cells, InsertMatrixEntries);
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

}  // namespace dare::Matrix
