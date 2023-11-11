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

#ifndef MATRIXSYSTEM_TRILINOS_H_
#define MATRIXSYSTEM_TRILINOS_H_

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
// Tpetra  -- Vectors and Matrices
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Version.hpp>
// Xpetra  -- Wrapper for dual use of Tpetra and Epetra (required by MueLu)
#include <Xpetra_CrsMatrix.hpp>

#include "../MPI/ExecutionManager.h"
#include "../Data/GridVector.h"

namespace dare::Matrix {

template <typename SC, typename LO, typename GO>
class Trilinos {
public:
    typedef ScalarType SC;
    typedef LocalOrdinalType LO;
    typedef GlobalOrdinalType GO;
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> MatrixType;
    typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal> OpType;
    typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> VecType;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal> MultiVecType;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal> MapType;
    typedef typename MatrixType::nonconst_local_inds_host_view_type LOViewType;
    typedef typename MatrixType::nonconst_global_inds_host_view_type GOViewType;
    typedef typename MatrixType::nonconst_values_host_view_type SViewType;

    explicit Trilinos(dare::mpi::ExecutionManager* exman);
    virtual ~Trilinos();

    /*!
     * @brief constructs matrix system according to functor
     * @tparam Grid type
     * @tparam Field type of field
     * @tparam Lambda functor with instruction
     * @param grid grid representation
     * @param field instance of the field
     * @param functor instructions for building matrix
     * @param rebuild identifies if matrix graph changes
     */
    template<typename Grid, std::size_t N, typename Lambda>
    void Build(const typename Grid::Representation& grid,
               const dare::Data::GridVector<Grid, SC, N>& field,
               Lambda functor,
               bool rebuild);

    const Teuchos::RCP<MatrixType>& GetA() const;
    const Teuchos::RCP<VectorType>& GetB() const;
    Teuchos::RCP<VectorType>& GetX();  // NOLINT
    const Teuchos::RCP<VectorType>& GetX() const;
    const Teuchos::RCP<OpType>& GetM() const;

private:
    dare::mpi::ExecutionManager* exec_man;  //!< pointer to execution manager
    Teuchos::RCP<MatrixType> A;             //!< Matrix
    Teuchos::RCP<VecType> x;                //!< Solution vector
    Teuchos::RCP<VecType> B;                //!< Right hand side vector
    Teuchos::RCP<OpType> M;                 //!< Preconditioner
};

}  // namespace dare::Matrix

#include "Trilinos.inl"

#endif  // MATRIXSYSTEM_TRILINOS_H_
