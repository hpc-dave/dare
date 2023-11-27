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

#include "mpi.h"
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
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
// MueLu   -- Multigrid solvers & preconditioners
#include <MueLu.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_TpetraOperator.hpp>

#include "../Data/GridVector.h"
#include "../MPI/ExecutionManager.h"
#include "../Utilities/InitializationTracker.h"
#include "MatrixBlock.h"

namespace dare::Matrix {

template <typename SC, typename LO, typename GO>
class Trilinos : public dare::utils::InitializationTracker {
public:
    using ScalarType = SC;
    using LocalOrdinalType = LO;
    using GlobalOrdinalType = GO;
    using MatrixType = Tpetra::CrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using OpType = Tpetra::Operator<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using VecType = Tpetra::Vector<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using MultiVecType = Tpetra::MultiVector<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using MapType = Tpetra::Map<LocalOrdinalType, GlobalOrdinalType>;
    using Communicator = Teuchos::Comm<int>;
    using LOViewType = typename MatrixType::nonconst_local_inds_host_view_type;
    using GOViewType = typename MatrixType::nonconst_global_inds_host_view_type;
    using SViewType =  typename MatrixType::nonconst_values_host_view_type;

    /*!
     * @brief default constructor
     */
    Trilinos();

    /*!
     * @brief initializing constructor
     * @param exman pointer to execution manager
     */
    explicit Trilinos(dare::mpi::ExecutionManager* exman);

    /*!
     * @brief default destructor
     */
    virtual ~Trilinos();

    /*!
     * @brief Initializes the object
     * @param exman reference to execution manager
     */
    void Initialize(dare::mpi::ExecutionManager* exman);

    /*!
     * @brief constructs matrix system according to functor
     * @tparam Grid type
     * @tparam Field type of field
     * @tparam Lambda functor with instruction
     * @param grid grid representation
     * @param field instance of the field
     * @param functor instructions for building matrix
     * @param rebuild identifies if matrix graph changes
     * \note synchronize the device with the host view prior to this call!
     */
    template<typename Grid, std::size_t N, typename Lambda>
    void Build(const typename Grid::Representation& grid,
               const dare::Data::GridVector<Grid, SC, N>& field,
               Lambda functor,
               bool rebuild);

    /*!
     * @brief sets rhs vector according to functor
     * @tparam Grid grid type
     * @tparam Lambda functor of form (MatrixBlock*):void
     * @tparam N number of components
     * @param grid grid representation
     * @param field reference to field
     * @param functor functor for matrix assembly
     */
    template <typename Grid, std::size_t N, typename Lambda>
    void SetB(const typename Grid::Representation& grid,
              const dare::Data::GridVector<Grid, SC, N>& field,
              Lambda functor);

    void SetInitialGuess(const SC value);

    const Teuchos::RCP<MatrixType>& GetA() const;
    const Teuchos::RCP<VecType>& GetB() const;
    Teuchos::RCP<VecType>& GetX();  // NOLINT
    const Teuchos::RCP<VecType>& GetX() const;
    const Teuchos::RCP<OpType>& GetM() const;

private:
    template <typename Grid, std::size_t N, typename Lambda>
    void BuildNew(const typename Grid::Representation& grid,
                  const dare::Data::GridVector<Grid, SC, N>& field,
                  Lambda functor);

    template <typename Grid, std::size_t N, typename Lambda>
    void BuildReplace(const typename Grid::Representation& grid,
                      const dare::Data::GridVector<Grid, SC, N>& field,
                      Lambda functor);

    template <typename Grid, std::size_t N, typename Lambda>
    void BuildUpdate(const typename Grid::Representation& grid,
                     const dare::Data::GridVector<Grid, SC, N>& field,
                     Lambda functor);

    template <typename Grid, std::size_t N>
    void AllocateMap(const typename Grid::Representation& grid);

    dare::mpi::ExecutionManager* exec_man;  //!< pointer to execution manager
    Teuchos::RCP<const Communicator> comm;        //!< mpi communicator for Trilinos
    Teuchos::RCP<const MapType> map;        //!< map of row distribution
    Teuchos::RCP<MatrixType> A;             //!< Matrix
    Teuchos::RCP<VecType> x;                //!< Solution vector
    Teuchos::RCP<VecType> B;                //!< Right hand side vector
    Teuchos::RCP<OpType> M;                 //!< Preconditioner
    Kokkos::View<LO*> g_stencil;   //!< cells with global stencil
    Kokkos::View<LO*> l_stencil;   //!< cells with local stencil
};

}  // namespace dare::Matrix

#include "Trilinos.inl"

#endif  // MATRIXSYSTEM_TRILINOS_H_
