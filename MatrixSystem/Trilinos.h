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

#ifndef MATRIXSYSTEM_TRILINOS_H_
#define MATRIXSYSTEM_TRILINOS_H_
#include <vector>

#include "mpi.h"
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
// Tpetra  -- Vectors and Matrices
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Version.hpp>

#include "Data/DefaultTypes.h"
#include "Data/GridVector.h"
#include "MPI/ExecutionManager.h"
#include "Utilities/InitializationTracker.h"
#include "MatrixBlock.h"

namespace dare::Matrix {

/*!
 * @brief holds matrices, vectors and preconditioners for solving a LES
 * @tparam SC type of scalar value
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 */
template <typename SC>
class Trilinos : public dare::utils::InitializationTracker {
public:
    using ScalarType = SC;
    using LocalOrdinalType = dare::defaults::LocalOrdinalType;
    using GlobalOrdinalType = dare::defaults::GlobalOrdinalType;
    using GO = GlobalOrdinalType;
    using LO = LocalOrdinalType;
    using MatrixType = Tpetra::CrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using OpType = Tpetra::Operator<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using VecType = Tpetra::Vector<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using MultiVecType = Tpetra::MultiVector<ScalarType, LocalOrdinalType, GlobalOrdinalType>;
    using MapType = Tpetra::Map<LocalOrdinalType, GlobalOrdinalType>;
    using Communicator = Teuchos::Comm<int>;
    using LOViewType = typename MatrixType::nonconst_local_inds_host_view_type;
    using GOViewType = typename MatrixType::nonconst_global_inds_host_view_type;
    using SViewType =  typename MatrixType::nonconst_values_host_view_type;
    using constLOViewType = typename MatrixType::local_inds_host_view_type;
    using constGOViewType = typename MatrixType::global_inds_host_view_type;
    using constSViewType = typename MatrixType::values_host_view_type;

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
     * @tparam Lambda functor with instruction of type (auto*):void
     * @param grid grid representation
     * @param field instance of the field
     * @param functor instructions for building matrix
     * @param rebuild identifies if matrix graph changes
     * This function ONLY works on the host view and uses OpenMP directives
     * for acceleration, if you want to use the device defined by Trilinos,
     * write a custom build function.
     */
    template<typename Grid, std::size_t N, typename Lambda>
    void Build(const typename Grid::Representation& grid,
               const dare::Data::GridVector<Grid, SC, N>& field,
               Lambda functor,
               bool rebuild);

    /*!
     * @brief sets rhs vector according to functor
     * @tparam Grid grid type
     * @tparam Lambda functor of form (auto*):void
     * @tparam N number of components
     * @param grid grid representation
     * @param field reference to field
     * @param functor functor for matrix assembly
     */
    template <typename Grid, std::size_t N, typename Lambda>
    void SetB(const typename Grid::Representation& grid,
              const dare::Data::GridVector<Grid, SC, N>& field,
              Lambda functor);

    /*!
     * @brief Sets the initial guess to the specified value
     * @param value specific value
     */
    void SetInitialGuess(const SC value);

    /*!
     * @brief returns the matrix
     * @return Tpetra::CrsMatrix
     */
    Teuchos::RCP<MatrixType>& GetA();

    /*!
     * @brief returns RHS vector
     * @return Tpetra::Vector
     */
    Teuchos::RCP<VecType>& GetB();

    /*!
     * @brief returns initial guess and solution vector
     * @return Tpetra::Vector
     */
    Teuchos::RCP<VecType>& GetX();

    /*!
     * @brief returns preconditioner
     * @return Tpetra::Operator
     */
    Teuchos::RCP<OpType>& GetM();

    /*!
     * @brief returns map
     * @return const Tpetra::Map
     */
    Teuchos::RCP<const MapType>& GetMap();

    /*!
     * @brief returns the matrix
     * @return Tpetra::CrsMatrix
     */
    Teuchos::RCP<const MatrixType> GetA() const;

    /*!
     * @brief returns RHS vector
     * @return Tpetra::Vector
     */
    Teuchos::RCP<const VecType> GetB() const;

    /*!
     * @brief returns initial guess and solution vector
     * @return Tpetra::Vector
     */
    Teuchos::RCP<const VecType> GetX() const;

    /*!
     * @brief returns preconditioner
     * @return Tpetra::Operator
     */
    Teuchos::RCP<const OpType> GetM() const;

    /*!
     * @brief returns map
     * @return const Tpetra::Map
     */
    Teuchos::RCP<const MapType> GetMap() const;

private:
    /*!
     * @brief constructs new matrix from scratch and overwrites
     * @tparam Grid type of grid
     * @tparam Lambda functor with matrix block
     * @tparam N number of components
     * @param grid reference to grid
     * @param field reference to field with data
     * @param functor functor
     */
    template <typename Grid, std::size_t N, typename Lambda>
    void BuildNew(const typename Grid::Representation& grid,
                  const dare::Data::GridVector<Grid, SC, N>& field,
                  Lambda functor);

    /*!
     * @brief builds matrix, but reuses the existing one
     * @tparam Grid type of grid
     * @tparam Lambda type of functor
     * @tparam N number of components
     * @param grid reference to grid
     * @param field field with data
     * @param functor actual functor
     * You can only overwrite existing entries, the graph of
     * the matrix is constant!
     */
    template <typename Grid, std::size_t N, typename Lambda>
    void BuildReplace(const typename Grid::Representation& grid,
                      const dare::Data::GridVector<Grid, SC, N>& field,
                      Lambda functor);

    /*!
     * @brief updates matrix graph and entries
     * @tparam Grid type of grid
     * @tparam Lambda type of functor
     * @tparam N number of components
     * @param grid reference to grid
     * @param field reference to field
     * @param functor actual functor
     * currently call the BuildNew function
     */
    template <typename Grid, std::size_t N, typename Lambda>
    void BuildUpdate(const typename Grid::Representation& grid,
                     const dare::Data::GridVector<Grid, SC, N>& field,
                     Lambda functor);

    /*!
     * @brief allocates the map depending on the provided grid
     * @tparam Grid type of grid
     * @tparam N number of components
     * @param grid reference to grid
     * Right now, the local maps contain ONLY internal cells, excluding halo cells.
     * In principle, it might be preferable to include the halo cells to the local map,
     * sharing the entries between processes. Then looping through the internal cells
     * would only require local information and the build process for a matrix might be
     * sped up.
     */
    template <typename Grid, std::size_t N>
    void AllocateMap(const typename Grid::Representation& grid);

    dare::mpi::ExecutionManager* exec_man;  //!< pointer to execution manager
    Teuchos::RCP<const Communicator> comm;  //!< mpi communicator for Trilinos
    Teuchos::RCP<const MapType> map;        //!< map of row distribution
    Teuchos::RCP<MatrixType> A;             //!< Matrix
    Teuchos::RCP<VecType> x;                //!< Solution vector
    Teuchos::RCP<VecType> B;                //!< Right hand side vector
    Teuchos::RCP<OpType> M;                 //!< Preconditioner
    std::vector<GO> g_stencil;              //!< cells with global stencil
    std::vector<LO> l_stencil;              //!< cells with local stencil
};

}  // namespace dare::Matrix

#include "Trilinos.inl"

#endif  // MATRIXSYSTEM_TRILINOS_H_
