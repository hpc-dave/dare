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

#ifndef MATRIXSYSTEM_TRILINOSSOLVER_H_
#define MATRIXSYSTEM_TRILINOSSOLVER_H_

#include <string>
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

#include "BiCGStab2.h"
namespace dare::Matrix {

template <typename SC, typename LO, typename GO>
class TrilinosSolver {
public:
    enum class PreCondPackage {
        None,
        Ifpack2,
        MueLu
    };
    enum class SolverPackage {
        Amesos,
        Belos,
        BumbleBee
    };
    using ScalarType = SC;
    using LocalOrdinalType = LO;
    using GlobalOrdinalType = GO;
    using OperatorType = Tpetra::Operator<SC, LO, GO>;
    using MatrixType = Tpetra::CrsMatrix<SC, LO, GO>;
    using VectorType = Tpetra::Vector<SC, LO, GO>;
    using MultiVectorType = Tpetra::MultiVector<SC, LO, GO>;
    using ParameterList = Teuchos::ParameterList;
    using ReturnType = Belos::ReturnType;
    using SolverManager = Belos::SolverManager<ScalarType, MultiVectorType, OperatorType>;
    using ProblemType = Belos::LinearProblem<ScalarType, MultiVectorType, OperatorType>;

    /*!
     * @brief default constructor
     */
    TrilinosSolver() = default;

    /*!
     * @brief default destructor
     */
    virtual ~TrilinosSolver() = default;

    ReturnType Solve(SolverPackage solver_pack,
                     const std::string& type,
                     Teuchos::RCP<MatrixType> A,
                     Teuchos::RCP<VectorType> x,
                     Teuchos::RCP<VectorType> B,
                     Teuchos::RCP<ParameterList> param);

    ReturnType Solve(SolverPackage solver_pack,
                     const std::string& type,
                     Teuchos::RCP<OperatorType> M,
                     Teuchos::RCP<MatrixType> A,
                     Teuchos::RCP<VectorType> x,
                     Teuchos::RCP<VectorType> B,
                     Teuchos::RCP<ParameterList> param);

    ReturnType Solve(SolverPackage solver_pack,
                     const std::string& type,
                     Teuchos::RCP<MatrixType> A,
                     Teuchos::RCP<MultiVectorType> x,
                     Teuchos::RCP<MultiVectorType> B,
                     Teuchos::RCP<ParameterList> param);

    ReturnType Solve(SolverPackage solver_pack,
                     const std::string& type,
                     Teuchos::RCP<OperatorType> M,
                     Teuchos::RCP<MatrixType> A,
                     Teuchos::RCP<MultiVectorType> x,
                     Teuchos::RCP<MultiVectorType> B,
                     Teuchos::RCP<ParameterList> param);

    Teuchos::RCP<OperatorType> BuildPreconditioner(PreCondPackage precond_packag,
                                                   const std::string& type,
                                                   Teuchos::RCP<ParameterList> param,
                                                   Teuchos::RCP<MatrixType> A);

private:
    Teuchos::RCP<SolverManager> CreateSolver(SolverPackage solver_pack,
                                             const std::string& type,
                                             Teuchos::RCP<ParameterList> param);

    Teuchos::RCP<OperatorType> CreatePreconditionerIfPack2(const std::string& type,
                                                           Teuchos::RCP<ParameterList> param,
                                                           Teuchos::RCP<const MatrixType> A);

    Teuchos::RCP<OperatorType> CreatePreconditionerMueLu(const std::string& type,
                                                         Teuchos::RCP<ParameterList> param,
                                                         Teuchos::RCP<MatrixType> A);
};

}  // namespace dare::Matrix

#include "TrilinosSolver.inl"

#endif  // MATRIXSYSTEM_TRILINOSSOLVER_H_
