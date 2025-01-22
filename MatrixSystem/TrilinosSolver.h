/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder

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

#include <Amesos2.hpp>

#include "Utilities/Errors.h"
#include "Data/DefaultTypes.h"
#include "BiCGStab2.h"
namespace dare {

enum class PreCondPackage {
    None,
    Ifpack2,
    MueLu
};
enum class SolverPackage {
    Amesos2,
    Belos,
    BumbleBee
};

namespace detail {

/*!
 * @brief provides a default for the solver properties
 * @param type name of the solving algorithm
 */
inline Teuchos::RCP<Teuchos::ParameterList> GetDefaultSolverPropertiesTrilinos(std::string type = "BICGSTAB") {
    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList());
    p->set("Convergence Tolerance", 1e-16);
    p->set("Maximum Iterations", 1000);
    p->set("Num Blocks", 100);  // for GMRES
    return p;
}

/*!
 * @brief provdiesa a default for the preconditioner properties
 * @param type name fo the preconditioning algorithm
 */
inline Teuchos::RCP<Teuchos::ParameterList> GetDefaultPreconditionerPropertiesTrilinos(std::string type = "ILU") {
    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList());
    if (type.compare("ILU") == 0 || type.compare("ILUT") == 0 || type.compare("RILUK") == 0) {
        p->set("fact: drop tolerance", 1e-9);
        p->set("fact: level of fill", 1);
        p->set("schwarz: combine mode", "Add");
    } else if (type.compare("AMG") == 0) {
        p->set("problem: type", "MHD");  // works best in our cases
        p->set("verbosity", "none");
    } else if (type.compare("NONE") != 0) {
        ERROR << "Unknown type of preconditioner provided: " << type << ERROR_CLOSE;
    }
    return p;
}

/*!
 * @brief information struct for better transfer of information
 */
struct TrilinosNumericalProperties {
    using PropertyType = Teuchos::RCP<Teuchos::ParameterList>;
    TrilinosNumericalProperties()
        : solver_package(SolverPackage::Belos),
          solver_type("BICGSTAB"),
          precond_package(PreCondPackage::Ifpack2),
          precond_type("ILUT") {
        solver_properties = GetDefaultSolverPropertiesTrilinos(solver_type);
        precond_properties = GetDefaultPreconditionerPropertiesTrilinos(precond_type);
    }
    SolverPackage solver_package;       //!< the solver package of Trilinos
    std::string solver_type;            //!< solver type
    PropertyType solver_properties;     //!< detailed properties provided to the solver
    PreCondPackage precond_package;     //!< the preconditioner package
    std::string precond_type;           //!< preconditioner type
    PropertyType precond_properties;    //!< detailed properties of the preconditioner
};

}  // namespace detail

template <typename SC>
class TrilinosSolver {
public:
    using ScalarType = SC;
    using LocalOrdinalType = dare::defaults::LocalOrdinalType;
    using GlobalOrdinalType = dare::defaults::GlobalOrdinalType;
    using LO = LocalOrdinalType;
    using GO = GlobalOrdinalType;
    using OperatorType = Tpetra::Operator<SC, LO, GO>;
    using MatrixType = Tpetra::CrsMatrix<SC, LO, GO>;
    using VectorType = Tpetra::Vector<SC, LO, GO>;
    using MultiVectorType = Tpetra::MultiVector<SC, LO, GO>;
    using ParameterList = Teuchos::ParameterList;
    using PropertyType = Teuchos::RCP<ParameterList>;
    using ReturnType = Belos::ReturnType;
    using SolverManager = Belos::SolverManager<ScalarType, MultiVectorType, OperatorType>;
    using ProblemType = Belos::LinearProblem<ScalarType, MultiVectorType, OperatorType>;
    using NumericalPropertiesType = detail::TrilinosNumericalProperties;
    static const ReturnType Converged = ReturnType::Converged;

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

    ReturnType Solve(NumericalPropertiesType prop,
                     Teuchos::RCP<OperatorType> M,
                     Teuchos::RCP<MatrixType> A,
                     Teuchos::RCP<MultiVectorType> x,
                     Teuchos::RCP<MultiVectorType> B);

    ReturnType Solve(NumericalPropertiesType prop,
                     Teuchos::RCP<MatrixType> A,
                     Teuchos::RCP<MultiVectorType> x,
                     Teuchos::RCP<MultiVectorType> B);


    /*!
     * @brief builds preconditioner
     * @param precond_packag Trilinos specific package for the preconditioner
     * @param type the name of the preconditioner
     * @param param a parameter list for the specific preconditioner
     * @param A matrix on which the preconditioner is based
     */
    Teuchos::RCP<OperatorType> BuildPreconditioner(PreCondPackage precond_packag,
                                                   const std::string& type,
                                                   Teuchos::RCP<ParameterList> param,
                                                   Teuchos::RCP<MatrixType> A);

    /*!
     * @brief builds preconditioner with info structure
     * @param prop numerical property struct
     * @param A matrix for which the preconditioner will be build
     */
    Teuchos::RCP<OperatorType> BuildPreconditioner(NumericalPropertiesType prop,
                                                   Teuchos::RCP<MatrixType> A);

    /*!
     * @brief return number of iterations of last solving call
     */
    int GetNumIterations() const;

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

    ReturnType SolveWithAmesos2(const std::string& type,
                                Teuchos::RCP<MatrixType> A,
                                Teuchos::RCP<MultiVectorType> x,
                                Teuchos::RCP<MultiVectorType> B,
                                Teuchos::RCP<ParameterList> param);
    int num_iter{-1};
};

}  // namespace dare

#include "TrilinosSolver.inl"

#endif  // MATRIXSYSTEM_TRILINOSSOLVER_H_
