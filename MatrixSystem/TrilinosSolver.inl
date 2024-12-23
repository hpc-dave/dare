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
typename TrilinosSolver<SC>::ReturnType
TrilinosSolver<SC>::Solve(SolverPackage solver_pack,
                                  const std::string& type,
                                  Teuchos::RCP<MatrixType> A,
                                  Teuchos::RCP<VectorType> x,
                                  Teuchos::RCP<VectorType> B,
                                  Teuchos::RCP<ParameterList> param) {
    // adding empty preconditioner
    Teuchos::RCP<OperatorType> M;
    return Solve(solver_pack, type, M, A, x, B, param);
}

template <typename SC>
typename TrilinosSolver<SC>::ReturnType
TrilinosSolver<SC>::Solve(SolverPackage solver_pack,
                                  const std::string& type,
                                  Teuchos::RCP<OperatorType> M,
                                  Teuchos::RCP<MatrixType> A,
                                  Teuchos::RCP<VectorType> x,
                                  Teuchos::RCP<VectorType> B,
                                  Teuchos::RCP<ParameterList> param) {
    // need to convert the Vectors into MultiVectors to avoid template
    // deduction issues
    Teuchos::RCP<MultiVectorType> x_m(x);
    Teuchos::RCP<MultiVectorType> B_m(B);
    return Solve(solver_pack, type, M, A, x_m, B_m, param);
}

template <typename SC>
typename TrilinosSolver<SC>::ReturnType
TrilinosSolver<SC>::Solve(SolverPackage solver_pack,
                                  const std::string& type,
                                  Teuchos::RCP<MatrixType> A,
                                  Teuchos::RCP<MultiVectorType> x,
                                  Teuchos::RCP<MultiVectorType> B,
                                  Teuchos::RCP<ParameterList> param) {
    // adding empty preconditioner
    Teuchos::RCP<OperatorType> M;
    return Solve(solver_pack, type, M, A, x, B, param);
}

template <typename SC>
typename TrilinosSolver<SC>::ReturnType
TrilinosSolver<SC>::Solve(SolverPackage solver_pack,
                                  const std::string& type,
                                  Teuchos::RCP<OperatorType> M,
                                  Teuchos::RCP<MatrixType> A,
                                  Teuchos::RCP<MultiVectorType> x,
                                  Teuchos::RCP<MultiVectorType> B,
                                  Teuchos::RCP<ParameterList> param) {
    if (solver_pack == SolverPackage::Amesos2) {
        return SolveWithAmesos2(type, A, x, B, param);
    } else {
        Teuchos::RCP<SolverManager> solver = CreateSolver(solver_pack, type, param);
        Teuchos::RCP<ProblemType> problem = Teuchos::rcp(new ProblemType(A, x, B));

        if (!M.is_null())
            problem->setRightPrec(M);

        problem->setProblem();
        solver->setProblem(problem);
        auto ret = solver->solve();
        num_iter = solver->getNumIters();
        return ret;
    }
}

template <typename SC>
Teuchos::RCP<typename TrilinosSolver<SC>::OperatorType>
TrilinosSolver<SC>::BuildPreconditioner(PreCondPackage precond_package,
                                                const std::string& type,
                                                Teuchos::RCP<ParameterList> param,
                                                Teuchos::RCP<MatrixType> A) {
    Teuchos::RCP<OperatorType> M;
    if (precond_package == PreCondPackage::None)
        return M;

    switch (precond_package) {
    case PreCondPackage::Ifpack2:
        M = CreatePreconditionerIfPack2(type, param, A);
        break;
    case PreCondPackage::MueLu:
        M = CreatePreconditionerMueLu(type, param, A);
        break;
    case PreCondPackage::None:
        break;
    }

    if (M.is_null()) {
        int rank{-1};
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        bool am_i_root{rank == 0};
        if (am_i_root)
            ERROR << "Preconditioner of type " << type
                  << "is unknown or provided with wrong parameters!" << ERROR_CLOSE;
        MPI_Abort(MPI_COMM_WORLD, -5);
    }
    return M;
}

template <typename SC>
Teuchos::RCP<typename TrilinosSolver<SC>::SolverManager>
TrilinosSolver<SC>::CreateSolver(SolverPackage solver_pack,
                                         const std::string& type,
                                         Teuchos::RCP<ParameterList> param) {
    int rank{-1};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool am_i_root{rank == 0};

    Teuchos::RCP<SolverManager> sm;
    switch (solver_pack) {
        case SolverPackage::Amesos2:
            if (am_i_root)
                std::cerr << "Amesos solver should not end up here!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -5);
            break;
        case SolverPackage::Belos:
        {
            Belos::SolverFactory<SC, MultiVectorType, OperatorType> factory;
            sm = factory.create(type, param);
        }
            break;
        case SolverPackage::BumbleBee:
            if (type.compare("BICGSTAB2") == 0) {
                sm = Teuchos::rcp(new dare::Matrix::BiCGStab2<SC, MultiVectorType, OperatorType>(param));
            } else {
                if (am_i_root)
                    std::cerr << "Type: " << type << " is not recognized by BumbleBee!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, -5);
            }
            break;
        }

    return sm;
}

template <typename SC>
Teuchos::RCP<typename TrilinosSolver<SC>::OperatorType>
TrilinosSolver<SC>::CreatePreconditionerIfPack2(const std::string& type,
                                                        Teuchos::RCP<ParameterList> param,
                                                        Teuchos::RCP<const MatrixType> A) {
    using PrecType = Ifpack2::Preconditioner<SC>;
    Ifpack2::Factory factory;
    Teuchos::RCP<PrecType> prec = factory.create(type, A);
    prec->setParameters(*param);
    prec->initialize();
    prec->compute();
    return prec;
}

template <typename SC>
Teuchos::RCP<typename TrilinosSolver<SC>::OperatorType>
TrilinosSolver<SC>::CreatePreconditionerMueLu(const std::string& type,
                                                      Teuchos::RCP<ParameterList> param,
                                                      Teuchos::RCP<MatrixType> A) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    /*
     * Since MueLu is based on Xpetra, we need to wrap the matrix and nullspace
     * vector in an Xpetra-adapters.
     * Once the preconditioner is created, the Hierarchy is returned as a
     * Tpetra-operator
     */
    using MueLuTpetraAdapater = MueLu::TpetraOperator<SC>;
    using Hierarchy = MueLu::Hierarchy<SC>;
    using HierarchyManager = MueLu::HierarchyManager<SC>;
    using XpetraMatrix = Xpetra::Matrix<SC, LO, GO>;
    // typedef XpetraMultiVector = Xpetra::MultiVector<SC>;
    // typedef XpetraTpetraMVecAdapter = Xpetra::TpetraMultiVector<SC>;

    RCP<XpetraMatrix> A_xpetra = MueLu::TpetraCrs_To_XpetraMatrix<SC>(A);
    RCP<HierarchyManager> factory = rcp(new MueLu::ParameterListInterpreter<SC>(*param));
    RCP<Hierarchy> H = factory->CreateHierarchy();

    // // for now we just set it to one, later optimization may change this (then mulitvector is required)
    // LO num_dim = 1;
    // RCP<MultiVecType> nullspace = rcp(new MultiVecType(A_tpetra->getDomainMap(), 1));
    // RCP<XpetraMultiVector> nullspace_xpetra = rcp(new XpetraTpetraMVecAdapter(nullspace));

    // if (!filter.is_null()) {
    //     // here we apply a filter for the nullspace, should be revised eventually
    //     // RCP<NullspaceFilter> filter = properties.get<RCP<NullspaceFilter>>("AMG: NullspaceFilter");

    //     for (LO i = 0; i < num_dim; ++i) {
    //         Teuchos::ArrayRCP<SC> values = nullspace->getDataNonConst(i);
    //         LO num_blocks = values.size() / num_dim;
    //         for (LO j = 0; j < num_blocks; ++j) {
    //             if (filter->find(j * num_dim + i) == filter->end())
    //                 values[j * num_dim + i] = 1.0;
    //             else
    //                 values[j * num_dim + i] = 0.0;
    //         }
    //     }
    //     // remember, we need to provide the nullspace filter again at this point
    //     properties.set("AMG: NullspaceFilter", filter);
    // } else {
    //     // without filter we just set the nullspace to 1.0 everywhere
    //     // however, for later optimization, this algorithm was copied from the MueLu tutorial
    //     for (LO i = 0; i < num_dim; ++i) {
    //         Teuchos::ArrayRCP<SC> values = nullspace->getDataNonConst(i);
    //         LO num_blocks = values.size() / num_dim;
    //         for (LO j = 0; j < num_blocks; ++j) {
    //             values[j * num_dim + i] = 1.0;
    //         }
    //     }
    // }

    H->GetLevel(0)->Set("A", A_xpetra);
    // H->GetLevel(0)->Set("Nullspace", nullspace_xpetra);

    factory->SetupHierarchy(*H);

    H->IsPreconditioner(true);

    Teuchos::RCP<OperatorType> M = rcp(new MueLuTpetraAdapater(H));

    return M;
}

template <typename SC>
typename TrilinosSolver<SC>::ReturnType
TrilinosSolver<SC>::SolveWithAmesos2(const std::string& type,
                                             Teuchos::RCP<MatrixType> A,
                                             Teuchos::RCP<MultiVectorType> x,
                                             Teuchos::RCP<MultiVectorType> B,
                                             Teuchos::RCP<ParameterList> param) {
    Teuchos::RCP<Amesos2::Solver<MatrixType, MultiVectorType>> solver;
    solver = Amesos2::create<MatrixType, MultiVectorType>(type, A, x, B);
    solver->setParameters(param);
    solver->solve();
    return Belos::ReturnType::Converged;
}

template <typename SC>
int TrilinosSolver<SC>::GetNumIterations() const {
    return num_iter;
}

}  // end namespace dare::Matrix
