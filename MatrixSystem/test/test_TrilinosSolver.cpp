/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder
 *
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

#include <gtest/gtest.h>

#include "../../Data/GridVector.h"
#include "../../Grid/DefaultTypes.h"
#include "../Trilinos.h"
#include "../TrilinosSolver.h"
#include "test_TrilinosTestGrid.h"

namespace dare::Matrix::test {
struct _solverTestParam {
    using SC = double;
    using LO = dare::Grid::details::LocalOrdinalType;
    using GO = dare::Grid::details::GlobalOrdinalType;
    using PT = dare::Matrix::SolverPackage;
    using PTM = dare::Matrix::PreCondPackage;
    PT package;
    std::string type;
    PTM package_pre;
    std::string type_pre;

    _solverTestParam(PT pack, std::string typ, PTM pack_pre = PTM::None, std::string typ_pre = "none")
        : package(pack), type(typ), package_pre(pack_pre), type_pre(typ_pre) {}
};

Teuchos::RCP<Teuchos::ParameterList> GetPreconditionerParameters(_solverTestParam::PTM pack, std::string type) {
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::rcp(new Teuchos::ParameterList());
    switch (pack) {
    case _solverTestParam::PTM::Ifpack2:
        // paramters for ILU
        parameters->set("fact: drop tolerance", 1e-9);
        parameters->set("fact: level of fill", 1);
        parameters->set("schwarz: combine mode", "Add");
        break;
    case _solverTestParam::PTM::MueLu:
        parameters->set("problem: type", "MHD");  // works best in our cases
        parameters->set("verbosity", "none");
        break;
    otherwise:
        {}
    }
    return parameters;
}

}  // end namespace dare::Matrix::test

/*!
 * @brief templated fixture for testing multiple solvers
 */
class TrilinosSolverTest : public testing::TestWithParam<dare::Matrix::test::_solverTestParam> {
public:
    using GridType = dare::Matrix::test::TrilinosTestGrid;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using GridRepresentation = typename GridType::Representation;
    static const std::size_t N = dare::Matrix::test::N;
    using FieldType = dare::Data::GridVector<GridType, SC, N>;
    using GOViewType = typename dare::Matrix::Trilinos<SC>::GOViewType;
    using LOViewType = typename dare::Matrix::Trilinos<SC>::LOViewType;
    using SViewType = typename dare::Matrix::Trilinos<SC>::SViewType;
    GridType grid;
    FieldType field;
    dare::mpi::ExecutionManager exec_man;
    Teuchos::RCP<Teuchos::ParameterList> solver_param;
    bool multiSolves = false;
    void SetUp() {
        grid.Initialize(&exec_man);
        field = FieldType("test field", grid.GetRepresentation());
    }
};

TEST_P(TrilinosSolverTest, SolveLaplace) {
    dare::Matrix::test::_solverTestParam test_param = GetParam();

    if (test_param.package == dare::Matrix::test::_solverTestParam::PT::Amesos2) {
        solver_param = Teuchos::rcp(new Teuchos::ParameterList("Amesos2"));
        // using standard solver parameters
        // The KLU solver only works in serial
        if (!exec_man.IsSerial())
            return;
    } else {
        solver_param = Teuchos::rcp(new Teuchos::ParameterList());
        solver_param->set("Convergence Tolerance", 1e-12);
        solver_param->set("Maximum Iterations", 5000);
        solver_param->set("Num Blocks", 100);  // for GMRES
    }

    GridRepresentation g_rep{grid.GetRepresentation()};

    for (LO node = 0; node < grid.local_size; node++) {
        LO row = node * N;
        for (LO i{0}; i < N; i++) {
            field.At(node, i) = 0.5;
        }
    }

    // assembles a pseudo-2D problem, since the preconditioners of Ifpack and MueLu seem
    // to have issues with fully 1D problems
    auto functor = [&](auto mblock) {
        const std::size_t num_rows = grid.size_global * N;
        GO node_g = mblock->GetNode();
        LO node_l = g_rep.MapInternalToLocal(g_rep.MapGlobalToLocalInternal(node_g));
        bool is_left_edge = node_g == 0;
        bool is_right_edge = node_g == (g_rep.GetNumberGlobalCellsInternal()-1);
        for (std::size_t n{0}; n < N; n++) {
            std::size_t size_stencil{0};
            if (is_left_edge || is_right_edge) {
                size_stencil += 4;
            } else {
                size_stencil += 5;
            }
            if (n == 0)
                size_stencil--;
            if (n == (N - 1))
                size_stencil--;

            mblock->Resize(n, size_stencil);
            if (is_left_edge) {
                mblock->SetCoefficient(n, mblock->GetRow(n), 3.);
                mblock->SetCoefficient(n, mblock->GetRow(n) + N, -1.);
                mblock->SetRhs(n, 0.);
                mblock->SetInitialGuess(n, 0.);
            } else if (is_right_edge) {
                mblock->SetCoefficient(n, mblock->GetRow(n) - N, -1.);
                mblock->SetCoefficient(n, mblock->GetRow(n), 3.);
                mblock->SetRhs(n, 2.);
                mblock->SetInitialGuess(n, 1.);
            } else {
                mblock->SetCoefficient(n, mblock->GetRow(n) - N, -1.);
                mblock->SetCoefficient(n, mblock->GetRow(n), 2.);
                mblock->SetCoefficient(n, mblock->GetRow(n) + N, -1.);
                mblock->SetRhs(n, 0.);
            }
            if (n != 0) {
                mblock->SetCoefficient(n, mblock->GetRow(n) - 1, -1.);
                mblock->GetCoefficientByOrdinal(n, mblock->GetRow(n)) += 1;
            }
            if (n != (N-1)) {
                mblock->SetCoefficient(n, mblock->GetRow(n) + 1, -1.);
                mblock->GetCoefficientByOrdinal(n, mblock->GetRow(n)) += 1;
            }
        }
    };

    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    trilinos.Build(g_rep, field, functor, false);

    dare::Matrix::TrilinosSolver<SC> solver;

    auto precond_param = dare::Matrix::test::GetPreconditionerParameters(test_param.package_pre, test_param.type_pre);

    trilinos.GetM() = solver.BuildPreconditioner(test_param.package_pre,
                                                 test_param.type_pre,
                                                 precond_param,
                                                 trilinos.GetA());

    Belos::ReturnType ret = solver.Solve(test_param.package, test_param.type, trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(), solver_param);
    bool is_converged = ret == Belos::ReturnType::Converged;
    Tpetra::Vector<SC> r(trilinos.GetMap());
    Tpetra::Vector<SC> v(trilinos.GetMap());
    r.assign(*trilinos.GetB());

    trilinos.GetA()->apply(*trilinos.GetX(), v);
    r.update(-1., v, 1.);
    SC norm2 = r.norm2();

    ASSERT_TRUE(is_converged) << "Did not converge with solver: " << test_param.type
                              << ", preconditioner " << test_param.type_pre
                              << " and Euclidian norm: " << norm2;

    double d_phi = 1. / grid.size_global;
    for (LO node_l{0}; node_l < g_rep.GetNumberLocalCellsInternal(); node_l++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node_l);
        for (std::size_t n{0}; n < N; n++) {
            double phi = (node_g + 0.5) * d_phi;
            double val = trilinos.GetX()->getData()[node_l * N + n];
            EXPECT_TRUE(std::fabs(phi - val) < 1e-8) << "Expected: " << phi << " and found: " << val
                                                     << " deviation is: " << std::fabs(phi - val);
        }
    }
}

using dare::Matrix::test::_solverTestParam;

INSTANTIATE_TEST_SUITE_P(SolverConvergence,
                         TrilinosSolverTest,
                         testing::Values(
                             _solverTestParam(_solverTestParam::PT::BumbleBee, "BICGSTAB2"),
                             _solverTestParam(_solverTestParam::PT::BumbleBee, "BICGSTAB2",
                                              _solverTestParam::PTM::Ifpack2, "ILUT"),
                             _solverTestParam(_solverTestParam::PT::BumbleBee, "BICGSTAB2",
                                              _solverTestParam::PTM::MueLu, "AMG"),
                             _solverTestParam(_solverTestParam::PT::Belos, "CG"),
                             _solverTestParam(_solverTestParam::PT::Belos, "CG",
                                              _solverTestParam::PTM::Ifpack2, "ILUT"),
                             //  _solverTestParam(_solverTestParam::PT::Belos, "CG",
                             //                   _solverTestParam::PTM::MueLu, "AMG"),
                             _solverTestParam(_solverTestParam::PT::Belos, "BICGSTAB"),
                             _solverTestParam(_solverTestParam::PT::Belos, "BICGSTAB",
                                              _solverTestParam::PTM::Ifpack2, "ILUT"),
                             //  _solverTestParam(_solverTestParam::PT::Belos, "BICGSTAB",
                             //                   _solverTestParam::PTM::MueLu, "AMG"),
                             _solverTestParam(_solverTestParam::PT::Belos, "GMRES"),
                             _solverTestParam(_solverTestParam::PT::Belos, "GMRES",
                                              _solverTestParam::PTM::Ifpack2, "ILUT"),
                             _solverTestParam(_solverTestParam::PT::Belos, "GMRES",
                                              _solverTestParam::PTM::MueLu, "AMG"),
                             _solverTestParam(_solverTestParam::PT::Amesos2, "KLU")));
