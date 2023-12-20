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

#include <gtest/gtest.h>

#include "../Grid/Cartesian.h"
#include "../Grid/DefaultTypes.h"
#include "../MatrixSystem/Trilinos.h"
#include "../MatrixSystem/TrilinosSolver.h"
#include "../Utilities/Vector.h"

namespace dare::test {

template <std::size_t Dim, typename GO>
dare::utils::Vector<Dim, GO> GetResolutionIntegrationTestCartesianTrilinos() {
    dare::utils::Vector<Dim, GO> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 5 + n;
    return res;
}
template <std::size_t Dim, typename SC>
dare::utils::Vector<Dim, SC> GetSizeIntegrationTestCartesianTrilinos() {
    dare::utils::Vector<Dim, SC> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

}  // namespace dare::test

template <std::size_t Dimension>
class IntegrationCartesianTrilinos
    : public testing::TestWithParam<typename dare::Grid::Cartesian<Dimension>::Options> {
public:
    static const std::size_t N{3};
    using GridType = dare::Grid::Cartesian<Dimension>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using SC = double;
    using GridVector = dare::Data::GridVector<GridType, SC, N>;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(
            &exec_man,
            dare::test::GetResolutionIntegrationTestCartesianTrilinos<Dimension, GO>(),
            dare::test::GetSizeIntegrationTestCartesianTrilinos<Dimension, SC>(),
            num_ghost);
    }

    std::unique_ptr<GridType> grid;
    dare::mpi::ExecutionManager exec_man;
};

using IntegrationCartesianTrilinos1D = IntegrationCartesianTrilinos<1>;
using IntegrationCartesianTrilinos2D = IntegrationCartesianTrilinos<2>;
using IntegrationCartesianTrilinos3D = IntegrationCartesianTrilinos<3>;

TEST_F(IntegrationCartesianTrilinos1D, SolveScalar) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC, LO, GO> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC, LO, GO> solver;
    const double value_west{0.};
    const double value_east{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";
    const dare::Matrix::PreCondPackage precond_pack = dare::Matrix::PreCondPackage::MueLu;
    const std::string precond_type = "AMG";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    Teuchos::RCP<Teuchos::ParameterList> param_precond = Teuchos::rcp(new Teuchos::ParameterList());
    param_precond->set("problem: type", "MHD");  // works best in our cases
    param_precond->set("verbosity", "none");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        GO node_g = mblock->GetNode();
        bool is_left_edge{node_g == 0};
        bool is_right_edge{node_g == (g_rep.GetGlobalResolutionInternal()[0] - 1)};
        std::size_t stencil_size{3};
        if (is_left_edge)
            stencil_size--;
        else if (is_right_edge)
            stencil_size--;
        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            mblock->GetRhs(n) = 0.;
            mblock->template Get<CN::CENTER>(n) = 2.;
            if (!is_left_edge) {
                mblock->template Get<CN::WEST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_west;
                mblock->SetInitialGuess(n, value_west);
            }

            if (!is_right_edge) {
                mblock->template Get<CN::EAST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_east;
                mblock->SetInitialGuess(n, value_east);
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

    trilinos.GetM() = solver.BuildPreconditioner(precond_pack,
                                                 precond_type,
                                                 param_precond,
                                                 trilinos.GetA());

    Belos::ReturnType ret = solver.Solve(solver_pack, solver_type,
                                         trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(),
                                         param_solver);
    bool is_converged = ret == Belos::ReturnType::Converged;
    EXPECT_TRUE(is_converged);
    double d_phi = (value_east - value_west) / g_rep.GetGlobalResolutionInternal()[0];
    for (LO node_l{0}; node_l < g_rep.GetNumberLocalCellsInternal(); node_l++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node_l);
        IndexGlobal ind = g_rep.MapOrdinalToIndexGlobalInternal(node_g);
        for (std::size_t n{0}; n < N; n++) {
            double phi = value_west + (ind.i() + 0.5) * d_phi;
            double val = trilinos.GetX()->getData()[node_l * N + n];
            EXPECT_NEAR(phi, val, 1e-8) << "Expected: " << phi << " and found: " << val
                                                     << " deviation is: " << std::fabs(phi - val);
        }
    }
}

TEST_F(IntegrationCartesianTrilinos1D, SolveStaggered) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(1);  // staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC, LO, GO> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC, LO, GO> solver;
    const double value_west{0.};
    const double value_east{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";
    const dare::Matrix::PreCondPackage precond_pack = dare::Matrix::PreCondPackage::MueLu;
    const std::string precond_type = "AMG";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    Teuchos::RCP<Teuchos::ParameterList> param_precond = Teuchos::rcp(new Teuchos::ParameterList());
    param_precond->set("problem: type", "MHD");  // works best in our cases
    param_precond->set("verbosity", "none");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        GO node_g = mblock->GetNode();
        bool is_left_edge{node_g == 0};
        bool is_right_edge{node_g == (g_rep.GetGlobalResolutionInternal()[0] - 1)};
        std::size_t stencil_size{3};
        if (is_left_edge)
            stencil_size = 1;
        else if (is_right_edge)
            stencil_size = 1;
        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            if (is_left_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_west;
                mblock->SetInitialGuess(n, value_west);
            } else if (is_right_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_east;
                mblock->SetInitialGuess(n, value_east);
            } else {
                mblock->template Get<CN::CENTER>(n) = 2.;
                mblock->template Get<CN::WEST>(n) = -1.;
                mblock->template Get<CN::EAST>(n) = -1.;
                mblock->GetRhs(n) = 0.;
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

    trilinos.GetM() = solver.BuildPreconditioner(precond_pack,
                                                 precond_type,
                                                 param_precond,
                                                 trilinos.GetA());

    Belos::ReturnType ret = solver.Solve(solver_pack, solver_type,
                                         trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(),
                                         param_solver);
    bool is_converged = ret == Belos::ReturnType::Converged;
    EXPECT_TRUE(is_converged);
    double d_phi = (value_east - value_west) / (g_rep.GetGlobalResolutionInternal()[0] - 1);
    for (LO node_l{0}; node_l < g_rep.GetNumberLocalCellsInternal(); node_l++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node_l);
        IndexGlobal ind = g_rep.MapOrdinalToIndexGlobalInternal(node_g);
        for (std::size_t n{0}; n < N; n++) {
            double phi = value_west + ind.i()* d_phi;
            double val = trilinos.GetX()->getData()[node_l * N + n];
            EXPECT_NEAR(phi, val, 1e-8) << "Expected: " << phi << " and found: " << val
                                        << " deviation is: " << std::fabs(phi - val);
        }
    }
}
