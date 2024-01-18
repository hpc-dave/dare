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

#include <gtest/gtest.h>

#include "Grid/Cartesian.h"
#include "Grid/DefaultTypes.h"
#include "MatrixSystem/Trilinos.h"
#include "MatrixSystem/TrilinosSolver.h"
#include "Utilities/Vector.h"

namespace dare::test {

template <std::size_t Dim, typename GO>
dare::utils::Vector<Dim, GO> GetResolutionIntegrationTestCartesianTrilinos() {
    dare::utils::Vector<Dim, GO> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 20 + n;
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
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_west{0.};
    const double value_east{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "high");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        GO node_g = mblock->GetNode();
        bool is_west_edge{node_g == 0};
        bool is_east_edge{node_g == (g_rep.GetGlobalResolutionInternal()[0] - 1)};
        std::size_t stencil_size{3};
        if (is_west_edge)
            stencil_size--;
        else if (is_east_edge)
            stencil_size--;
        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            mblock->GetRhs(n) = 0.;
            mblock->template Get<CN::CENTER>(n) = 2.;
            if (!is_west_edge) {
                mblock->template Get<CN::WEST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_west;
                mblock->SetInitialGuess(n, value_west);
            }

            if (!is_east_edge) {
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
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_west{0.};
    const double value_east{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "low");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        GO node_g = mblock->GetNode();
        bool is_west_edge{node_g == 0};
        bool is_east_edge{node_g == (g_rep.GetGlobalResolutionInternal()[0] - 1)};
        std::size_t stencil_size{3};
        if (is_west_edge)
            stencil_size = 1;
        else if (is_east_edge)
            stencil_size = 1;
        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            if (is_west_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_west;
                mblock->SetInitialGuess(n, value_west);
            } else if (is_east_edge) {
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

TEST_F(IntegrationCartesianTrilinos2D, SolveScalarX) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(0, 0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_west{0.};
    const double value_east{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "low");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        IndexGlobal ind = mblock->GetIndexInternal();
        bool is_west_edge{ind.i() == 0};
        bool is_east_edge{ind.i() == (g_rep.GetGlobalResolutionInternal().i() - 1)};
        bool is_south_edge{ind.j() == 0};
        bool is_north_edge{ind.j() == (g_rep.GetGlobalResolutionInternal().j() - 1)};

        std::size_t stencil_size{5};
        stencil_size -= is_west_edge;
        stencil_size -= is_east_edge;
        stencil_size -= is_south_edge;
        stencil_size -= is_north_edge;

        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            mblock->GetRhs(n) = 0.;
            mblock->template Get<CN::CENTER>(n) = 4.;

            // Dirichlet condition
            if (!is_west_edge) {
                mblock->template Get<CN::WEST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_west;
                mblock->SetInitialGuess(n, value_west);
            }

            // Dirichlet condition
            if (!is_east_edge) {
                mblock->template Get<CN::EAST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_east;
                mblock->SetInitialGuess(n, value_east);
            }

            // Neumann condition
            if (!is_south_edge) {
                mblock->template Get<CN::SOUTH>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }

            // Neumann condition
            if (!is_north_edge) {
                mblock->template Get<CN::NORTH>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

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

TEST_F(IntegrationCartesianTrilinos2D, SolveScalarY) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(0, 0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_south{0.};
    const double value_north{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "low");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        IndexGlobal ind = mblock->GetIndexInternal();
        bool is_west_edge{ind.i() == 0};
        bool is_east_edge{ind.i() == (g_rep.GetGlobalResolutionInternal().i() - 1)};
        bool is_south_edge{ind.j() == 0};
        bool is_north_edge{ind.j() == (g_rep.GetGlobalResolutionInternal().j() - 1)};

        std::size_t stencil_size{5};
        stencil_size -= is_west_edge;
        stencil_size -= is_east_edge;
        stencil_size -= is_south_edge;
        stencil_size -= is_north_edge;

        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            mblock->GetRhs(n) = 0.;
            mblock->template Get<CN::CENTER>(n) = 4.;

            // Dirichlet condition
            if (!is_south_edge) {
                mblock->template Get<CN::SOUTH>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_south;
                mblock->SetInitialGuess(n, value_south);
            }

            // Dirichlet condition
            if (!is_north_edge) {
                mblock->template Get<CN::NORTH>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_north;
                mblock->SetInitialGuess(n, value_north);
            }

            // Neumann condition
            if (!is_west_edge) {
                mblock->template Get<CN::WEST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }

            // Neumann condition
            if (!is_east_edge) {
                mblock->template Get<CN::EAST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

    Belos::ReturnType ret = solver.Solve(solver_pack, solver_type,
                                         trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(),
                                         param_solver);
    bool is_converged = ret == Belos::ReturnType::Converged;
    EXPECT_TRUE(is_converged);
    double d_phi = (value_north - value_south) / g_rep.GetGlobalResolutionInternal().j();
    for (LO node_l{0}; node_l < g_rep.GetNumberLocalCellsInternal(); node_l++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node_l);
        IndexGlobal ind = g_rep.MapOrdinalToIndexGlobalInternal(node_g);
        for (std::size_t n{0}; n < N; n++) {
            double phi = value_south + (ind.j() + 0.5) * d_phi;
            double val = trilinos.GetX()->getData()[node_l * N + n];
            EXPECT_NEAR(phi, val, 1e-8) << "Expected: " << phi << " and found: " << val
                                        << " deviation is: " << std::fabs(phi - val);
        }
    }
}

TEST_F(IntegrationCartesianTrilinos2D, SolveStaggeredX) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(1, 0);  // staggered in X direction
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_west{0.};
    const double value_east{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "low");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        IndexGlobal ind = mblock->GetIndexInternal();
        bool is_west_edge{ind.i() == 0};
        bool is_east_edge{ind.i() == (g_rep.GetGlobalResolutionInternal().i() - 1)};
        bool is_south_edge{ind.j() == 0};
        bool is_north_edge{ind.j() == (g_rep.GetGlobalResolutionInternal().j() - 1)};
        std::size_t stencil_size{5};
        if (is_west_edge || is_east_edge) {
            stencil_size = 1;
        } else {
            stencil_size -= is_south_edge;
            stencil_size -= is_north_edge;
        }
        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            if (is_west_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_west;
                mblock->SetInitialGuess(n, value_west);
            } else if (is_east_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_east;
                mblock->SetInitialGuess(n, value_east);
            } else {
                if (is_south_edge) {
                    mblock->template Get<CN::CENTER>(n) = 3.;
                    mblock->template Get<CN::WEST>(n) = -1.;
                    mblock->template Get<CN::EAST>(n) = -1.;
                    mblock->template Get<CN::NORTH>(n) = -1.;
                } else if (is_north_edge) {
                    mblock->template Get<CN::CENTER>(n) = 3.;
                    mblock->template Get<CN::WEST>(n) = -1.;
                    mblock->template Get<CN::EAST>(n) = -1.;
                    mblock->template Get<CN::SOUTH>(n) = -1.;
                } else {
                    mblock->template Get<CN::CENTER>(n) = 4.;
                    mblock->template Get<CN::WEST>(n) = -1.;
                    mblock->template Get<CN::EAST>(n) = -1.;
                    mblock->template Get<CN::SOUTH>(n) = -1.;
                    mblock->template Get<CN::NORTH>(n) = -1.;
                }
                mblock->GetRhs(n) = 0.;
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

    Belos::ReturnType ret = solver.Solve(solver_pack, solver_type,
                                         trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(),
                                         param_solver);
    bool is_converged = ret == Belos::ReturnType::Converged;
    EXPECT_TRUE(is_converged);
    double d_phi = (value_east - value_west) / (g_rep.GetGlobalResolutionInternal().i() - 1);
    for (LO node_l{0}; node_l < g_rep.GetNumberLocalCellsInternal(); node_l++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node_l);
        IndexGlobal ind = g_rep.MapOrdinalToIndexGlobalInternal(node_g);
        for (std::size_t n{0}; n < N; n++) {
            double phi = value_west + ind.i() * d_phi;
            double val = trilinos.GetX()->getData()[node_l * N + n];
            EXPECT_NEAR(phi, val, 1e-8) << "Expected: " << phi << " and found: " << val
                                        << " deviation is: " << std::fabs(phi - val);
        }
    }
}

TEST_F(IntegrationCartesianTrilinos2D, SolveStaggeredY) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(0, 1);  // staggered in Y direction
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_south{0.};
    const double value_north{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "low");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        IndexGlobal ind = mblock->GetIndexInternal();
        bool is_west_edge{ind.i() == 0};
        bool is_east_edge{ind.i() == (g_rep.GetGlobalResolutionInternal().i() - 1)};
        bool is_south_edge{ind.j() == 0};
        bool is_north_edge{ind.j() == (g_rep.GetGlobalResolutionInternal().j() - 1)};
        std::size_t stencil_size{5};
        if (is_south_edge || is_north_edge) {
            stencil_size = 1;
        } else {
            stencil_size -= is_west_edge;
            stencil_size -= is_east_edge;
        }
        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            if (is_south_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_south;
                mblock->SetInitialGuess(n, value_south);
            } else if (is_north_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_north;
                mblock->SetInitialGuess(n, value_north);
            } else {
                if (is_west_edge) {
                    mblock->template Get<CN::CENTER>(n) = 3.;
                    mblock->template Get<CN::EAST>(n) = -1.;
                    mblock->template Get<CN::SOUTH>(n) = -1.;
                    mblock->template Get<CN::NORTH>(n) = -1.;
                } else if (is_east_edge) {
                    mblock->template Get<CN::CENTER>(n) = 3.;
                    mblock->template Get<CN::WEST>(n) = -1.;
                    mblock->template Get<CN::SOUTH>(n) = -1.;
                    mblock->template Get<CN::NORTH>(n) = -1.;
                } else {
                    mblock->template Get<CN::CENTER>(n) = 4.;
                    mblock->template Get<CN::WEST>(n) = -1.;
                    mblock->template Get<CN::EAST>(n) = -1.;
                    mblock->template Get<CN::SOUTH>(n) = -1.;
                    mblock->template Get<CN::NORTH>(n) = -1.;
                }
                mblock->GetRhs(n) = 0.;
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

    Belos::ReturnType ret = solver.Solve(solver_pack, solver_type,
                                         trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(),
                                         param_solver);
    bool is_converged = ret == Belos::ReturnType::Converged;
    EXPECT_TRUE(is_converged);
    double d_phi = (value_north - value_south) / (g_rep.GetGlobalResolutionInternal().j() - 1);
    for (LO node_l{0}; node_l < g_rep.GetNumberLocalCellsInternal(); node_l++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node_l);
        IndexGlobal ind = g_rep.MapOrdinalToIndexGlobalInternal(node_g);
        for (std::size_t n{0}; n < N; n++) {
            double phi = value_south + ind.j() * d_phi;
            double val = trilinos.GetX()->getData()[node_l * N + n];
            EXPECT_NEAR(phi, val, 1e-8) << "Expected: " << phi << " and found: " << val
                                        << " deviation is: " << std::fabs(phi - val);
        }
    }
}

TEST_F(IntegrationCartesianTrilinos3D, SolveScalarX) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(0, 0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_west{0.};
    const double value_east{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "high");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        IndexGlobal ind = mblock->GetIndexInternal();
        bool is_west_edge{ind.i() == 0};
        bool is_east_edge{ind.i() == (g_rep.GetGlobalResolutionInternal().i() - 1)};
        bool is_south_edge{ind.j() == 0};
        bool is_north_edge{ind.j() == (g_rep.GetGlobalResolutionInternal().j() - 1)};
        bool is_bottom_edge{ind.k() == 0};
        bool is_top_edge{ind.k() == (g_rep.GetGlobalResolutionInternal().k() - 1)};

        std::size_t stencil_size{7};
        stencil_size -= is_west_edge;
        stencil_size -= is_east_edge;
        stencil_size -= is_south_edge;
        stencil_size -= is_north_edge;
        stencil_size -= is_bottom_edge;
        stencil_size -= is_top_edge;

        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            mblock->GetRhs(n) = 0.;
            mblock->template Get<CN::CENTER>(n) = 6.;

            // Dirichlet condition
            if (!is_west_edge) {
                mblock->template Get<CN::WEST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_west;
                mblock->SetInitialGuess(n, value_west);
            }

            // Dirichlet condition
            if (!is_east_edge) {
                mblock->template Get<CN::EAST>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) += 1.;
                mblock->GetRhs(n) = 2. * value_east;
                mblock->SetInitialGuess(n, value_east);
            }

            // Neumann condition
            if (!is_south_edge) {
                mblock->template Get<CN::SOUTH>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }

            // Neumann condition
            if (!is_north_edge) {
                mblock->template Get<CN::NORTH>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }

            // Neumann condition
            if (!is_bottom_edge) {
                mblock->template Get<CN::BOTTOM>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }

            // Neumann condition
            if (!is_top_edge) {
                mblock->template Get<CN::TOP>(n) = -1.;
            } else {
                mblock->template Get<CN::CENTER>(n) -= 1.;
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

    Belos::ReturnType ret = solver.Solve(solver_pack, solver_type,
                                         trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(),
                                         param_solver);

    bool is_converged = ret == Belos::ReturnType::Converged;
    EXPECT_TRUE(is_converged);
    double d_phi = (value_east - value_west) / g_rep.GetGlobalResolutionInternal().i();
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

// We only test the configuration in z-direction to save some time in the unit testing
TEST_F(IntegrationCartesianTrilinos3D, SolveStaggeredZ) {
    using CN = dare::Matrix::CartesianNeighbor;
    GridType::Options opt(1, 0);  // staggered in X direction
    auto g_rep = grid->GetRepresentation(opt);
    GridVector data("test", g_rep);
    dare::Matrix::Trilinos<SC> trilinos(&exec_man);
    dare::Matrix::TrilinosSolver<SC> solver;
    const double value_bottom{0.};
    const double value_top{1.};

    const dare::Matrix::SolverPackage solver_pack = dare::Matrix::SolverPackage::BumbleBee;
    const std::string solver_type = "BICGSTAB2";

    Teuchos::RCP<Teuchos::ParameterList> param_solver = Teuchos::rcp(new Teuchos::ParameterList());
    param_solver->set("Convergence Tolerance", 1e-14);
    param_solver->set("Maximum Iterations", 5000);
    param_solver->set("Verbosity", "low");

    auto functor = [&](auto mblock) {
        EXPECT_TRUE(mblock->IsGlobal());
        IndexGlobal ind = mblock->GetIndexInternal();
        bool is_west_edge{ind.i() == 0};
        bool is_east_edge{ind.i() == (g_rep.GetGlobalResolutionInternal().i() - 1)};
        bool is_south_edge{ind.j() == 0};
        bool is_north_edge{ind.j() == (g_rep.GetGlobalResolutionInternal().j() - 1)};
        bool is_bottom_edge{ind.k() == 0};
        bool is_top_edge{ind.k() == (g_rep.GetGlobalResolutionInternal().k() - 1)};
        std::size_t stencil_size{7};
        if (is_bottom_edge || is_top_edge) {
            stencil_size = 1;
        } else {
            stencil_size -= is_west_edge;
            stencil_size -= is_east_edge;
            stencil_size -= is_south_edge;
            stencil_size -= is_north_edge;
        }
        for (std::size_t n{0}; n < mblock->GetNumComponents(); n++) {
            mblock->Resize(n, stencil_size);
            if (is_bottom_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_bottom;
                mblock->SetInitialGuess(n, value_bottom);
            } else if (is_top_edge) {
                mblock->template Get<CN::CENTER>(n) = 1.;
                mblock->GetRhs(n) = value_top;
                mblock->SetInitialGuess(n, value_top);
            } else {
                mblock->template Get<CN::CENTER>(n) = 6.;
                if (is_west_edge)
                    mblock->template Get<CN::CENTER>(n) -= 1.;
                else
                    mblock->template Get<CN::WEST>(n) = -1.;

                if (is_east_edge)
                    mblock->template Get<CN::CENTER>(n) -= 1.;
                else
                    mblock->template Get<CN::EAST>(n) = -1.;

                if (is_south_edge)
                    mblock->template Get<CN::CENTER>(n) -= 1.;
                else
                    mblock->template Get<CN::SOUTH>(n) = -1.;

                if (is_north_edge)
                    mblock->template Get<CN::CENTER>(n) -= 1.;
                else
                    mblock->template Get<CN::NORTH>(n) = -1.;

                mblock->template Get<CN::BOTTOM>(n) = -1.;
                mblock->template Get<CN::TOP>(n) = -1.;
                mblock->GetRhs(n) = 0.;
            }
        }
        mblock->Finalize();
    };

    trilinos.Build(g_rep, data, functor, false);

    Belos::ReturnType ret = solver.Solve(solver_pack, solver_type,
                                         trilinos.GetM(),
                                         trilinos.GetA(), trilinos.GetX(), trilinos.GetB(),
                                         param_solver);

    bool is_converged = ret == Belos::ReturnType::Converged;
    EXPECT_TRUE(is_converged);
    double d_phi = (value_top - value_bottom) / (g_rep.GetGlobalResolutionInternal().k() - 1);
    for (LO node_l{0}; node_l < g_rep.GetNumberLocalCellsInternal(); node_l++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node_l);
        IndexGlobal ind = g_rep.MapOrdinalToIndexGlobalInternal(node_g);
        for (std::size_t n{0}; n < N; n++) {
            double phi = value_bottom + ind.k() * d_phi;
            double val = trilinos.GetX()->getData()[node_l * N + n];
            EXPECT_NEAR(phi, val, 1e-8) << "Expected: " << phi << " and found: " << val
                                        << " deviation is: " << std::fabs(phi - val);
        }
    }
}
