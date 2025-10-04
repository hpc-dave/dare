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

#include <iostream>

#include "AnalyticalSolutions/Diffusion.h"
#include "Data/DefaultTypes.h"
#include "Data/Field.h"
#include "Grid/Cartesian.h"
#include "IO/FileSystemManager.h"
#include "IO/VTKWriter.h"
#include "MPI/ExecutionManager.h"
#include "MatrixSystem/Trilinos.h"
#include "MatrixSystem/TrilinosSolver.h"
#include "ScopeGuard/ScopeGuard.h"

int main(int argc, char* argv[]) {
    using SC = dare::defaults::ScalarType;
    using GO = dare::defaults::GlobalOrdinalType;
    using LO = dare::defaults::LocalOrdinalType;
    using Grid = dare::Cartesian<1>;
    using GridVector = dare::GridVector<Grid, SC, 1>;
    using Field = dare::Field<Grid, SC, 1>;
    using Writer = dare::VTKWriter<Grid>;
    using IndexGlobal = typename IndexGlobal;
    using IndexLocal = typename Index;
    using VecSC = typename VecSC;
    using CNB = typename NeighborID;

    dare::ScopeGuard scope_guard(&argc, &argv);
    {
            GO nx{10};
            SC L{1};
            LO num_ghost = 2;
            int freq_write = 1;
            SC dt = 1e-3;

            IndexGlobal resolution_global(nx);
            VecSC size_global(L);

            dare::ExecutionManager exman;
            dare::FileSystemManager fman(&exman, "verification");
            fman.CheckWithUser(false);

            Grid grid("scalar_1D",
                      &exman,
                      resolution_global,
                      size_global,
                      num_ghost);
            IndexLocal scalar(0), staggered_x(1);  // grid is not staggered
            auto grid_s = grid.GetRepresentation(scalar);
            auto grid_x = grid.GetRepresentation(staggered_x);

            Field pressure("pressure", grid_s, 2);
            pressure.SetValues(1.);

            Field dp("dp", grid_s, 1);

            Field mom_x("momentum_x", grid_x, 2);
            mom_x.SetValues(0.);

            Field rho("density", grid_s, 2);
            rho.SetValues(0.25);

            Field beta("beta", grid_s, 1);
            beta.SetValues(1.);
            Field epsilon("porosity", grid_s, 2);
            epsilon.SetValues(0.5);
            // beta.GetDataVector().At(0, 0) = 0.;
            // beta.GetDataVector().At(1, 0) = 0.;
            // beta.GetDataVector().At(2, 0) = 0.;

            auto build_coef_x = [&](auto mblock) {
                // Normalize!
                using TimeScheme = dare::EULER_BACKWARD;
                using Divergence = dare::Divergence<Grid, TimeScheme>;
                using DDT = dare::DDT<Grid>;
                using FluxLimiter = dare::MINMOD;
                using TVD = dare::TVD<Grid, SC, FluxLimiter>;

                auto g_r = mblock->GetRepresentation();
                LO loc_o{mblock->GetLocalOrdinal()};
                IndexLocal ind{mblock->GetIndex()};
                const GridVector& u = mom_x.GetDataVector();

                Divergence div(*g_r, loc_o);
                DDT ddt(*g_r, loc_o, dt);
                TVD tvd(*g_r, loc_o, dare::Vector<1, const GridVector*>(&u));

                // time derivative
                (*mblock) = ddt(epsilon, rho, mom_x);

                // advection
                (*mblock) += div(epsilon, rho, tvd * u);

                // pressure force
                SC grad_p{0.};
                {
                    // This one here could be more generic by interpolating the pressure drop across all faces
                    // and treat this one here as an optimization for Cartesian grids
                    IndexLocal ind_nb(ind);
                    ind_nb.i() -= 1;
                    grad_p = pressure.GetDataVector(0).At(ind, 0) - pressure.GetDataVector().At(ind_nb, 0);
                    grad_p /= g_r->GetDistances().x();
                    grad_p *= g_r->GetCellVolume();
                }
                SC poro = dare::InterpolateToFace(*g_r, ind, CNB::WEST, epsilon.GetDataVector(), 0);

                mblock->GetRhs(0) -= poro * grad_p;

                // Apply BC
                IndexLocal ind_l = g_r->MapOrdinalToIndexLocalInternal(loc_o);
                IndexGlobal ind_g = g_r->MapLocalToGlobal(ind_l);
                if (ind_g.i() == 0) {
                    // WEST boundary
                    // Dirichlet
                    // std::cout << *mblock;
                    SC bc_0 = 1.;
                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->Get(0, 0, CNB::EAST) = 0.;
                    mblock->GetRhs(0) = bc_0;
                    mblock->Remove(0, 0, CNB::WEST);
                } else if (ind_g.i() == nx) {
                    // EAST boundary
                    // Neumann
                    mblock->Get(0, 0, CNB::WEST) = -1.;
                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->GetRhs(0) = 0;
                    mblock->Remove(0, 0, CNB::EAST);
                }
            };

            auto compute_defect = [&](const auto& grep, LO loc_o, IndexLocal ind) {
                using Divergence = dare::Divergence<Grid, dare::EULER_BACKWARD>;
                using FluxLimiter = dare::CDS;
                using TVD = dare::TVD<Grid, SC, FluxLimiter>;
                using FVStencil = dare::FaceValueStencil<Grid, SC, 1>;
                using CVStencil = dare::CenterValueStencil<Grid, SC, 1>;
                const GridVector& u = mom_x.GetDataVector();
                Divergence div(grep, loc_o);
                TVD tvd(grep, loc_o, dare::Vector<1, const GridVector*>(&u));

                CVStencil ONES;
                FVStencil uloc = dare::InterpolateToFaceStencil(grep, ind, u);
                ONES.SetAll(1.);
                SC defect = div(tvd.Interpolate(ONES, ONES) * uloc)[0];
                return defect;
            };

            auto build_coef_p = [&](auto mblock) {
                // Normalize
                // using Divergence = dare::Divergence<Grid, dare::EULER_BACKWARD>;
                // using FluxLimiter = dare::CDS;
                // using TVD = dare::TVD<Grid, SC, FluxLimiter>;
                // using FVStencil = dare::FaceValueStencil<Grid, SC, 1>;
                // using CVStencil = dare::CenterValueStencil<Grid, SC, 1>;

                auto g_r = mblock->GetRepresentation();
                LO loc_o{mblock->GetLocalOrdinal()};
                IndexLocal ind{mblock->GetIndex()};
                // const GridVector& u = mom_x.GetDataVector();

                // Divergence div(*g_r, loc_o);
                // TVD tvd(*g_r, loc_o, dare::Vector<1, const GridVector*>(&u));

                // subsitute with interpolation to face and then gradient
                SC rho_w = dare::InterpolateToFace(*g_r, ind, CNB::WEST, rho.GetDataVector(), 0);
                SC rho_e = dare::InterpolateToFace(*g_r, ind, CNB::EAST, rho.GetDataVector(), 0);
                SC beta_w = dare::InterpolateToFace(*g_r, ind, CNB::WEST, beta.GetDataVector(), 0);
                SC beta_e = dare::InterpolateToFace(*g_r, ind, CNB::EAST, beta.GetDataVector(), 0);
                SC poro_w = dare::InterpolateToFace(*g_r, ind, CNB::WEST, epsilon.GetDataVector(), 0);
                SC poro_e = dare::InterpolateToFace(*g_r, ind, CNB::EAST, epsilon.GetDataVector(), 0);

                mblock->Get(0, 0, CNB::WEST) = -poro_w * dt /
                            (poro_w * rho_w + beta_w * dt) * g_r->GetFaceArea().x() / g_r->GetDistances().x();
                mblock->Get(0, 0, CNB::EAST) = -poro_e * dt /
                            (poro_e * rho_e + beta_e * dt) * g_r->GetFaceArea().x() / g_r->GetDistances().x();

                mblock->Get(0, 0, CNB::CENTER) = -mblock->Get(0, 0, CNB::WEST) - mblock->Get(0, 0, CNB::EAST);

                // populate stencil, FOR HIGHER DIMENSIONS THAT SHOULD BE A DEDICATED SETUP
                // FVStencil u_close = dare::InterpolateToFaceStencil(*g_r, ind, u, 0);
                // FVStencil u_far = dare::InterpolateToFaceStencil(*g_r, ind, u, 1);
                // CVStencil ONES;
                // FVStencil uloc = dare::InterpolateToFaceStencil(*g_r, ind, u);
                // ONES.SetAll(1.);
                // SC defect = div(tvd.Interpolate(ONES, ONES) * uloc)[0];
                // SC u_w = dare::InterpolateToFace(*g_r, ind, CNB::WEST, u, 0);
                // SC u_e = dare::InterpolateToFace(*g_r, ind, CNB::EAST, u, 0);
                // std::cout << ind << " u_west: " << u_w << " u_east: " << u_e << std::endl;
                SC defect = compute_defect(*g_r, loc_o, ind);
                mblock->GetRhs(0) = -defect;
                // Apply BC
                IndexLocal ind_l = g_r->MapOrdinalToIndexLocalInternal(loc_o);
                IndexGlobal ind_g = g_r->MapLocalToGlobal(ind_l);

                if (ind_g.i() == 0) {
                    // WEST boundary
                    // Inflow
                    // std::cout << *mblock;
                    mblock->Get(0, 0, CNB::CENTER) += mblock->Get(0, 0, CNB::WEST);
                    mblock->Remove(0, 0, CNB::WEST);
                } else if (ind_g.i() == (nx - 1)) {
                    // EAST boundary
                    // prescribed pressure
                    SC bc_1 = 1.;

                    mblock->Get(0, 0, CNB::WEST) = 0.;
                    mblock->Get(0, 0, CNB::EAST) = 0.;
                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->GetRhs(0) = bc_1 - pressure.GetDataVector().At(ind, 0);
                    mblock->Remove(0, 0, CNB::EAST);
                }
            };

            dare::Trilinos<SC> msystem_x(&exman);
            dare::Trilinos<SC> msystem_p(&exman);
            Teuchos::RCP<Teuchos::ParameterList> p_ilu = Teuchos::rcp(new Teuchos::ParameterList());
            // parameters for ILU
            p_ilu->set("fact: drop tolerance", 1e-9);
            p_ilu->set("fact: level of fill", 1);
            p_ilu->set("schwarz: combine mode", "Add");
            Teuchos::RCP<Teuchos::ParameterList> p_amg = Teuchos::rcp(new Teuchos::ParameterList());
            // parameters for AMG
            p_amg->set("problem: type", "MHD");  // works best in our cases
            p_amg->set("verbosity", "none");
            Teuchos::RCP<Teuchos::ParameterList> p_solver = Teuchos::rcp(new Teuchos::ParameterList());
            p_solver->set("Convergence Tolerance", 1e-16);
            p_solver->set("Maximum Iterations", 5000);

            SC t_end = 1;
            SC time = 0.;
            int timestep = 0;

            {
                // in a scope, so it won't persist longer than required
                Writer writer(&exman, time, timestep);
                writer.Write(fman, &mom_x.GetDataVector(), &pressure.GetDataVector(), &dp.GetDataVector());
            }

            while (time < t_end) {
                timestep++;
                time += dt;

                mom_x.CopyDataVectorsToOldTimeStep();
                pressure.CopyDataVectorsToOldTimeStep();
                // dp.CopyDataVectorsToOldTimeStep();
                rho.CopyDataVectorsToOldTimeStep();
                epsilon.CopyDataVectorsToOldTimeStep();

                exman.Print(dare::Verbosity::Low) << "t: " << time << "\tstep: " << timestep << "\n";
                msystem_x.Build(grid_x, mom_x.GetDataVector(), build_coef_x, false);

                // msystem.PrintMatrix();
                // msystem.PrintX();
                // msystem.PrintB();

                dare::TrilinosSolver<SC> solver_x;
                msystem_x.GetM() = solver_x.BuildPreconditioner(dare::PreCondPackage::Ifpack2,
                                                            "ILUT",
                                                            p_ilu,
                                                            msystem_x.GetA());
                auto ret_x = solver_x.Solve(dare::SolverPackage::Belos,
                                          "BICGSTAB",
                                          msystem_x.GetM(),
                                          msystem_x.GetA(), msystem_x.GetX(), msystem_x.GetB(),
                                          p_solver);
                msystem_x.CopyTo(&mom_x.GetDataVector());
                if (ret_x != Belos::ReturnType::Converged) {
                    exman.Print(dare::Verbosity::Low) << "momentum did not converge\n";
                }

                mom_x.ExchangeHaloCells();

                // Update boundary cells
                mom_x.GetDataVector().At(num_ghost - 1, 0) = mom_x.GetDataVector().At(num_ghost, 0);
                mom_x.GetDataVector().At(nx + num_ghost + 1, 0) = mom_x.GetDataVector().At(nx + num_ghost, 0);

                exman.Print(dare::Verbosity::Low)
                    << "U: " << solver_x.GetNumIterations() << '\n';

                int it{0};
                for (; it < 10; it++) {
                    msystem_p.Build(grid_s, pressure.GetDataVector(), build_coef_p, false);

                    dare::TrilinosSolver<SC> solver_p;
                    msystem_p.GetM() = solver_p.BuildPreconditioner(dare::PreCondPackage::MueLu,
                                                                    "AMG",
                                                                    p_amg,
                                                                    msystem_p.GetA());

                    auto ret_p = solver_p.Solve(dare::SolverPackage::Belos,
                                                "BICGSTAB",
                                                msystem_p.GetM(),
                                                msystem_p.GetA(), msystem_p.GetX(), msystem_p.GetB(),
                                                p_solver);

                    msystem_p.CopyTo(&dp.GetDataVector());
                    msystem_p.AddTo(&pressure.GetDataVector());

                    if (ret_p != Belos::ReturnType::Converged) {
                        exman.Print(dare::Verbosity::Low) << "pressure did not converge\n";
                    }
                    exman.Print(dare::Verbosity::Low)
                        << "P: " << solver_p.GetNumIterations() << '\n';
                    dp.ExchangeHaloCells();
                    pressure.ExchangeHaloCells();

                    for (LO lo_int{1}; lo_int < grid_x.GetNumberLocalCellsInternal() - 1; lo_int++) {
                        LO lo = grid_x.MapInternalToLocal(lo_int);
                        IndexLocal ind = grid_x.MapOrdinalToIndexLocal(lo);
                        IndexLocal ind_prev(ind);

                        SC rho_x = dare::InterpolateToFace(grid_s, ind, CNB::WEST, rho.GetDataVector(), 0);
                        SC beta_x = dare::InterpolateToFace(grid_s, ind, CNB::WEST, beta.GetDataVector(), 0);
                        SC poro_x = dare::InterpolateToFace(grid_s, ind, CNB::WEST, epsilon.GetDataVector(), 0);
                        SC beta_xit = it == 0 ? beta_x : 0.;
                        ind_prev.i() -= 1;
                        SC dp_e = dp.GetDataVector().At(ind, 0);
                        SC dp_w = dp.GetDataVector().At(ind_prev, 0);
                        SC grad_dp = dp_e - dp_w;
                        grad_dp /= grid_s.GetDistances().x();
                        SC u_n = mom_x.GetDataVector().At(lo);
                        SC denominator_1 = 1. / (poro_x * rho_x + beta_xit * dt);
                        SC denominator_2 = 1. / (poro_x * rho_x + beta_x * dt);
                        SC f1 = denominator_1 * poro_x * rho_x * u_n;
                        SC f2 = -denominator_2 * poro_x * grad_dp * dt;
                        // SC delta_u = -dt / rho_x * grad_dp;
                        SC u_nn = f1 + f2;
                        mom_x.GetDataVector().At(lo) = u_nn;
                    }
                    // check for defect
                    SC defect_max{0.};
                    for (LO lo_int{0}; lo_int < grid_s.GetNumberLocalCellsInternal() - 1; lo_int++) {
                        LO lo = grid_x.MapInternalToLocal(lo_int);
                        IndexLocal ind = grid_x.MapOrdinalToIndexLocal(lo);
                        SC defect_loc = std::abs(compute_defect(pressure.GetGridRepresentation(), lo, ind));
                        defect_max = std::max(defect_max, defect_loc);
                    }
                    exman.Print(dare::Verbosity::Low) << "it: " << it << " Defect: " << defect_max << std::endl;
                    if (defect_max < 1e-14)
                        break;
                }
                // Update boundary cells
                mom_x.GetDataVector().At(num_ghost - 1, 0) = mom_x.GetDataVector().At(num_ghost, 0);
                mom_x.GetDataVector().At(nx + num_ghost, 0) = mom_x.GetDataVector().At(nx + num_ghost - 1, 0);
                mom_x.GetDataVector().At(nx + num_ghost + 1, 0) = mom_x.GetDataVector().At(nx + num_ghost, 0);
                if ((timestep % freq_write) == 0) {
                    Writer writer(&exman, time, timestep);
                    writer.Write(fman, &mom_x.GetDataVector(), &pressure.GetDataVector(), &dp.GetDataVector());
                }
            }
    }
    return 0;
}
