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

template <typename GridVector, typename GridRepresentation>
void ComputeAnalyticalSolution(GridVector* field, typename GridVector::DataType tau, const GridRepresentation& grep) {
    using LO = typename GridRepresentation::LocalOrdinalType;
    using SC = typename GridVector::DataType;
    using Index = typename GridRepresentation::Index;
    auto res = grep.GetLocalResolution();
    for (LO i{0}; i < res.i(); i++) {
        Index ind = grep.MapOrdinalToIndexLocal(i);
        SC zeta_b = grep.GetCoordinatesCenter(ind).x();
        SC zeta_f = 1 - zeta_b;
        if (zeta_f <= 0. || zeta_f > 1.)
            continue;

        SC phi_f = dare::analytical::Diffusion_N0D1(tau, zeta_f);
        SC phi_b = dare::analytical::Diffusion_N0D1(tau, zeta_b);
        field->At(i, 0) = phi_f;
        field->At(i, 1) = phi_b;
    }
}

int main(int argc, char* argv[]) {
    using SC = dare::defaults::ScalarType;
    using GO = dare::defaults::GlobalOrdinalType;
    using LO = dare::defaults::LocalOrdinalType;
    using Grid = dare::Grid::Cartesian<1>;
    using GridVector = dare::Data::GridVector<Grid, SC, 2>;
    using Field = dare::Data::Field<Grid, SC, 2>;
    using Writer = dare::io::VTKWriter<Grid>;
    using IndexGlobal = typename Grid::IndexGlobal;
    using IndexLocal = typename Grid::Index;
    using VecSC = typename Grid::VecSC;
    using CNB = typename Grid::NeighborID;

    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        GO nx = 10;
        SC L = 1;
        LO num_ghost = 2;
        int freq_write = 10;
        double dt = 1e-3;

        IndexGlobal resolution_global(nx);
        VecSC size_global(L);

        dare::mpi::ExecutionManager exman;
        dare::io::FileSystemManager fman(&exman, "verification");
        fman.CheckWithUser(false);

        if (exman.GetNumberProcesses() > 1) {
            exman.Terminate(__func__, "Trilinos has issues with distributed 1D, use the 2D test instead");
        }

        Grid grid("scalar",
                  &exman,
                  resolution_global,
                  size_global,
                  num_ghost);
        IndexLocal staggered(0);  // grid is not staggered
        auto grep = grid.GetRepresentation(staggered);

        Field field("concentration", grep, 2);
        field.SetComponentName(0, "forward");
        field.SetComponentName(1, "backward");
        field.SetValues(0.);

        GridVector analytical("analytical", grep);
        analytical.SetComponentName(0, "forward");
        analytical.SetComponentName(1, "backward");

        auto build_coef = [&](auto mblock) {
            using Gradient = dare::Matrix::Gradient<Grid>;
            using TimeScheme = dare::Matrix::EULER_BACKWARD;
            using Divergence = dare::Matrix::Divergence<Grid, TimeScheme>;
            using DDT = dare::Matrix::DDT<Grid>;

            auto g_r = mblock->GetRepresentation();
            LO loc_o{mblock->GetLocalOrdinal()};

            Gradient grad(*g_r, loc_o);
            Divergence div(*g_r, loc_o);
            DDT ddt(*g_r, loc_o, dt);

            (*mblock) = ddt(field);
            (*mblock) -= div(grad(*mblock));

            // Apply BC
            GO glob_o{mblock->GetGlobalOrdinal()};
            if (glob_o == 0) {
                // WEST boundary
                // Dirichlet for component 0
                // std::cout << *mblock;
                SC bc_0 = 1.;
                mblock->Get(0, CNB::CENTER) -= mblock->Get(0, CNB::WEST);
                mblock->GetRhs(0) -= 2. * bc_0 * mblock->Get(0, CNB::WEST);
                mblock->Remove(0, CNB::WEST);
                // Neumann for component 1
                mblock->Get(1, CNB::CENTER) += mblock->Get(1, CNB::WEST);
                mblock->Remove(1, CNB::WEST);
            } else if (glob_o == (nx - 1)) {
                // EAST boundary
                // Dirichlet for component 1
                SC bc_1 = 1.;
                mblock->Get(1, CNB::CENTER) -= mblock->Get(1, CNB::EAST);
                mblock->GetRhs(1) -= 2. * bc_1 * mblock->Get(1, CNB::EAST);
                mblock->Remove(1, CNB::EAST);
                // Neumann for component 0
                mblock->Get(0, CNB::CENTER) += mblock->Get(0, CNB::EAST);
                mblock->Remove(0, CNB::EAST);
            }
            // std::cout << *mblock;
        };

        dare::Matrix::Trilinos<SC> msystem(&exman);
        Teuchos::RCP<Teuchos::ParameterList> p_ilu = Teuchos::rcp(new Teuchos::ParameterList());
        // parameters for ILU
        p_ilu->set("fact: drop tolerance", 1e-9);
        p_ilu->set("fact: level of fill", 1);
        p_ilu->set("schwarz: combine mode", "Add");
        Teuchos::RCP<Teuchos::ParameterList> p_solver = Teuchos::rcp(new Teuchos::ParameterList());
        p_solver->set("Convergence Tolerance", 1e-14);
        p_solver->set("Maximum Iterations", 5000);

        SC t_end = 1;
        SC time = 0.;
        int timestep = 0;

        {
            // in a scope, so it won't persist longer than required
            Writer writer(&exman, time, timestep);
            writer.Write(fman, &field.GetDataVector(), &analytical);
        }

        while (time < t_end) {
            timestep++;
            time += dt;

            msystem.Build(grep, field.GetDataVector(), build_coef, false);

            dare::Matrix::TrilinosSolver<SC> solver;
            msystem.GetM() = solver.BuildPreconditioner(dare::Matrix::PreCondPackage::Ifpack2,
                                                        "ILUT",
                                                        p_ilu,
                                                        msystem.GetA());
            auto ret = solver.Solve(dare::Matrix::SolverPackage::Belos,
                                    "BICGSTAB",
                                    msystem.GetM(),
                                    msystem.GetA(), msystem.GetX(), msystem.GetB(),
                                    p_solver);
            msystem.CopyTo(&field.GetDataVector());
            if (ret != Belos::ReturnType::Converged) {
                exman.Terminate(__func__, "solver did not converge");
            }
            exman.Print(dare::mpi::Verbosity::Low)
                << "t: " << time << "\tstep: " << timestep << "\tit: " << solver.GetNumIterations() << '\n';

            // Verification
            ComputeAnalyticalSolution(&analytical, time, grep);

            SC err_f{0.}, err_b{0.};
            for (LO i{0}; i < grep.GetNumberLocalCellsInternal(); i++) {
                LO i_loc = grep.MapInternalToLocal(i);
                SC phi_af = analytical.At(i_loc, 0);
                SC phi_ab = analytical.At(i_loc, 1);
                SC phi_ff = field.GetDataVector().At(i_loc, 0);
                SC phi_fb = field.GetDataVector().At(i_loc, 1);
                err_f += phi_af - phi_ff;
                err_b += phi_ab - phi_fb;
            }
            err_f = exman.Allsum(err_f) / nx;
            err_b = exman.Allsum(err_b) / nx;
            exman(dare::mpi::Verbosity::Low)
                << "error (f | b): " << err_f << " | " << err_b << '\n';

            field.CopyDataVectorsToOldTimeStep();

            if ((timestep % freq_write) == 0) {
                Writer writer(&exman, time, timestep);
                writer.Write(fman, &field.GetDataVector(), &analytical);
            }
        }
    }
    return 0;
}
