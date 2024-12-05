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

#include "ScopeGuard/ScopeGuard.h"
#include "Data/DefaultTypes.h"
#include "Data/Field.h"
#include "Grid/Cartesian.h"
#include "MatrixSystem/Trilinos.h"
#include "MatrixSystem/TrilinosSolver.h"
#include "IO/VTKWriter.h"
#include "MPI/ExecutionManager.h"
#include "IO/FileSystemManager.h"

int main(int argc, char* argv[]) {
    using SC = dare::defaults::ScalarType;
    using GO = dare::defaults::GlobalOrdinalType;
    using LO = dare::defaults::LocalOrdinalType;
    using Grid = dare::Grid::Cartesian<1>;
    using Field = dare::Data::Field<Grid, SC, 2>;
    using Writer = dare::io::VTKWriter<Grid>;
    using IndexGlobal = typename Grid::IndexGlobal;
    using IndexLocal = typename Grid::Index;
    using VecSC = typename Grid::VecSC;
    using CNB = typename Grid::NeighborID;

    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        GO nx = 100;
        SC L = 1;
        LO num_ghost = 2;
        int freq_write = 1;
        double dt = 1e-5;

        IndexGlobal resolution_global(nx);
        VecSC size_global(L);

        dare::mpi::ExecutionManager exman;
        dare::io::FileSystemManager fman(&exman, "verification");
        fman.CheckWithUser(false);

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

        auto build_coef = [&](auto mblock) {
            using Gradient = dare::Matrix::Gradient<Grid>;
            using Divergence = dare::Matrix::Divergence<Grid>;
            using TimeDerivative = dare::Matrix::EULER_BACKWARD;
            using DDT = dare::Matrix::DDT<Grid, TimeDerivative>;

            auto g_r = mblock->GetRepresentation();
            LO loc_o{mblock->GetLocalOrdinal()};

            Gradient grad(*g_r, loc_o);
            Divergence div(*g_r);
            DDT ddt(*g_r, loc_o, dt);

            (*mblock) = ddt(field) + div(grad(*mblock));

            // Apply BC
            GO glob_o{mblock->GetGlobalOrdinal()};
            if (glob_o == 0) {
                // WEST boundary
                // Dirichlet for component 0
                SC bc_0 = 1.;
                mblock->Get(0, CNB::CENTER) -= mblock->Get(0, CNB::WEST);
                mblock->GetRhs(0) -= 2. * bc_0 * mblock->Get(0, CNB::WEST);
                mblock->Remove(0, CNB::WEST);
                // Neumann for component 1
                mblock->Get(1, CNB::CENTER) += mblock->Get(1, CNB::WEST);
                mblock->Remove(1, CNB::WEST);
            } else if (glob_o == (nx-1)) {
                // EAST boundary
                // Dirichlet for component 1
                SC bc_1 = 1.;
                mblock->Get(1, CNB::CENTER) -= mblock->Get(1, CNB::WEST);
                mblock->GetRhs(1) -= 2. * bc_1 * mblock->Get(1, CNB::WEST);
                mblock->Remove(1, CNB::WEST);
                // Neumann for component 0
                mblock->Get(0, CNB::CENTER) += mblock->Get(0, CNB::WEST);
                mblock->Remove(0, CNB::WEST);
            }
        };

        dare::Matrix::Trilinos<SC> msystem(&exman);
        Teuchos::RCP<Teuchos::ParameterList> p_ilu = Teuchos::rcp(new Teuchos::ParameterList());
        // paramters for ILU
        p_ilu->set("fact: drop tolerance", 1e-9);
        p_ilu->set("fact: level of fill", 1);
        p_ilu->set("schwarz: combine mode", "Add");
        Teuchos::RCP<Teuchos::ParameterList> p_solver = Teuchos::rcp(new Teuchos::ParameterList());
        p_solver->set("Convergence Tolerance", 1e-14);
        p_solver->set("Maximum Iterations", 5000);

        double t_end = 1;
        double time = 0.;
        int timestep = 0;

        {
            // in a scope, so it won't persist longer than required
            Writer writer(&exman, time, timestep);
            writer.Write(fman, &field.GetDataVector());
            }

            while (time < t_end) {
                timestep++;
                time += dt;

                msystem.Build(grep, field.GetDataVector(), build_coef);
                dare::Matrix::TrilinosSolver<SC> solver;
                msystem.GetM() = solver.BuildPreconditioner(dare::Matrix::PreCondPackage::Ifpack2,
                                                    "RILUK",
                                                    p_ilu,
                                                    msystem.GetA());
                auto ret = solver.Solve(dare::Matrix::SolverPackage::Belos,
                                   "BICGSTAB",
                                   msystem.GetM(), msystem.GetA(), msystem.GetX(), msystem.GetB(),
                                   p_solver);
                if (ret != Belos::ReturnType::Converged) {
                    exman.Terminate(__func__, "solver did not converge");
                }

                field.CopyDataVectorsToOldTimeStep();

                if (timestep % freq_write) {
                    Writer writer(&exman, time, timestep);
                    writer.Write(fman, &field.GetDataVector());
                }
            }
        }
}
