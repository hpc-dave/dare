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
#include "Equations/FluxLimiter.h"
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
    // using GridVector = dare::Data::GridVector<Grid, SC, 2>;
    using Field = dare::Data::Field<Grid, SC, 2>;
    using Writer = dare::io::VTKWriter<Grid>;
    using IndexGlobal = typename Grid::IndexGlobal;
    using IndexLocal = typename Grid::Index;
    using VecLO = typename Grid::VecLO;
    using VecSC = typename Grid::VecSC;
    using VecSC2 = dare::utils::Vector<2, double>;
    using SCHEME = typename dare::Matrix::MINMOD;

    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        GO nx = 50;
        SC L = 1;
        LO num_ghost = 2;
        int freq_write = 10;
        double dt = 1e-3;
        VecLO periodic(1);
        VecSC velocity_0(1.), velocity_1(-1.);

        IndexGlobal resolution_global(nx);
        VecSC size_global(L);

        dare::mpi::ExecutionManager exman;
        dare::io::FileSystemManager fman(&exman, "dancing_waves");
        fman.CheckWithUser(false);

        if (exman.GetNumberProcesses() > 1) {
            exman.Terminate(__func__, "Trilinos has issues with distributed 1D, use the 2D test instead");
        }

        Grid grid("scalar",
                  &exman,
                  resolution_global,
                  size_global,
                  num_ghost,
                  periodic);
        IndexLocal staggered(0);  // grid is not staggered
        auto grep = grid.GetRepresentation(staggered);

        Field field("dancers", grep, 2);
        field.SetComponentName(0, "first");
        field.SetComponentName(1, "second");
        field.SetValues(0.);
        for (LO i{0}; i < grep.GetNumberLocalCells(); i++) {
            SC x = grep.GetCoordinatesCenter(IndexLocal(i)).x();
            if (x > 0.35 && x < 0.65) {
                field.GetDataVector().At(IndexLocal(i), 0) = 1.;
                field.GetDataVector().At(IndexLocal(i), 1) = 1.;
            }
        }
        field.CopyDataVectorsToOldTimeStep();

        VecSC2 mass_init(0., 0.);
        for (LO i{0}; i < grep.GetNumberLocalCellsInternal(); i++) {
            mass_init[0] += field.GetDataVector().At(IndexLocal(i), 0);
            mass_init[1] += field.GetDataVector().At(IndexLocal(i), 1);
        }
        mass_init[0] = exman.Allsum(mass_init[0]);
        mass_init[1] = exman.Allsum(mass_init[1]);

        auto build_coef = [&](auto mblock) {
            using TimeScheme = dare::Matrix::EULER_BACKWARD;
            using Divergence = dare::Matrix::Divergence<Grid, TimeScheme>;
            using TVD = dare::Matrix::TVD<Grid, SC, SCHEME>;
            using DDT = dare::Matrix::DDT<Grid>;

            auto g_r = mblock->GetRepresentation();
            LO loc_o{mblock->GetLocalOrdinal()};

            TVD u_0(*g_r, loc_o, velocity_0);
            TVD u_1(*g_r, loc_o, velocity_1);
            Divergence div(*g_r, loc_o);
            DDT ddt(*g_r, loc_o, dt);

            (*mblock) = ddt(field);
            // (*mblock) += div(u * field.GetDataVector());
            mblock->Add(0, div(u_0 * field.GetDataVector()).GetSlice(0));
            mblock->Add(1, div(u_1 * field.GetDataVector()).GetSlice(1));
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

        SC t_end = 4.;
        SC time = 0.;
        int timestep = 0;

        {
            // in a scope, so it won't persist longer than required
            Writer writer(&exman, time, timestep);
            writer.Write(fman, &field.GetDataVector());
        }

        while (time < t_end) {
            timestep++;
            time += dt;

            double base_velocity = static_cast<int>(time) % 2 == 1 ? -1. : 1.;
            velocity_0.x() = -base_velocity;
            velocity_1.x() = base_velocity;

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

            // exchange halo cells
            field.ExchangeHaloCells();

            VecSC2 mass(0., 0.);
            for (LO i{0}; i < grep.GetNumberLocalCellsInternal(); i++) {
                mass[0] += field.GetDataVector().At(IndexLocal(i), 0);
                mass[1] += field.GetDataVector().At(IndexLocal(i), 1);
            }
            mass[0] = exman.Allsum(mass[0]);
            mass[1] = exman.Allsum(mass[1]);

            VecSC2 err = (mass - mass_init) / mass_init;
            exman.Print(dare::mpi::Verbosity::Low)
                << "t: " << time << "\tstep: " << timestep << "\tit: " << solver.GetNumIterations()
                << "\tmass-error: " << err << '\n';

            field.CopyDataVectorsToOldTimeStep();

            if ((timestep % freq_write) == 0) {
                Writer writer(&exman, time, timestep);
                writer.Write(fman, &field.GetDataVector());
            }
        }
    }
    return 0;
}
