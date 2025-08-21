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

#include <concepts>
#include <iostream>
#include <type_traits>

#include "Algorithm/ConstantTimeStep.h"
#include "Algorithm/NavierStokes/ProjectionMethod.h"
#include "Algorithm/NavierStokes/ProjectionMethod_Cartesian.h"
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

using SC = dare::defaults::ScalarType;
using GO = dare::defaults::GlobalOrdinalType;
using LO = dare::defaults::LocalOrdinalType;
using Grid = dare::Cartesian<2>;
using GridVector = dare::GridVector<Grid, SC, 1>;
using Field = dare::Field<Grid, SC, 1>;
using Writer = dare::VTKWriter<Grid>;
using IndexGlobal = typename Grid::IndexGlobal;
using IndexLocal = typename Grid::Index;
using VecSC = typename Grid::VecSC;
using CNB = typename Grid::NeighborID;

struct PDict {
    using density = double;
    using viscosity = double;
};

struct SDict {
    using tvd = dare::MINMOD;
    using viscous_stress = dare::PMDijkhuizenStressTensor;
    using time_scheme_convective = dare::EULER_BACKWARD;
};

using BoundaryStrategy = dare::PMDefaultBoundaryType<Grid, SC>;
using ProjectionMethod = dare::ProjectionMethod<Grid, BoundaryStrategy, PDict, SDict>;

struct BCPressure {
public:
    using CRefType = std::unique_ptr<typename ProjectionMethod::ContinuityType>;
    explicit BCPressure(const CRefType& c) : continuity(c) {}

    /*!
     * @brief Applies the boundary condititions/boundary update for the dP field
     * @tparam T type of the input argument (local or global MatrixBlock or Field)
     * @param o the object for which the boundary condition is applied
     */
    template <typename T>
    void Apply(T* o) const {
        if constexpr (dare::FieldType<T>) {
            // ghost cell values
            IndexLocal extent_l = o->GetGridRepresentation().GetLocalResolutionInternal();
            IndexGlobal extent_g = o->GetGridRepresentation().GetGlobalResolutionInternal();
            IndexLocal indl_low;
            for (auto& i : indl_low)
                i = 0;
            IndexGlobal indg_low = o->GetGridRepresentation().MapLocalToGlobal(indl_low);
            IndexGlobal indg_up = o->GetGridRepresentation().MapLocalToGlobal(extent_l);
            if (indg_low.i() == 0) {
                // WEST - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(indl_low)};
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind_nb.i() -= 1;
                for (LO j{0}; j < extent_l.j(); j++) {
                    ind.j() = ind_orig.j() + j;
                    ind_nb.j() = ind_orig.j() + j;
                    o->GetDataVector(0).At(ind_nb, 0) = o->GetDataVector(0).At(ind, 0);
                }
            }
            if (indg_up.i() == extent_g.i()) {
                // EAST - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(extent_l)};
                ind_orig.j() = o->GetGridRepresentation().GetNumberGhostCells();
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind.i() -= 1;
                for (LO j{0}; j < extent_l.j(); j++) {
                    ind.j() = ind_orig.j() + j;
                    ind_nb.j() = ind_orig.j() + j;
                    o->GetDataVector(0).At(ind_nb, 0) = o->GetDataVector(0).At(ind, 0);
                }
            }
            if (indg_low.j() == 0) {
                // SOUTH - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(indl_low)};
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind_nb.j() -= 1;
                for (LO i{0}; i < extent_l.i(); i++) {
                    ind.i() = ind_orig.i() + i;
                    ind_nb.i() = ind_orig.i() + i;
                    o->GetDataVector(0).At(ind_nb, 0) = o->GetDataVector(0).At(ind, 0);
                }
            }
            if (indg_up.j() == extent_g.j()) {
                // NORTH - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(extent_l)};
                ind_orig.i() = o->GetGridRepresentation().GetNumberGhostCells();
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind.j() -= 1;
                for (LO i{0}; i < extent_l.i(); i++) {
                    ind.i() = ind_orig.i() + i;
                    ind_nb.i() = ind_orig.i() + i;
                    o->GetDataVector(0).At(ind_nb, 0) = o->GetDataVector(0).At(ind, 0);
                }
            }
        } else if constexpr (o->IsGlobal()) {
            auto g_r = o->GetRepresentation();
            LO loc_o = o->GetLocalOrdinal();
            IndexGlobal extent = g_r->GetGlobalResolutionInternal();
            IndexLocal ind_l = g_r->MapOrdinalToIndexLocalInternal(loc_o);
            IndexGlobal ind_g = g_r->MapLocalToGlobal(ind_l);
            // pressure
            if (ind_g.j() == 0) {
                // SOUTH - Wall
                o->Get(0, 0, CNB::CENTER) += o->Get(0, 0, CNB::SOUTH);
            } else if (ind_g.j() + 1 == extent.j()) {
                // NORTH - Wall
                o->Get(0, 0, CNB::CENTER) += o->Get(0, 0, CNB::NORTH);
            }
            if (ind_g.i() == 0) {
                // WEST - Inlet
                o->Get(0, 0, CNB::CENTER) += o->Get(0, 0, CNB::WEST);
            } else if (ind_g.i() + 1 == extent.i()) {
                // EAST - Wall
                o->Get(0, 0, CNB::CENTER) += o->Get(0, 0, CNB::EAST);
            }

            // Removal of coefficients
            if (ind_g.i() == 0)
                o->Remove(0, 0, CNB::WEST);
            if (ind_g.i() == (extent.i() - 1))
                o->Remove(0, 0, CNB::EAST);
            if (ind_g.j() == 0)
                o->Remove(0, 0, CNB::SOUTH);
            if (ind_g.j() == (extent.j() - 1))
                o->Remove(0, 0, CNB::NORTH);
        }
    }
    SC p_bc{0};
    const CRefType& continuity;  //!< reference to the continuity equation
};

template <int staggered, int order = 1>
struct BCMom {
public:
    SC ux_top;
    SC alpha{1. + static_cast<SC>(order == 2)};
    SC beta{static_cast<SC>(order == 2) / 3.};
    SC gamma{2. + static_cast<SC>(order == 2) * 2. / 3.};
    explicit BCMom(SC u_top = 0) : ux_top{u_top} {}

    template <typename T>
    void Apply(T* o) const {
        if constexpr (dare::FieldType<T>) {
            // ghost cell values
            IndexLocal extent_l = o->GetGridRepresentation().GetLocalResolutionInternal();
            IndexGlobal extent_g = o->GetGridRepresentation().GetGlobalResolutionInternal();
            IndexLocal indl_low;
            for (auto& i : indl_low)
                i = 0;
            IndexGlobal indg_low = o->GetGridRepresentation().MapLocalToGlobal(indl_low);
            IndexGlobal indg_up = o->GetGridRepresentation().MapLocalToGlobal(extent_l);
            if (indg_low.i() == 0) {
                // WEST - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(indl_low)};
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind.i() += static_cast<LO>(staggered == 0);
                ind_nb.i() -= 1;
                for (LO j{0}; j < extent_l.j(); j++) {
                    ind.j() = ind_orig.j() + j;
                    ind_nb.j() = ind_orig.j() + j;
                    o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                }
            }
            if (indg_up.i() == extent_g.i()) {
                // EAST - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(extent_l)};
                ind_orig.j() = o->GetGridRepresentation().GetNumberGhostCells();
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind.i() -= (1 + static_cast<LO>(staggered == 0));
                for (LO j{0}; j < extent_l.j(); j++) {
                    ind.j() = ind_orig.j() + j;
                    ind_nb.j() = ind_orig.j() + j;
                    o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                }
            }
            if (indg_low.j() == 0) {
                // SOUTH - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(indl_low)};
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                IndexLocal ind_far{ind};
                ind_nb.j() -= 1;
                ind.j() += static_cast<LO>(staggered == 1);
                for (LO i{0}; i < extent_l.i(); i++) {
                    ind.i() = ind_orig.i() + i;
                    ind_nb.i() = ind_orig.i() + i;
                    o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                }
            }
            if (indg_up.j() == extent_g.j()) {
                // NORTH - Wall
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(extent_l)};
                ind_orig.i() = o->GetGridRepresentation().GetNumberGhostCells();
                IndexLocal ind{ind_orig};
                ind_orig.j() += static_cast<LO>(staggered == 1);
                IndexLocal ind_nb{ind_orig};
                ind.j() -= 1;
                for (LO i{0}; i < extent_l.i(); i++) {
                    ind.i() = ind_orig.i() + i;
                    ind_nb.i() = ind_orig.i() + i;
                    if constexpr (staggered == 0) {
                        o->GetDataVector(0).At(ind_nb, 0) = 2. * ux_top - o->GetDataVector(0).At(ind, 0);
                    } else {
                        o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                    }
                }
            }
        } else if constexpr (o->IsGlobal()) {
            auto g_r = o->GetRepresentation();
            LO loc_o = o->GetLocalOrdinal();
            IndexGlobal extent = g_r->GetGlobalResolutionInternal();
            IndexLocal ind_l = g_r->MapOrdinalToIndexLocalInternal(loc_o);
            IndexGlobal ind_g = g_r->MapLocalToGlobal(ind_l);

            if constexpr (staggered == 0) {
                if (ind_g.j() == 0) {
                    // SOUTH - Wall
                    SC bc_wall{0.};
                    o->Get(0, 0, CNB::CENTER) -= alpha * o->Get(0, 0, CNB::SOUTH);
                    o->Get(0, 0, CNB::NORTH) += beta * o->Get(0, 0, CNB::SOUTH);
                    o->GetRhs(0) -= gamma * bc_wall * o->Get(0, 0, CNB::SOUTH);
                } else if (ind_g.j() + 1 == extent.j()) {
                    // NORTH - Wall
                    SC bc_wall{ux_top};
                    o->Get(0, 0, CNB::CENTER) -= alpha * o->Get(0, 0, CNB::NORTH);
                    o->Get(0, 0, CNB::SOUTH) += beta * o->Get(0, 0, CNB::NORTH);
                    o->GetRhs(0) -= gamma * bc_wall * o->Get(0, 0, CNB::NORTH);
                }
                if (ind_g.i() == 0) {
                    // WEST - Wall
                    // Dirichlet
                    for (CNB face : g_r->GetFaces())
                        o->Get(0, 0, face) = 0.;
                    o->Get(0, 0, CNB::CENTER) = 1.;
                    o->GetRhs(0) = 0.;
                } else if (ind_g.i() + 1 == extent.i()) {
                    // EAST
                    // Dirichlet
                    for (CNB face : g_r->GetFaces())
                        o->Get(0, 0, face) = 0.;
                    o->Get(0, 0, CNB::CENTER) = 1.;
                    o->GetRhs(0) = 0.;
                }
            } else if constexpr (staggered == 1) {
                if (ind_g.i() == 0) {
                    // WEST - Wall
                    // Dirichlet
                    SC bc_0 = 0.;
                    o->Get(0, 0, CNB::CENTER) -= alpha * o->Get(0, 0, CNB::WEST);
                    o->Get(0, 0, CNB::EAST) += beta * o->Get(0, 0, CNB::WEST);
                    o->GetRhs(0) -= gamma * bc_0 * o->Get(0, 0, CNB::WEST);
                } else if (ind_g.i() + 1 == extent.i()) {
                    // EAST - Wall
                    // Dirichlet
                    SC bc_0 = 0.;
                    o->Get(0, 0, CNB::CENTER) -= alpha * o->Get(0, 0, CNB::EAST);
                    o->Get(0, 0, CNB::WEST) += beta * o->Get(0, 0, CNB::EAST);
                    o->GetRhs(0) -= gamma * bc_0 * o->Get(0, 0, CNB::EAST);
                }
                if (ind_g.j() == 0) {
                    // SOUTH - Wall
                    SC bc_0 = 0.;
                    for (CNB face : g_r->GetFaces())
                        o->Get(0, 0, face) = 0.;
                    o->Get(0, 0, CNB::CENTER) = 1.;
                    o->GetRhs(0) = bc_0;
                } else if (ind_g.j() + 1 == extent.j()) {
                    // NORTH
                    // SOUTH - Wall
                    SC bc_0 = 0.;
                    for (CNB face : g_r->GetFaces())
                        o->Get(0, 0, face) = 0.;
                    o->Get(0, 0, CNB::CENTER) = 1.;
                    o->GetRhs(0) = bc_0;
                }
            }
            // Removal of coefficients
            if (ind_g.i() == 0)
                o->Remove(0, 0, CNB::WEST);
            if (ind_g.i() == (extent.i() - 1))
                o->Remove(0, 0, CNB::EAST);
            if (ind_g.j() == 0)
                o->Remove(0, 0, CNB::SOUTH);
            if (ind_g.j() == (extent.j() - 1))
                o->Remove(0, 0, CNB::NORTH);
        }
    }
};

int main(int argc, char* argv[]) {
    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        SC L{1}, H{1.};
        GO nx{128}, ny{128};
        SC Re = 20;
        LO num_ghost = 2;
        int num_tsteps = 2000;
        SC rho = 1000.;
        SC mu = 1e-3;
        SC utop = Re * mu / (rho * H);
        SC Co = 0.5;
        int freq_write = 100;
        dare::ConstantTimeStep<SC> dt{Co * L / nx / utop};
        SC sim_time = num_tsteps * dt;

        IndexGlobal resolution_global(nx, ny);
        VecSC size_global(L, H);

        dare::ExecutionManager exman;
        dare::FileSystemManager fman(&exman, "LidDrivenCavity2D");
        fman.CheckWithUser(false);

        Grid grid("Cartesian_2D",
                  &exman,
                  resolution_global,
                  size_global,
                  num_ghost);
        ProjectionMethod pm;

        IndexLocal scalar(0, 0), staggered_x(1, 0), staggered_y(0, 1);
        auto grid_s = grid.GetRepresentation(scalar);
        auto grid_x = grid.GetRepresentation(staggered_x);
        auto grid_y = grid.GetRepresentation(staggered_y);

        IndexLocal p_fix{nx / 2 + num_ghost, ny / 2 + num_ghost};
        pm.Initialize(&grid, &dt, BCPressure{pm.GetContinuity()}, BCMom<0, 2>{utop}, BCMom<1, 2>{});

        pm.SetMaxLoopIterations(1);
        pm.SetDensity(rho);
        pm.SetViscosity(mu);
        using ObsType = typename ProjectionMethod::ObserverType;
        using State = typename ProjectionMethod::StateChange;

        ObsType o_pm([&](const ProjectionMethod& pm_inst, State state) {
            if (state == State::PressureUpdated) {
                SC pref = pm.GetPressure()->GetDataVector().At(p_fix, 0);
                pm.GetPressure()->GetDataVector() -= pref;
            }
        });
        pm.Attach(&o_pm);
        auto cprop = pm.GetContinuity()->GetSolverNumericalProperties();
        cprop.solver_properties->set("Convergence Tolerance", dt / rho / L * nx * 1e-10);
        pm.GetContinuity()->SetSolverNumericalProperties(cprop);
        auto printer = [&]() {
            GridVector vec_u("u", grid_s);
            GridVector vec_v("v", grid_s);
            GridVector* u_staggered = &pm.GetMomentum(0)->GetField()->GetDataVector();
            GridVector* v_staggered = &pm.GetMomentum(1)->GetField()->GetDataVector();
            for (LO n{0}; n < grid_x.GetNumberLocalCellsInternal(); n++) {
                IndexLocal ind = grid_x.MapOrdinalToIndexLocalInternal(n);
                ind = grid_x.MapInternalToLocal(ind);
                IndexLocal ind_nb{ind};
                ind_nb.i() += 1;
                vec_u.At(ind, 0) = 0.5 * (u_staggered->At(ind_nb, 0) + u_staggered->At(ind, 0));
            }
            for (LO n{0}; n < grid_y.GetNumberLocalCellsInternal(); n++) {
                IndexLocal ind = grid_y.MapOrdinalToIndexLocalInternal(n);
                ind = grid_y.MapInternalToLocal(ind);
                IndexLocal ind_nb{ind};
                ind_nb.j() += 1;
                vec_v.At(ind, 0) = 0.5 * (v_staggered->At(ind_nb, 0) + v_staggered->At(ind, 0));
            }
            Writer writer(&exman, dt.GetTime(), dt.GetTimeStepCounter());
            writer.Write(fman,
                         &pm.GetContinuity()->GetField()->GetDataVector(),
                         &vec_u,
                         &vec_v,
                         &pm.GetMomentum(0)->GetField()->GetDataVector(),
                         &pm.GetMomentum(1)->GetField()->GetDataVector());
        };
        printer();
        // pm.SetTimer(dare::Timer{});
        dare::SetVerbosity(dare::Verbosity::Low);
        while (dt.GetTime() < sim_time) {
            pm.CopyToOld();
            dt.AdvanceTimeStep();
            pm.SolveFlowField();
            if (dt.GetTimeStepCounter() % freq_write == 0) {
                printer();
            }
        }
    }
    return 0;
}  // NOLINT
