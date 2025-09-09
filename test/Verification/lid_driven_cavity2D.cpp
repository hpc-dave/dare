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

using DataSet = std::vector<std::pair<double, double>>;
std::pair<DataSet, DataSet> GetLidDrivenCavityReference(int Re) {
    std::vector<double> y, u, x, v;
    x = std::vector<double>{1, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0};  // NOLINT
    y = std::vector<double>{1, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0};  // NOLINT
    switch (Re) {
    case 100:
        u = std::vector<double>{1, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.2109, -0.15662, -0.1015, -0.06424, -0.04775, -0.04192, -0.03717, 0};  // NOLINT
        v = std::vector<double>{0, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.1089, 0.10091, 0.09233, 0};   // NOLINT
        break;
    case 400:
        u = std::vector<double>{1, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256, 0.02135, -0.11477, -0.17119, -0.32726, -0.24299, -0.14612, -0.10338, -0.09266, -0.08186, 0};  // NOLINT
        v = std::vector<double>{0, -0.12146, -0.15663, -0.19254, -0.22847, -0.22827, -0.44993, -0.38598, 0.05186, 0.30174, 0.30203, 0.28124, 0.22965, 0.2092, 0.19713, 0.1836, 0};     // NOLINT
        break;
    case 1000:
        u = std::vector<double>{1, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.0608, -0.10648, -0.27805, -0.38289, -0.2973, -0.2222, -0.20196, -0.18109, 0};   // NOLINT
        v = std::vector<double>{0, -0.21388, -0.27669, -0.33714, -0.39188, -0.5155, -0.42665, -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0};  // NOLINT
        break;
    default:
        ERROR << "Cannot find reference data for Re=" << Re << ERROR_CLOSE;
    }

    DataSet yu, xv;
    for (std::size_t i{0}; i < y.size(); i++) {
        yu.push_back(std::make_pair(y[i], u[i]));
    }
    for (std::size_t i{0}; i < x.size(); i++) {
        yu.push_back(std::make_pair(x[i], v[i]));
    }
    return {yu, xv};
}

struct PDict {
    using density = double;
    using viscosity = double;
};

struct SDict {
    using tvd = dare::MINMOD;
    using viscous_stress = dare::PMDijkhuizenStressTensor;
    using time_scheme_convective = dare::EULER_BACKWARD;
    using momentum_iterations = dare::Newton;
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
        } else {
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
            if (ind_g.i() == extent.i() / 2 && ind_g.j() == extent.j() / 2) {
                for (auto face : g_r->GetFaces()) {
                    o->Get(0, 0, face) = 0.;
                }
                o->Get(0, 0, CNB::CENTER) = 1.;
                o->GetRhs(0) = 0;
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
        } else {
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
        SC Re = 100;
        LO num_ghost = 2;
        int num_tsteps = 2000;
        SC rho = 1.;
        SC utop = 1;
        SC mu = rho * utop * L / Re;
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

        pm.Initialize(&grid, &dt, BCPressure{pm.GetContinuity()}, BCMom<0, 1>{utop}, BCMom<1, 1>{});
        pm.GetParameterList()->set("momentum: Jacobian", "constant");
        pm.GetParameterList()->set("continuity: Jacobian", "constant");
        pm.SetDensity(rho);
        pm.SetViscosity(mu);
        using ObsType = typename ProjectionMethod::ObserverType;
        using State = typename ProjectionMethod::StateChange;

        ObsType o_pm([&](const ProjectionMethod& pm_inst, State state) {
            switch (state) {
            case State::SolvedContinuity: {
                // SC pref = pm.GetContinuity()->GetdP()->GetDataVector().At(p_fix, 0);
                // pm.GetContinuity()->GetdP()->GetDataVector() -= pref;
                // SC pref = pm.GetPressure()->GetDataVector().At(p_fix, 0);
                // pm.GetPressure()->GetDataVector() -= pref;
            } break;
            default: {
            }
            }
        });
        pm.Attach(&o_pm);
        auto cprop = pm.GetContinuity()->GetSolverNumericalProperties();
        // cprop.solver_type = "CG";
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

        auto [yu_ref, xv_ref] = GetLidDrivenCavityReference(static_cast<int>(Re));
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
