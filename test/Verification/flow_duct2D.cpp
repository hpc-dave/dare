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
#include <limits>
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

/*!
 * @brief boundary conditions for the continuity equation
 */
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
                // WEST - Inlet
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(indl_low)};
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind_nb.i() -= 1;
                for (LO j{0}; j < extent_l.j(); j++) {
                    ind.j() = ind_orig.j() + j;
                    ind_nb.j() = ind_orig.j() + j;
                    SC dp_c{o->GetDataVector(0).At(ind, 0)};
                    o->GetDataVector(0).At(ind_nb, 0) = dp_c;
                }
            }
            if (indg_up.i() == extent_g.i()) {
                // EAST - Outlet
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(extent_l)};
                ind_orig.j() = o->GetGridRepresentation().GetNumberGhostCells();
                IndexLocal ind{ind_orig};
                IndexLocal ind_e{ind};
                ind.i() -= 1;
                for (LO j{0}; j < extent_l.j(); j++) {
                    ind.j() = ind_orig.j() + j;
                    ind_e.j() = ind_orig.j() + j;
                    SC dp_c = o->GetDataVector(0).At(ind, 0);
                    SC p_c = continuity->GetPressure()->GetDataVector().At(ind, 0);
                    SC p_e = continuity->GetPressure()->GetDataVector().At(ind_e, 0);
                    p_c += dp_c;
                    SC p_e_new = 2. * p_bc - p_c;
                    SC dp_e = p_e_new - p_e;
                    o->GetDataVector(0).At(ind_e, 0) = dp_e;
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
                // EAST - Outlet
                // prescribed pressure

                SC bc_1 = p_bc;
                IndexLocal ind_c = g_r->MapInternalToLocal(ind_l);

                // for (CNB face : g_r->GetFaces())
                //     o->Get(0, 0, face) = 0.;

                // o->Get(0, 0, CNB::CENTER) = 1.;
                // o->GetRhs(0) = bc_1 - continuity->GetPressure()->GetDataVector().At(ind_c, 0);

                IndexLocal ind_w{ind_c};
                ind_w.i() -= 1;
                SC p_c = continuity->GetPressure()->GetDataVector().At(ind_c, 0);
                SC p_w = continuity->GetPressure()->GetDataVector().At(ind_w, 0);
                SC delta_p = p_c - p_w;
                SC dp_e = (bc_1 + 0.5 * delta_p) - (p_c + delta_p);
                SC alpha_e = o->Get(0, 0, CNB::EAST);
                o->Get(0, 0, CNB::EAST) = 0.;
                o->Get(0, 0, CNB::CENTER) -= alpha_e;
                o->GetRhs(0) -= alpha_e * dp_e;
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

/*!
 * @brief boundary conditions for the momentum equation
 * @tparam staggered indicator which dimension is staggered (0=x, 1=y)
 * @tparam order order of application (1=first order, 2=second order)
 */
template <int staggered, int order = 1>
struct BCMom {
public:
    SC ux_in;
    SC alpha{1. + static_cast<SC>(order == 2)};
    SC beta{static_cast<SC>(order == 2) / 3.};
    SC gamma{2. + static_cast<SC>(order == 2) * 2. / 3.};
    explicit BCMom(SC uin = 0) : ux_in{uin} {}

    /*!
     * @brief core function for applying boundary conditions and ghost cell updates
     * @tparam T type of incoming object
     * @param o object to adjust
     * the object can be a matrix block or a GridVector
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
                // WEST - Inlet
                IndexLocal ind_orig{o->GetGridRepresentation().MapInternalToLocal(indl_low)};
                IndexLocal ind{ind_orig};
                IndexLocal ind_nb{ind};
                ind_nb.i() -= 1;
                for (LO j{0}; j < extent_l.j(); j++) {
                    ind.j() = ind_orig.j() + j;
                    ind_nb.j() = ind_orig.j() + j;
                    if constexpr (staggered == 0)
                        o->GetDataVector(0).At(ind_nb, 0) = ux_in;
                    else
                        o->GetDataVector(0).At(ind_nb, 0) = 0.;
                }
            }
            if (indg_up.i() == extent_g.i()) {
                // EAST - Outlet
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
                IndexLocal ind_far{ind};
                ind_nb.j() -= 1;
                ind_far.j() += 1;
                ind.j() += static_cast<LO>(staggered == 1);
                for (LO i{0}; i < extent_l.i(); i++) {
                    ind.i() = ind_orig.i() + i;
                    ind_nb.i() = ind_orig.i() + i;
                    // o->GetDataVector(0).At(ind_nb, 0) = 0.;
                    if constexpr (staggered == 0) {
                        o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                        // o->GetDataVector(0).At(ind_nb, 0) += beta * o->GetDataVector(0).At(ind_far, 0);
                        // o->GetDataVector(0).At(ind_nb, 0) += gamma * 0.;
                    } else {
                        o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                    }
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
                    // o->GetDataVector(0).At(ind_nb, 0) = 0.;
                    if constexpr (staggered == 0) {
                        o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                    } else {
                        o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
                    }
                    // o->GetDataVector(0).At(ind_nb, 0) = -o->GetDataVector(0).At(ind, 0);
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
                    // NORTH
                    SC bc_wall{0.};
                    o->Get(0, 0, CNB::CENTER) -= alpha * o->Get(0, 0, CNB::NORTH);
                    o->Get(0, 0, CNB::SOUTH) += beta * o->Get(0, 0, CNB::NORTH);
                    o->GetRhs(0) -= gamma * bc_wall * o->Get(0, 0, CNB::NORTH);
                }
                if (ind_g.i() == 0) {
                    // WEST
                    // Dirichlet
                    for (CNB face : g_r->GetFaces())
                        o->Get(0, 0, face) = 0.;
                    o->Get(0, 0, CNB::CENTER) = 1.;
                    o->GetRhs(0) = ux_in;
                    o->GetInitialGuess(0) = ux_in;
                } else if (ind_g.i() + 1 == extent.i()) {
                    // EAST
                    // Neumann
                    o->Get(0, 0, CNB::WEST) = -1.;
                    o->Get(0, 0, CNB::CENTER) = 1.;
                    o->Get(0, 0, CNB::SOUTH) = 0.;
                    o->Get(0, 0, CNB::NORTH) = 0.;
                    o->GetRhs(0) = 0;
                }
            } else if constexpr (staggered == 1) {
                if (ind_g.i() == 0) {
                    // WEST
                    // Dirichlet
                    SC bc_0 = 0.;
                    o->Get(0, 0, CNB::CENTER) -= alpha * o->Get(0, 0, CNB::WEST);
                    o->Get(0, 0, CNB::EAST) += beta * o->Get(0, 0, CNB::WEST);
                    o->GetRhs(0) -= gamma * bc_0 * o->Get(0, 0, CNB::WEST);
                } else if (ind_g.i() + 1 == extent.i()) {
                    // EAST
                    // Neumann
                    o->Get(0, 0, CNB::CENTER) += o->Get(0, 0, CNB::EAST);
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

SC duct2d(SC Re, SC mu, SC rho, SC H, GO ny) {
    SC L{20 * H};
    GO nx = static_cast<GO>(L / H * ny);
    LO num_ghost = 2;
    int num_tsteps = 500;
    SC uin = Re * mu / (rho * H);
    SC Co = 0.25;
    int freq_write = 100;
    dare::ConstantTimeStep<SC> dt{Co * H / ny / uin};
    SC sim_time = num_tsteps * dt;
    GO i_beg = static_cast<GO>(nx * 2 / 3 + num_ghost);
    GO i_end = static_cast<GO>(nx - ny + num_ghost);
    GO j_mid = ny / 2 + num_ghost;

    SC dp_ana = 12 * mu * uin / (H * H);

    IndexGlobal resolution_global(nx, ny);
    VecSC size_global(L, H);

    dare::ExecutionManager exman;
    dare::FileSystemManager fman(&exman, "Verification_Duct2D");
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

    pm.Initialize(&grid, &dt, BCPressure{pm.GetContinuity()}, BCMom<0, 2>{uin}, BCMom<1, 2>{});

    pm.GetParameterList()->set("momentum: Jacobian", "constant");
    pm.GetParameterList()->set("continuity: Jacobian", "constant");

    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    pm.GetMomentum(0)->GetField()->SetValues(0.);
    pm.GetMomentum(1)->GetField()->SetValues(0.);
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
                     &vec_v);
    };
    printer();
    IndexGlobal ind_beg(i_beg, j_mid);
    IndexGlobal ind_end(i_end, j_mid);
    SC dx = grid_s.GetDistances().x() * (i_end - i_beg);
    bool p_beg_local = grid_s.IsLocalInternal(ind_beg);
    bool p_end_local = grid_s.IsLocalInternal(ind_end);
    dare::SetVerbosity(dare::Verbosity::Low);
    SC dp_dx{std::numeric_limits<SC>::lowest()};
    // pm.SetTimer(dare::Timer{});
    while (dt.GetTime() < sim_time) {
        pm.CopyToOld();
        dt.AdvanceTimeStep();
        pm.SolveFlowField();
        SC p_beg{std::numeric_limits<SC>::lowest()};
        SC p_end{std::numeric_limits<SC>::lowest()};
        if (p_beg_local)
            p_beg = pm.GetContinuity()->GetPressure()->GetDataVector().At(grid_s.MapGlobalToLocal(ind_beg), 0);
        if (p_end_local)
            p_end = pm.GetContinuity()->GetPressure()->GetDataVector().At(grid_s.MapGlobalToLocal(ind_end), 0);
        p_beg = exman.Allmax(p_beg);
        p_end = exman.Allmax(p_end);
        dp_dx = (p_beg - p_end) / dx;
        SC err = std::abs((dp_dx - dp_ana) / dp_ana);
        Print(dare::Verbosity::Low) << "Pressure drop: " << dp_dx << " -> Error: " << err << std::endl;
        if (dt.GetTimeStepCounter() % freq_write == 0) {
            printer();
        }
    }
    return dp_dx;
}

int main(int argc, char* argv[]) {
    bool all_second_order{true};
    dare::SetVerbosity(dare::Verbosity::High);
    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        SC H{0.2};
        std::array ny = {10, 20, 30};
        SC Re = 20;
        SC rho = 1000.;
        SC mu = 1e-3;
        std::vector<SC> dp_dx;
        for (auto n : ny)
            dp_dx.push_back(duct2d(Re, mu, rho, H, n));

        SC uin = Re * mu / (rho * H);
        SC dp_dx_ana = 12 * mu * uin / (H * H);

        std::vector<SC> err;
        for (auto dp : dp_dx)
            err.push_back(std::abs(dp / dp_dx_ana - 1.));

        std::vector<SC> G;
        for (auto& n : ny)
            G.push_back(std::log(static_cast<SC>(n)));
        for (auto& e : err)
            e = std::log(e);

        std::vector<SC> m;
        for (std::size_t i{0}; i < G.size() - 1; i++)
            m.push_back((err[i + 1] - err[i]) / (G[i + 1] - G[i]));
        SC order{10.};
        for (auto m_e : m)
            order = std::min(order, std::abs(m_e));

        all_second_order = order > 1.95;
        Print(dare::Verbosity::Low) << "Found following order of convergence: " << order << std::endl;
    }  // NOLINT

    return all_second_order ? 0 : -1;
}
