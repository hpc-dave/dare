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

#include "Algorithm/NavierStokes/ProjectionMethod.h"
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
#include "Algorithm/NavierStokes/ProjectionMethod.h"
#include "Algorithm/NavierStokes/ProjectionMethod_Cartesian.h"
#include "Algorithm/ConstantTimeStep.h"

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
    template<typename T>
    void Apply(T* o) const {
        if constexpr(dare::FieldType<T>) {
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
                    ind_nb.j() = ind_orig.j()+ j;
                    SC dp_c{o->GetDataVector(0).At(ind, 0)};
                    o->GetDataVector(0).At(ind_nb, 0) = dp_c;
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
                    SC p_1{o->GetDataVector(0).At(ind, 0)};
                    SC delta_p{p_bc - p_1};
                    SC delta_dp{(p_bc + 0.5 * delta_p) - (p_1 + delta_p)};
                    o->GetDataVector(0).At(ind_nb, 0) = delta_dp;
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
        } else if constexpr(o->IsGlobal()) {
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
    SC p_bc{1};
    const CRefType& continuity;  //!< reference to the continuity equation
};

template<int staggered>
struct BCMom {
public:
    BCMom() {}

    template<typename T>
    void Apply(T* mb) const {
        if constexpr(dare::FieldType<T>) {
            // ghost cell values
        } else if constexpr(mb->IsGlobal()) {
            auto g_r = mb->GetRepresentation();
            LO loc_o = mb->GetLocalOrdinal();
            IndexGlobal extent = g_r->GetGlobalResolutionInternal();
            IndexLocal ind_l = g_r->MapOrdinalToIndexLocalInternal(loc_o);
            IndexGlobal ind_g = g_r->MapLocalToGlobal(ind_l);

            if constexpr (staggered == 0) {
                if (ind_g.j() == 0) {
                    // SOUTH - Wall
                        SC bc_wall{0.};
                    mb->Get(0, 0, CNB::CENTER) -= mb->Get(0, 0, CNB::SOUTH);
                    mb->GetRhs(0) -= 2. * bc_wall * mb->Get(0, 0, CNB::SOUTH);
                } else if (ind_g.j() + 1 == extent.j()) {
                    // NORTH
                    // SOUTH - Wall
                    SC bc_wall{0.};
                    mb->Get(0, 0, CNB::CENTER) -= mb->Get(0, 0, CNB::NORTH);
                    mb->GetRhs(0) -= 2. * bc_wall * mb->Get(0, 0, CNB::NORTH);
                }
                if (ind_g.i() == 0) {
                    // WEST
                    // Dirichlet
                    SC bc_0{1.};
                    for (CNB face : g_r->GetFaces())
                        mb->Get(0, 0, face) = 0.;
                    mb->Get(0, 0, CNB::CENTER) = 1.;
                    mb->GetRhs(0) = bc_0;
                } else if (ind_g.i() + 1 == extent.i()) {
                    // EAST
                    // Neumann
                    mb->Get(0, 0, CNB::WEST) = -1.;
                    mb->Get(0, 0, CNB::CENTER) = 1.;
                    mb->Get(0, 0, CNB::SOUTH) = 0.;
                    mb->Get(0, 0, CNB::NORTH) = 0.;
                    mb->GetRhs(0) = 0;
                }
            } else if constexpr (staggered == 1) {
                if (ind_g.i() == 0) {
                    // WEST
                    // Dirichlet
                    SC bc_0 = 0.;
                    mb->Get(0, 0, CNB::CENTER) -= mb->Get(0, 0, CNB::WEST);
                    mb->GetRhs(0) -= 2. * bc_0 * mb->Get(0, 0, CNB::WEST);
                } else if (ind_g.i() + 1 == extent.i()) {
                    // EAST
                    // Neumann
                    mb->Get(0, 0, CNB::CENTER) += mb->Get(0, 0, CNB::EAST);
                }
                if (ind_g.j() == 0) {
                    // SOUTH - Wall
                    SC bc_0 = 0.;
                    for (CNB face : g_r->GetFaces())
                        mb->Get(0, 0, face) = 0.;
                    mb->Get(0, 0, CNB::CENTER) = 1.;
                    mb->GetRhs(0) = bc_0;
                } else if (ind_g.j() + 1 == extent.j()) {
                    // NORTH
                    // SOUTH - Wall
                    SC bc_0 = 0.;
                    for (CNB face : g_r->GetFaces())
                        mb->Get(0, 0, face) = 0.;
                    mb->Get(0, 0, CNB::CENTER) = 1.;
                    mb->GetRhs(0) = bc_0;
                }
            }
            // Removal of coefficients
            if (ind_g.i() == 0)
                mb->Remove(0, 0, CNB::WEST);
            if (ind_g.i() == (extent.i() - 1))
                mb->Remove(0, 0, CNB::EAST);
            if (ind_g.j() == 0)
                mb->Remove(0, 0, CNB::SOUTH);
            if (ind_g.j() == (extent.j() - 1))
                mb->Remove(0, 0, CNB::NORTH);
        }
    }
};

int main(int argc, char* argv[]) {
    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        GO nx{10}, ny{5};
        SC L{1}, H{0.5};
        LO num_ghost = 2;
        // int freq_write = 1;
        dare::ConstantTimeStep<SC> dt{1e-3};

        IndexGlobal resolution_global(nx, ny);
        VecSC size_global(L, H);

        dare::ExecutionManager exman;
        dare::FileSystemManager fman(&exman, "verification");
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

        SC rho = 1.;
        SC mu = 1.;

        pm.Initialize(&grid, &dt, BCPressure{pm.GetContinuity()}, BCMom<0>{}, BCMom<1>{});

        pm.SetDensity(rho);
        pm.SetViscosity(mu);

        pm.SolveFlowField();
    }
    return 0;
}   // NOLINT
