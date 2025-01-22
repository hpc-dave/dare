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
    using Grid = dare::Cartesian<2>;
    using GridVector = dare::GridVector<Grid, SC, 1>;
    using Field = dare::Field<Grid, SC, 1>;
    using Writer = dare::VTKWriter<Grid>;
    using IndexGlobal = typename Grid::IndexGlobal;
    using IndexLocal = typename Grid::Index;
    using VecSC = typename Grid::VecSC;
    using CNB = typename Grid::NeighborID;

    dare::ScopeGuard scope_guard(&argc, &argv);
    {
            GO nx{100}, ny{10};
            SC L{1}, H{0.1};
            LO num_ghost = 2;
            int freq_write = 1;
            SC dt = 1e-3;
            bool Dijkhuizen = false;

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
            IndexLocal scalar(0, 0), staggered_x(1, 0), staggered_y(0, 1);
            auto grid_s = grid.GetRepresentation(scalar);
            auto grid_x = grid.GetRepresentation(staggered_x);
            auto grid_y = grid.GetRepresentation(staggered_y);

            Field pressure("pressure", grid_s, 2);
            pressure.SetValues(0.);

            Field dp("dp", grid_s, 1);

            Field mom_x("u", grid_x, 2);
            mom_x.SetValues(0.);

            Field mom_y("v", grid_y, 2);
            mom_y.SetValues(0.);

            Field rho("density", grid_s, 2);
            rho.SetValues(1000);

            Field mu("viscosity", grid_s, 2);
            mu.SetValues(1e-3);

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
                using FVStencil = dare::FaceValueStencil<Grid, SC, 1>;
                using DDT = dare::DDT<Grid>;
                using FluxLimiter = dare::MINMOD;
                using TVD = dare::TVD<Grid, SC, FluxLimiter>;

                auto g_r = mblock->GetRepresentation();
                LO loc_o{mblock->GetLocalOrdinal()};
                IndexLocal ind{mblock->GetIndex()};
                const GridVector& u = mom_x.GetDataVector();
                const GridVector& v = mom_y.GetDataVector();

                Divergence div(*g_r, loc_o);
                DDT ddt(*g_r, loc_o, dt);
                TVD tvd(*g_r, loc_o, dare::Vector<2, const GridVector*>(&u, &v));

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
                    SC p_e = pressure.GetDataVector().At(ind, 0);
                    SC p_w = pressure.GetDataVector().At(ind_nb, 0);
                    grad_p = p_e - p_w;
                    grad_p /= g_r->GetDistances().x();
                    grad_p *= g_r->GetCellVolume();
                }
                SC poro = dare::InterpolateToFace(*g_r, ind, CNB::WEST, epsilon.GetDataVector(), 0);

                mblock->GetRhs(0) -= poro * grad_p;

                // stress tensor
                {   // again a scope to limit variable lifetime
                    FVStencil s_mu_eps = dare::InterpolateToFaceStencil(*g_r, ind, mu.GetDataVector());
                    s_mu_eps *= dare::InterpolateToFaceStencil(*g_r, ind, epsilon.GetDataVector());
                    // implicit components
                    auto dn = g_r->GetDistances();
                    auto A_dn = g_r->GetFaceArea() / dn;
                    FVStencil coef_faces;
                    coef_faces(CNB::WEST, 0) = -2. * s_mu_eps(CNB::WEST, 0) * A_dn.x();
                    coef_faces(CNB::EAST, 0) = -2. * s_mu_eps(CNB::EAST, 0) * A_dn.x();
                    if (Dijkhuizen) {
                        coef_faces(CNB::SOUTH, 0) = -s_mu_eps(CNB::SOUTH, 0) * A_dn.y();
                        coef_faces(CNB::NORTH, 0) = -s_mu_eps(CNB::NORTH, 0) * A_dn.y();
                    }
                    // SC coef_w_im = -2. * s_mu_eps(CNB::WEST, 0) * A_dn.x();
                    // SC coef_e_im = -2. * s_mu_eps(CNB::EAST, 0) * A_dn.x();
                    for (CNB face : g_r->GetFaces()) {
                        SC coef = coef_faces(face, 0);
                        mblock->Get(0, 0, face) += coef;
                        mblock->Get(0, 0, CNB::CENTER) -= coef;
                    }
                    // mblock->Get(0, 0, CNB::WEST) += coef_w_im;
                    // mblock->Get(0, 0, CNB::EAST) += coef_e_im;
                    // mblock->Get(0, 0, CNB::CENTER) -= (coef_w_im + coef_e_im);

                    // explicit components
                    // in Principle we could interpolate to points, but with a Cartesian grid, we can pretty
                    // much optimize that with a bit of manual work
                    IndexLocal ind_up(ind), ind_low(ind);
                    ind_low.j() -= 1;
                    SC du_y_s = u(ind_up, 0) - u(ind_low, 0);  // could be treated implicitly, see Dijkuizen option
                    ind_low.j() += 1;
                    ind_up.j() += 1;
                    SC du_y_n = u(ind_up, 0) - u(ind_low, 0);   // could be treated implicitly

                    ind_low = ind_up = ind;
                    ind_low.i() -= 1;
                    SC dv_x_s = v(ind_up, 0) - v(ind_low, 0);
                    ind_low.j() += 1;
                    ind_up.j() += 1;
                    SC dv_x_n = v(ind_up, 0) - v(ind_low, 0);

                    SC coeff_du_dy_ex = 0.;
                    if (!Dijkhuizen)
                        coeff_du_dy_ex = (s_mu_eps(CNB::NORTH, 0) * du_y_n - s_mu_eps(CNB::SOUTH, 0) * du_y_s) * A_dn.y();

                    SC coef_dv_dx_ex =
                        (s_mu_eps(CNB::NORTH, 0) * dv_x_n - s_mu_eps(CNB::SOUTH, 0) * dv_x_s) * A_dn.x();

                    // later add 3D case here!
                    mblock->GetRhs(0) += coeff_du_dy_ex + coef_dv_dx_ex;
                }
                // Apply BC
                IndexLocal ind_l = g_r->MapOrdinalToIndexLocalInternal(loc_o);
                IndexGlobal ind_g = g_r->MapLocalToGlobal(ind_l);

                if (ind_g.j() == 0) {
                    // SOUTH BOUNDARY
                    // Dirichlet
                    SC bc_0 = 0.;
                    mblock->Get(0, 0, CNB::CENTER) -= mblock->Get(0, 0, CNB::SOUTH);
                    mblock->GetRhs(0) -= 2. * bc_0 * mblock->Get(0, 0, CNB::SOUTH);

                    // Neumann
                    // mblock->Get(0, 0, CNB::CENTER) += mblock->Get(0, 0, CNB::SOUTH);
                    // mblock->Get(0, 0, CNB::SOUTH) = 0.;
                } else if (ind_g.j() == (ny-1)) {
                    // NORTH BOUNDARY
                    // Dirichlet
                    SC bc_0 = 0.;
                    mblock->Get(0, 0, CNB::CENTER) -= mblock->Get(0, 0, CNB::NORTH);
                    mblock->GetRhs(0) -= 2. * bc_0 * mblock->Get(0, 0, CNB::NORTH);

                    // NEUMANN
                    // mblock->Get(0, 0, CNB::CENTER) += mblock->Get(0, 0, CNB::NORTH);
                    // mblock->Get(0, 0, CNB::NORTH) = 0.;
                }
                if (ind_g.i() == 0) {
                    // WEST boundary
                    // Dirichlet
                    SC bc_0 = 1.;
                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->Get(0, 0, CNB::EAST) = 0.;
                    mblock->Get(0, 0, CNB::SOUTH) = 0.;
                    mblock->Get(0, 0, CNB::NORTH) = 0.;
                    mblock->GetRhs(0) = bc_0;
                }
                if (ind_g.i() == nx) {
                    // EAST boundary
                    // Neumann
                    mblock->Get(0, 0, CNB::WEST) = -1.;
                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->Get(0, 0, CNB::SOUTH) = 0.;
                    mblock->Get(0, 0, CNB::NORTH) = 0.;
                    mblock->GetRhs(0) = 0;
                }

                // Removal of coefficients
                if (ind_g.i() == 0)
                    mblock->Remove(0, 0, CNB::WEST);
                if (ind_g.i() == nx)
                    mblock->Remove(0, 0, CNB::EAST);
                if (ind_g.j() == 0)
                    mblock->Remove(0, 0, CNB::SOUTH);
                if (ind_g.j() == (ny - 1))
                    mblock->Remove(0, 0, CNB::NORTH);
            };

            auto build_coef_y = [&](auto mblock) {
                // Normalize!
                using TimeScheme = dare::EULER_BACKWARD;
                using Divergence = dare::Divergence<Grid, TimeScheme>;
                using FVStencil = dare::FaceValueStencil<Grid, SC, 1>;
                using DDT = dare::DDT<Grid>;
                using FluxLimiter = dare::MINMOD;
                using TVD = dare::TVD<Grid, SC, FluxLimiter>;

                auto g_r = mblock->GetRepresentation();
                LO loc_o{mblock->GetLocalOrdinal()};
                IndexLocal ind{mblock->GetIndex()};
                const GridVector& u = mom_x.GetDataVector();
                const GridVector& v = mom_y.GetDataVector();

                Divergence div(*g_r, loc_o);
                DDT ddt(*g_r, loc_o, dt);
                TVD tvd(*g_r, loc_o, dare::Vector<2, const GridVector*>(&u, &v));

                // time derivative
                (*mblock) = ddt(epsilon, rho, mom_y);

                // advection
                (*mblock) += div(epsilon, rho, tvd * v);

                // pressure force
                SC grad_p{0.};
                {
                    // This one here could be more generic by interpolating the pressure drop across all faces
                    // and treat this one here as an optimization for Cartesian grids
                    IndexLocal ind_nb(ind);
                    ind_nb.j() -= 1;
                    grad_p = pressure.GetDataVector(0).At(ind, 0) - pressure.GetDataVector().At(ind_nb, 0);
                    grad_p /= g_r->GetDistances().y();
                    grad_p *= g_r->GetCellVolume();
                }
                SC poro = dare::InterpolateToFace(*g_r, ind, CNB::SOUTH, epsilon.GetDataVector(), 0);

                mblock->GetRhs(0) -= poro * grad_p;

                // stress tensor
                {  // again a scope to limit variable lifetime
                    FVStencil s_mu_eps = dare::InterpolateToFaceStencil(*g_r, ind, mu.GetDataVector());
                    s_mu_eps *= dare::InterpolateToFaceStencil(*g_r, ind, epsilon.GetDataVector());
                    // implicit components
                    auto dn = g_r->GetDistances();
                    auto A_dn = g_r->GetFaceArea() / dn;
                    SC coef_s_im = -2. * s_mu_eps(CNB::SOUTH, 0) * A_dn.y();
                    SC coef_n_im = -2. * s_mu_eps(CNB::NORTH, 0) * A_dn.y();
                    mblock->Get(0, 0, CNB::SOUTH) += coef_s_im;
                    mblock->Get(0, 0, CNB::NORTH) += coef_n_im;
                    mblock->Get(0, 0, CNB::CENTER) -= (coef_s_im + coef_n_im);

                    FVStencil coef_faces;
                    coef_faces(CNB::WEST, 0) = -s_mu_eps(CNB::WEST, 0) * A_dn.x();
                    coef_faces(CNB::EAST, 0) = -s_mu_eps(CNB::EAST, 0) * A_dn.x();
                    if(Dijkhuizen) {
                        coef_faces(CNB::SOUTH, 0) = -2.*s_mu_eps(CNB::SOUTH, 0) * A_dn.y();
                        coef_faces(CNB::NORTH, 0) = -2.*s_mu_eps(CNB::NORTH, 0) * A_dn.y();
                    }
                    for (CNB face : g_r->GetFaces()) {
                        SC coef = coef_faces(face, 0);
                        mblock->Get(0, 0, face) += coef;
                        mblock->Get(0, 0, CNB::CENTER) -= coef;
                    }

                    // explicit components
                    // in Principle we could interpolate to points, but with a Cartesian grid, we can pretty
                    // much optimize that with a bit of manual work
                    IndexLocal ind_up(ind), ind_low(ind);

                    ind_low.i() -= 1;
                    SC dv_x_w = v(ind_up, 0) - v(ind_low, 0);  // could be treated implicitly
                    ind_low.i() += 1;
                    ind_up.i() += 1;
                    SC dv_x_e = v(ind_up, 0) - v(ind_low, 0);  // could be treated implicitly

                    ind_low = ind_up = ind;
                    ind_low.j() -= 1;
                    SC du_y_w = u(ind_up, 0) - u(ind_low, 0);
                    ind_low.i() += 1;
                    ind_up.i() += 1;
                    SC du_y_e = u(ind_up, 0) - u(ind_low, 0);

                    SC coef_dv_dx_ex = 0.;
                    if (!Dijkhuizen)
                        coef_dv_dx_ex = (s_mu_eps(CNB::EAST, 0) * dv_x_e - s_mu_eps(CNB::WEST, 0) * dv_x_w) * A_dn.x();

                    SC coef_du_dy_ex = (s_mu_eps(CNB::EAST, 0) * du_y_e - s_mu_eps(CNB::WEST, 0) * du_y_w) * A_dn.y();

                    // later add 3D case here!
                    mblock->GetRhs(0) += coef_dv_dx_ex + coef_du_dy_ex;
                }

                // Apply BC
                IndexLocal ind_l = g_r->MapOrdinalToIndexLocalInternal(loc_o);
                IndexGlobal ind_g = g_r->MapLocalToGlobal(ind_l);
                if (ind_g.i() == 0) {
                    // WEST boundary
                    // Dirichlet
                    // std::cout << *mblock;
                    SC bc_0 = 0.;
                    mblock->Get(0, 0, CNB::CENTER) -= mblock->Get(0, 0, CNB::WEST);
                    mblock->GetRhs(0) -= 2. * bc_0 * mblock->Get(0, 0, CNB::WEST);
                } else if (ind_g.i() == (nx-1)) {
                    // EAST boundary
                    // Neumann
                    mblock->Get(0, 0, CNB::CENTER) += mblock->Get(0, 0, CNB::EAST);
                }
                if (ind_g.j() == 0) {
                    // SOUTH BOUNDARY
                    // Dirichlet
                    SC bc_0 = 0.;
                    for (CNB face : g_r->GetFaces())
                        mblock->Get(0, 0, face) = 0.;
                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->GetRhs(0) = bc_0;

                    // NEUMANN
                    // mblock->Get(0, 0, CNB::NORTH) = -1.;
                    // mblock->Get(0, 0, CNB::CENTER) = 1.;
                    // mblock->Get(0, 0, CNB::EAST) = 0.;
                    // mblock->Get(0, 0, CNB::WEST) = 0.;
                    // mblock->GetRhs(0) = 0;
                } else if (ind_g.j() == ny) {
                    // NORTH BOUNDARY
                    SC bc_0 = 0.;
                    for (CNB face : g_r->GetFaces())
                        mblock->Get(0, 0, face) = 0.;
                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->GetRhs(0) = bc_0;

                    // mblock->Get(0, 0, CNB::SOUTH) = -1.;
                    // mblock->Get(0, 0, CNB::CENTER) = 1.;
                    // mblock->Get(0, 0, CNB::WEST) = 0.;
                    // mblock->Get(0, 0, CNB::EAST) = 0.;
                    // mblock->GetRhs(0) = 0;
                }

                // Removal of coefficients
                if (ind_g.i() == 0)
                    mblock->Remove(0, 0, CNB::WEST);
                if (ind_g.i() == (nx-1))
                    mblock->Remove(0, 0, CNB::EAST);
                if (ind_g.j() == 0)
                    mblock->Remove(0, 0, CNB::SOUTH);
                if (ind_g.j() == ny)
                    mblock->Remove(0, 0, CNB::NORTH);
            };

            auto compute_defect = [&](const auto& grep, LO loc_o, IndexLocal ind) {
                using Divergence = dare::Divergence<Grid, dare::EULER_BACKWARD>;
                using FluxLimiter = dare::CDS;
                using TVD = dare::TVD<Grid, SC, FluxLimiter>;
                using FVStencil = dare::FaceValueStencil<Grid, SC, 1>;
                using CVStencil = dare::CenterValueStencil<Grid, SC, 1>;
                const GridVector& u = mom_x.GetDataVector();
                const GridVector& v = mom_y.GetDataVector();
                Divergence div(grep, loc_o);
                TVD tvd(grep, loc_o, dare::Vector<2, const GridVector*>(&u, &v));

                CVStencil ONES;
                FVStencil uloc;
                uloc(CNB::WEST, 0) = dare::InterpolateToFace(grep, ind, CNB::WEST, u, 0);
                uloc(CNB::EAST, 0) = dare::InterpolateToFace(grep, ind, CNB::EAST, u, 0);
                uloc(CNB::SOUTH, 0) = dare::InterpolateToFace(grep, ind, CNB::SOUTH, v, 0);
                uloc(CNB::NORTH, 0) = dare::InterpolateToFace(grep, ind, CNB::NORTH, v, 0);
                ONES.SetAll(1.);
                SC defect = div(tvd.Interpolate(ONES, ONES) * uloc)[0];
                return defect;
            };

            auto build_coef_p = [&](auto mblock) {
                // Normalize
                // using Divergence = dare::Divergence<Grid, dare::EULER_BACKWARD>;
                // using FluxLimiter = dare::CDS;
                // using TVD = dare::TVD<Grid, SC, FluxLimiter>;
                using FVStencil = dare::FaceValueStencil<Grid, SC, 1>;
                // using CVStencil = dare::CenterValueStencil<Grid, SC, 1>;

                auto g_r = mblock->GetRepresentation();
                LO loc_o{mblock->GetLocalOrdinal()};
                IndexLocal ind{mblock->GetIndex()};
                // const GridVector& u = mom_x.GetDataVector();

                // Divergence div(*g_r, loc_o);
                // TVD tvd(*g_r, loc_o, dare::Vector<1, const GridVector*>(&u));

                // subsitute with interpolation to face and then gradient
                FVStencil s_rho = dare::InterpolateToFaceStencil(*g_r, ind, rho.GetDataVector());
                FVStencil s_beta = dare::InterpolateToFaceStencil(*g_r, ind, beta.GetDataVector());
                FVStencil s_eps = dare::InterpolateToFaceStencil(*g_r, ind, epsilon.GetDataVector());
                FVStencil A_dx;
                A_dx(CNB::WEST, 0) = A_dx(CNB::EAST, 0) = g_r->GetFaceArea().x() / g_r->GetDistances().x();
                A_dx(CNB::SOUTH, 0) = A_dx(CNB::NORTH, 0) = g_r->GetFaceArea().y() / g_r->GetDistances().y();
                // SC rho_w = dare::InterpolateToFace(*g_r, ind, CNB::WEST, rho.GetDataVector(), 0);
                // SC rho_e = dare::InterpolateToFace(*g_r, ind, CNB::EAST, rho.GetDataVector(), 0);
                // SC beta_w = dare::InterpolateToFace(*g_r, ind, CNB::WEST, beta.GetDataVector(), 0);
                // SC beta_e = dare::InterpolateToFace(*g_r, ind, CNB::EAST, beta.GetDataVector(), 0);
                // SC poro_w = dare::InterpolateToFace(*g_r, ind, CNB::WEST, epsilon.GetDataVector(), 0);
                // SC poro_e = dare::InterpolateToFace(*g_r, ind, CNB::EAST, epsilon.GetDataVector(), 0);

                for (CNB face : g_r->GetFaces()) {
                    // CNB face = dare::ToCartesianNeighbor(face + 1);
                    mblock->Get(0, 0, face) = -s_eps(face, 0) * dt
                            / (s_eps(face, 0) * s_rho(face, 0) + s_beta(face, 0) * dt)
                            * A_dx(face, 0);
                    mblock->Get(0, 0, CNB::CENTER) -= mblock->Get(0, 0, face);
                }

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

                if (ind_g.j() == 0) {
                    mblock->Get(0, 0, CNB::CENTER) += mblock->Get(0, 0, CNB::SOUTH);
                } else if (ind_g.j() == (ny - 1)) {
                    mblock->Get(0, 0, CNB::CENTER) += mblock->Get(0, 0, CNB::NORTH);
                }

                if (ind_g.i() == 0) {
                    mblock->Get(0, 0, CNB::CENTER) += mblock->Get(0, 0, CNB::WEST);
                } else if (ind_g.i() == (nx - 1)) {
                    // EAST boundary
                    // prescribed pressure
                    SC bc_1 = 1.;

                    for (CNB face : g_r->GetFaces())
                        mblock->Get(0, 0, face) = 0.;

                    mblock->Get(0, 0, CNB::CENTER) = 1.;
                    mblock->GetRhs(0) = bc_1 - pressure.GetDataVector().At(ind, 0);
                }

                if (ind_g.j() == 0)
                    mblock->Remove(0, 0, CNB::SOUTH);
                else if (ind_g.j() == (ny - 1))
                    mblock->Remove(0, 0, CNB::NORTH);
                if (ind_g.i() == 0)
                    mblock->Remove(0, 0, CNB::WEST);
                else if (ind_g.i() == (nx - 1))
                    mblock->Remove(0, 0, CNB::EAST);
            };

            auto update_boundaries_x = [&]() {
                // NOT SAFE FOR Parallel execution
                // Update boundary cells
                for (LO j{0}; j < grid_x.GetLocalResolutionInternal().j(); j++) {
                    // WEST
                    mom_x.GetDataVector().At(IndexLocal(num_ghost - 1, j + num_ghost), 0)
                        = mom_x.GetDataVector().At(IndexLocal(num_ghost, j + num_ghost), 0);
                    // EAST at prescribed boundary
                    mom_x.GetDataVector().At(IndexLocal(nx + num_ghost, j + num_ghost), 0)
                        = mom_x.GetDataVector().At(IndexLocal(nx + num_ghost - 1, j + num_ghost), 0);
                    mom_x.GetDataVector().At(IndexLocal(nx + num_ghost + 1, j + num_ghost), 0)
                        = mom_x.GetDataVector().At(IndexLocal(nx + num_ghost, j + num_ghost), 0);
                }
                for (LO i{0}; i < grid_x.GetLocalResolutionInternal().i(); i++) {
                    // SOUTH
                    mom_x.GetDataVector().At(IndexLocal(i + num_ghost, num_ghost - 1), 0)
                        = -mom_x.GetDataVector().At(IndexLocal(i + num_ghost, num_ghost), 0);
                    // NORTH
                    mom_x.GetDataVector().At(IndexLocal(i + num_ghost, ny + num_ghost), 0)
                        = -mom_x.GetDataVector().At(IndexLocal(i + num_ghost, ny + num_ghost - 1), 0);
                }
            };

            auto update_boundaries_y = [&]() {
                // NOT SAFE FOR Parallel execution
                // Update boundary cells
                for (LO j{0}; j < grid_y.GetLocalResolutionInternal().j(); j++) {
                    // WEST
                    mom_y.GetDataVector().At(IndexLocal(num_ghost - 1, j + num_ghost), 0)
                        = -mom_y.GetDataVector().At(IndexLocal(num_ghost, j + num_ghost), 0);
                    // EAST (prescribed)
                    mom_y.GetDataVector().At(IndexLocal(nx + num_ghost, j + num_ghost), 0)
                        = mom_y.GetDataVector().At(IndexLocal(nx + num_ghost - 1, j + num_ghost), 0);
                }
                for (LO i{0}; i < grid_y.GetLocalResolutionInternal().i(); i++) {
                    // SOUTH
                    mom_y.GetDataVector().At(IndexLocal(i + num_ghost, num_ghost - 1), 0)
                        = -mom_y.GetDataVector().At(IndexLocal(i + num_ghost, num_ghost), 0);
                    // NORTH
                    mom_y.GetDataVector().At(IndexLocal(i + num_ghost, ny + num_ghost+1), 0)
                        = -mom_y.GetDataVector().At(IndexLocal(i + num_ghost, ny + num_ghost), 0);
                }
            };

            dare::Trilinos<SC> msystem_x(&exman);
            dare::Trilinos<SC> msystem_y(&exman);
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

            SC t_end = 0.1;
            SC time = 0.;
            int timestep = 0;

            {
                // in a scope, so it won't persist longer than required
                Writer writer(&exman, time, timestep);
                writer.Write(fman, &mom_x.GetDataVector(), &mom_y.GetDataVector(), &pressure.GetDataVector());
            }

            while (time < t_end) {
                timestep++;
                time += dt;

                mom_x.CopyDataVectorsToOldTimeStep();
                mom_y.CopyDataVectorsToOldTimeStep();
                pressure.CopyDataVectorsToOldTimeStep();
                rho.CopyDataVectorsToOldTimeStep();
                epsilon.CopyDataVectorsToOldTimeStep();

                exman.Print(dare::Verbosity::Low) << "t: " << time << "\tstep: " << timestep << "\n";
                msystem_x.Build(grid_x, mom_x.GetDataVector(), build_coef_x, false);
                // msystem_x.PrintMatrix();
                // msystem_x.PrintX();
                // msystem_x.PrintB();
                msystem_y.Build(grid_y, mom_y.GetDataVector(), build_coef_y, false);

                // msystem.PrintMatrix();
                // msystem_y.PrintX();
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

                update_boundaries_x();

                // mom_x.GetDataVector().At(num_ghost - 1, 0) = mom_x.GetDataVector().At(num_ghost, 0);
                // mom_x.GetDataVector().At(nx + num_ghost + 1, 0) = mom_x.GetDataVector().At(nx + num_ghost, 0);
                exman.Print(dare::Verbosity::Low)
                    << "U: " << solver_x.GetNumIterations() << "\t";

                dare::TrilinosSolver<SC> solver_y;
                msystem_y.GetM() = solver_y.BuildPreconditioner(dare::PreCondPackage::Ifpack2,
                                                                "ILUT",
                                                                p_ilu,
                                                                msystem_y.GetA());
                auto ret_y = solver_y.Solve(dare::SolverPackage::Belos,
                                            "BICGSTAB",
                                            msystem_y.GetM(),
                                            msystem_y.GetA(), msystem_y.GetX(), msystem_y.GetB(),
                                            p_solver);
                msystem_y.CopyTo(&mom_y.GetDataVector());
                if (ret_y != Belos::ReturnType::Converged) {
                    exman.Print(dare::Verbosity::Low) << "momentum did not converge\n";
                }

                mom_y.ExchangeHaloCells();

                update_boundaries_y();

                exman.Print(dare::Verbosity::Low)
                    << "V: " << solver_y.GetNumIterations() << "\n";

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

                    for (LO lo_int{0}; lo_int < grid_x.GetNumberLocalCellsInternal(); lo_int++) {
                        LO lo = grid_x.MapInternalToLocal(lo_int);
                        IndexLocal ind = grid_x.MapOrdinalToIndexLocal(lo);

                        // skip at inlet
                        if (ind.i() == num_ghost)
                            continue;

                        // skip at prescribed boundary
                        if (ind.i() == (nx + num_ghost))
                            continue;

                        IndexLocal ind_prev(ind);
                        ind_prev.i() -= 1;

                        SC rho_x = dare::InterpolateToFace(grid_s, ind, CNB::WEST, rho.GetDataVector(), 0);
                        SC beta_x = dare::InterpolateToFace(grid_s, ind, CNB::WEST, beta.GetDataVector(), 0);
                        SC poro_x = dare::InterpolateToFace(grid_s, ind, CNB::WEST, epsilon.GetDataVector(), 0);
                        SC beta_xit = it == 0 ? beta_x : 0.;
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

                    for (LO lo_int{0}; lo_int < grid_y.GetNumberLocalCellsInternal(); lo_int++) {
                        LO lo = grid_y.MapInternalToLocal(lo_int);
                        IndexLocal ind = grid_y.MapOrdinalToIndexLocal(lo);

                        // skip at South
                        if (ind.j() == num_ghost)
                            continue;

                        // skip at North
                        if (ind.j() == (ny + num_ghost))
                            continue;

                        IndexLocal ind_prev(ind);
                        ind_prev.j() -= 1;

                        SC rho_y = dare::InterpolateToFace(grid_s, ind, CNB::SOUTH, rho.GetDataVector(), 0);
                        SC beta_y = dare::InterpolateToFace(grid_s, ind, CNB::SOUTH, beta.GetDataVector(), 0);
                        SC poro_y = dare::InterpolateToFace(grid_s, ind, CNB::SOUTH, epsilon.GetDataVector(), 0);
                        SC beta_yit = it == 0 ? beta_y : 0.;
                        SC dp_n = dp.GetDataVector().At(ind, 0);
                        SC dp_s = dp.GetDataVector().At(ind_prev, 0);
                        SC grad_dp = dp_n - dp_s;
                        grad_dp /= grid_s.GetDistances().y();
                        SC u_n = mom_y.GetDataVector().At(lo);
                        SC denominator_1 = 1. / (poro_y * rho_y + beta_yit * dt);
                        SC denominator_2 = 1. / (poro_y * rho_y + beta_y * dt);
                        SC f1 = denominator_1 * poro_y * rho_y * u_n;
                        SC f2 = -denominator_2 * poro_y * grad_dp * dt;
                        // SC delta_u = -dt / rho_x * grad_dp;
                        SC u_nn = f1 + f2;
                        mom_y.GetDataVector().At(lo) = u_nn;
                    }

                    // check for defect
                    SC defect_max{0.};
                    for (LO lo_int{0}; lo_int < grid_s.GetNumberLocalCellsInternal(); lo_int++) {
                        LO lo = grid_x.MapInternalToLocal(lo_int);
                        IndexLocal ind = grid_x.MapOrdinalToIndexLocal(lo);
                        // skip at prescribed boundary
                        if (ind.i() == (num_ghost + nx - 1))
                            continue;
                        SC defect_loc = std::abs(compute_defect(pressure.GetGridRepresentation(), lo, ind));
                        defect_max = std::max(defect_max, defect_loc);
                    }
                    exman.Print(dare::Verbosity::Low) << "it: " << it << " Defect: " << defect_max << std::endl;
                    if (defect_max < 1e-14)
                        break;
                }
                // Update boundary cells
                update_boundaries_x();
                update_boundaries_y();


                if ((timestep % freq_write) == 0) {
                    Writer writer(&exman, time, timestep);
                    writer.Write(fman,
                                 &mom_x.GetDataVector(),
                                 &mom_y.GetDataVector(),
                                 &pressure.GetDataVector());
                }
            }
    }
    return 0;
}
