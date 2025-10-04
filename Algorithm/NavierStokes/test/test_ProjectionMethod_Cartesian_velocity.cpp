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

#include "test_ProjectionMethod_Cartesian.h"

TEST_F(ProjectionMethodCartesian1DTest, VelocityUpdate_Incompressible_NoForce) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                    - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[0];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian1DTest, VelocityUpdate_Incompressible) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = Field;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.AddImplicitForce(&beta);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        beta.GetDataVector().At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = 1./(1.+ beta_face * dt/(rho_face * eps_face)) * (u_ex - dt / rho_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    GridVector u_0{pm.GetMomentum(0)->GetField()->GetDataVector()};
    free_pm_update_velocity(&pm, 1);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{u_0.At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[0];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian2DTest, VelocityUpdate_Incompressible_NoForce) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[0];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[1];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian2DTest, VelocityUpdate_Incompressible) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = Field;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.AddImplicitForce(&beta);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        beta.GetDataVector().At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = 1. / (1. + beta_face * dt / (rho_face * eps_face)) * (u_ex - dt / rho_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = 1. / (1. + beta_face * dt / (rho_face * eps_face)) * (u_ex - dt / rho_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    GridVector u_0{pm.GetMomentum(0)->GetField()->GetDataVector()};
    GridVector v_0{pm.GetMomentum(1)->GetField()->GetDataVector()};
    free_pm_update_velocity(&pm, 1);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{u_0.At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[0];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{v_0.At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[1];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian3DTest, VelocityUpdate_Incompressible_NoForce) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[0];

        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[1];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_z = pm.GetMomentum(2)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_z.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[2] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(2)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[2];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(2)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian3DTest, VelocityUpdate_Incompressible) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = Field;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.AddImplicitForce(&beta);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        beta.GetDataVector().At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = 1. / (1. + beta_face * dt / (rho_face * eps_face)) * (u_ex - dt / rho_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = 1. / (1. + beta_face * dt / (rho_face * eps_face)) * (u_ex - dt / rho_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_z = pm.GetMomentum(2)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_z.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[2] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = 1. / (1. + beta_face * dt / (rho_face * eps_face)) * (u_ex - dt / rho_face * dP * dn_r[2]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(2)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    GridVector u_0{pm.GetMomentum(0)->GetField()->GetDataVector()};
    GridVector v_0{pm.GetMomentum(1)->GetField()->GetDataVector()};
    GridVector w_0{pm.GetMomentum(2)->GetField()->GetDataVector()};
    free_pm_update_velocity(&pm, 1);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{u_0.At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[0];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{v_0.At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[1];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    for (LO n_loc = 0; n_loc < g_z.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[2] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{w_0.At(ind, 0)};
        u_ex = u_ex - dt / rho_face * dP * dn_r[2];
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(2)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian1DTest, VelocityUpdate_Compressible_NoForce) {
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return v; }
        SC v;
    };
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
        using density_derivative = density_derivative_functor;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    GridVector rho_prev("rho_prev", g_s);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the compressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);

    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    // Note that we have to set the values of the density field before it is given to the
    // instance of the projection method, since the density at the previous iteration
    // is only copied at this step, since the pressure update step is missing in this artificial
    // situation
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    rho_prev = rho.GetDataVector();
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.SetDensityDerivative(density_derivative_functor{1.});

    // add an offset in the density
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) += 1.;
    }

    free_pm_update_velocity(&pm, 0);
    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian1DTest, VelocityUpdate_Compressible) {
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return v; }
        SC v;
    };
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = Field;
        using density_derivative = density_derivative_functor;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    GridVector rho_prev("rho_prev", g_s);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the compressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);

    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    // Note that we have to set the values of the density field before it is given to the
    // instance of the projection method, since the density at the previous iteration
    // is only copied at this step, since the pressure update step is missing in this artificial
    // situation
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        beta.GetDataVector().At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.AddImplicitForce(&beta);

    rho_prev = rho.GetDataVector();
    // add an offset in the density
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) += 1.;
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / (rho_face + beta_face * dt / eps_face) * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    GridVector u_0{pm.GetMomentum(0)->GetField()->GetDataVector()};
    free_pm_update_velocity(&pm, 1);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{u_0.At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian2DTest, VelocityUpdate_Compressible_NoForce) {
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return v; }
        SC v;
    };
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
        using density_derivative = density_derivative_functor;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    GridVector rho_prev("rho_prev", g_s);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the compressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);

    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    // Note that we have to set the values of the density field before it is given to the
    // instance of the projection method, since the density at the previous iteration
    // is only copied at this step, since the pressure update step is missing in this artificial
    // situation
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    rho_prev = rho.GetDataVector();
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.SetDensityDerivative(density_derivative_functor{1.});

    // add an offset in the density
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) += 1.;
    }

    free_pm_update_velocity(&pm, 0);
    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian2DTest, VelocityUpdate_Compressible) {
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return v; }
        SC v;
    };
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = Field;
        using density_derivative = density_derivative_functor;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    GridVector rho_prev("rho_prev", g_s);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the compressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);

    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    // Note that we have to set the values of the density field before it is given to the
    // instance of the projection method, since the density at the previous iteration
    // is only copied at this step, since the pressure update step is missing in this artificial
    // situation
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        beta.GetDataVector().At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.AddImplicitForce(&beta);

    rho_prev = rho.GetDataVector();
    // add an offset in the density
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) += 1.;
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / (rho_face + beta_face * dt / eps_face) * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / (rho_face + beta_face * dt / eps_face) * (u_ex - dt / rho_prev_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    GridVector u_0{pm.GetMomentum(0)->GetField()->GetDataVector()};
    GridVector v_0{pm.GetMomentum(1)->GetField()->GetDataVector()};
    free_pm_update_velocity(&pm, 1);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{u_0.At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{v_0.At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian3DTest, VelocityUpdate_Compressible_NoForce) {
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return v; }
        SC v;
    };
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
        using density_derivative = density_derivative_functor;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    GridVector rho_prev("rho_prev", g_s);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the compressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);

    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    // Note that we have to set the values of the density field before it is given to the
    // instance of the projection method, since the density at the previous iteration
    // is only copied at this step, since the pressure update step is missing in this artificial
    // situation
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    rho_prev = rho.GetDataVector();
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.SetDensityDerivative(density_derivative_functor{1.});

    // add an offset in the density
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) += 1.;
    }

    free_pm_update_velocity(&pm, 0);
    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_z = pm.GetMomentum(2)->GetField()->GetGridRepresentation();
    for (LO n_loc = 0; n_loc < g_z.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[2] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(2)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[2]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(2)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}

TEST_F(ProjectionMethodCartesian3DTest, VelocityUpdate_Compressible) {
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return v; }
        SC v;
    };
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = Field;
        using density_derivative = density_derivative_functor;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e4;
    Field rho("rho", g_s, 2);
    GridVector rho_prev("rho_prev", g_s);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the compressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);

    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    // Note that we have to set the values of the density field before it is given to the
    // instance of the projection method, since the density at the previous iteration
    // is only copied at this step, since the pressure update step is missing in this artificial
    // situation
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        beta.GetDataVector().At(i) = rd.Generate();
        pm.GetContinuity()->GetdP()->GetDataVector().At(i) = rd.Generate();
    }

    for (std::size_t d{0}; d < Dim; d++) {
        pm.GetMomentum(d)->GetField()->CopyDataVectorsToOldTimeStep();
    }

    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    pm.AddImplicitForce(&beta);

    rho_prev = rho.GetDataVector();
    // add an offset in the density
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) += 1.;
    }

    auto g_x = pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    free_pm_update_velocity(&pm, 0);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / (rho_face + beta_face * dt / eps_face) * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_y = pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / (rho_face + beta_face * dt / eps_face) * (u_ex - dt / rho_prev_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    auto g_z = pm.GetMomentum(2)->GetField()->GetGridRepresentation();
    for (LO n_loc = 0; n_loc < g_z.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[2] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC beta_face = 0.5 * (beta.GetDataVector().At(ind, 0) + beta.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{pm.GetMomentum(2)->GetField()->GetDataVector(1).At(ind, 0)};
        u_ex = rho_prev_face / (rho_face + beta_face * dt / eps_face) * (u_ex - dt / rho_prev_face * dP * dn_r[2]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(2)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    GridVector u_0{pm.GetMomentum(0)->GetField()->GetDataVector()};
    GridVector v_0{pm.GetMomentum(1)->GetField()->GetDataVector()};
    GridVector w_0{pm.GetMomentum(2)->GetField()->GetDataVector()};
    free_pm_update_velocity(&pm, 1);

    for (LO n_loc = 0; n_loc < g_x.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[0] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{u_0.At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[0]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(0)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    for (LO n_loc = 0; n_loc < g_y.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[1] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{v_0.At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[1]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(1)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }

    for (LO n_loc = 0; n_loc < g_z.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z.MapInternalToLocal(ind_loc);

        Index ind_lo{ind};
        ind_lo[2] -= 1;
        SC rho_face = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_lo, 0));
        SC rho_prev_face = 0.5 * (rho_prev.At(ind, 0) + rho_prev.At(ind_lo, 0));
        SC eps_face = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_lo, 0));
        SC dP = pm.GetContinuity()->GetdP()->GetDataVector().At(ind, 0)
                - pm.GetContinuity()->GetdP()->GetDataVector().At(ind_lo, 0);
        SC u_ex{w_0.At(ind, 0)};
        u_ex = rho_prev_face / rho_face * (u_ex - dt / rho_prev_face * dP * dn_r[2]);
        if (std::abs(rho_face) > 1e-15 && std::abs(eps_face) > 1e-15 && std::abs(rho_prev_face) > 1e-15) {
            EXPECT_NEAR(pm.GetMomentum(2)->GetField()->GetDataVector(0).At(ind, 0), u_ex,
                        tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(u_ex));
        }
    }
}
