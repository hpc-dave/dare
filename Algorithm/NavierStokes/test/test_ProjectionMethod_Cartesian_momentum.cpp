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

TEST_F(ProjectionMethodCartesian1DTest, Initialization) {
    dare::ProjectionMethod<GridType, BStrat, PDefault, NDefault> pm;
    dare::test::BStrat bstrat;
    dare::ConstantTimeStep dt(1.);
    EXPECT_FALSE(pm.IsInitialized());
    pm.Initialize(grid, &dt, bstrat, bstrat);
    EXPECT_TRUE(pm.IsInitialized());
    EXPECT_FALSE(pm.CheckStatus());
    for (std::size_t d{0}; d < Dim; d++) {
        EXPECT_TRUE(pm.GetMomentum(d)->GetGridRepresentation()->GetOptions() == opt_m[d]);
    }
}

TEST_F(ProjectionMethodCartesian2DTest, Initialization) {
    dare::ProjectionMethod<GridType, BStrat, PDefault, NDefault> pm;
    dare::test::BStrat bstrat;
    dare::ConstantTimeStep dt(1.);
    EXPECT_FALSE(pm.IsInitialized());
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    EXPECT_TRUE(pm.IsInitialized());
    EXPECT_FALSE(pm.CheckStatus());
    for (std::size_t d{0}; d < Dim; d++) {
        EXPECT_TRUE(pm.GetMomentum(d)->GetGridRepresentation()->GetOptions() == opt_m[d]);
    }
}

TEST_F(ProjectionMethodCartesian3DTest, Initialization) {
    dare::ProjectionMethod<GridType, BStrat, PDefault, NDefault> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    EXPECT_FALSE(pm.IsInitialized());
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    EXPECT_TRUE(pm.IsInitialized());
    EXPECT_FALSE(pm.CheckStatus());
    for (std::size_t d{0}; d < Dim; d++) {
        EXPECT_TRUE(pm.GetMomentum(d)->GetGridRepresentation()->GetOptions() == opt_m[d]);
    }
}

TEST_F(ProjectionMethodCartesian1DTest, FinalizeWithForce) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using explicit_force = Field;
        using implicit_force = Field;
    };
    dare::ConstantTimeStep dt(1.);
    Field rho("rho", grid->GetRepresentation(opt_s), 2);
    double mu = 1.;
    Field beta_im("beta_im", grid->GetRepresentation(opt_s), 1);
    Field beta_ex("beta_ex", grid->GetRepresentation(opt_s), 1);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDefault> pm;
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.AddImplicitForce(&beta_im);
    EXPECT_FALSE(pm.CheckStatus());
    for (std::size_t d{0}; d < Dim; d++) {
        pm.AddExplicitForce(&beta_ex, d);
        EXPECT_FALSE(pm.CheckStatus());
    }
    pm.SetDensity(&rho);
    EXPECT_FALSE(pm.CheckStatus());
    pm.SetViscosity(mu);
    EXPECT_TRUE(pm.CheckStatus());
}

TEST_F(ProjectionMethodCartesian2DTest, FinalizeWithForce) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using explicit_force = Field;
        using implicit_force = Field;
    };
    dare::ConstantTimeStep dt(1.);
    Field rho("rho", grid->GetRepresentation(opt_s), 2);
    double mu = 1.;
    Field beta_im("beta_im", grid->GetRepresentation(opt_s), 1);
    Field beta_ex("beta_ex", grid->GetRepresentation(opt_s), 1);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDefault> pm;
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.AddImplicitForce(&beta_im);
    EXPECT_FALSE(pm.CheckStatus());
    for (std::size_t d{0}; d < Dim; d++) {
        pm.AddExplicitForce(&beta_ex, d);
        EXPECT_FALSE(pm.CheckStatus());
    }
    pm.SetDensity(&rho);
    EXPECT_FALSE(pm.CheckStatus());
    pm.SetViscosity(mu);
    EXPECT_TRUE(pm.CheckStatus());
}

TEST_F(ProjectionMethodCartesian3DTest, FinalizeWithForce) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using explicit_force = Field;
        using implicit_force = Field;
    };
    dare::ConstantTimeStep dt(1.);
    Field rho("rho", grid->GetRepresentation(opt_s), 2);
    double mu = 1.;
    Field beta_im("beta_im", grid->GetRepresentation(opt_s), 1);
    Field beta_ex("beta_ex", grid->GetRepresentation(opt_s), 1);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDefault> pm;
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.AddImplicitForce(&beta_im);
    EXPECT_FALSE(pm.CheckStatus());
    for (std::size_t d{0}; d < Dim; d++) {
        pm.AddExplicitForce(&beta_ex, d);
        EXPECT_FALSE(pm.CheckStatus());
    }
    pm.SetDensity(&rho);
    EXPECT_FALSE(pm.CheckStatus());
    pm.SetViscosity(mu);
    EXPECT_TRUE(pm.CheckStatus());
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_ddt_test) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDefault> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++)
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
    }
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        LO o_loc = g_x->MapIndexToOrdinalLocalInternal(ind_loc);
        auto s = dare::free_pm_ddt_Cartesian(&pm, dare::ZERO, *g_x, o_loc);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        SC dV_dt = g_x->GetCellVolume() / dt;
        SC u = pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0);
        Index ind_nb(ind);
        ind_nb.i() -= 1;
        SC rho_0 = 0.5 * (rho.GetDataVector(0).At(ind, 0) + rho.GetDataVector(0).At(ind_nb, 0));
        SC rho_1 = 0.5 * (rho.GetDataVector(1).At(ind, 0) + rho.GetDataVector(1).At(ind_nb, 0));
        SC eps_0 = 0.5 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_nb, 0));
        SC eps_1 = 0.5 * (epsilon.GetDataVector(1).At(ind, 0) + epsilon.GetDataVector(1).At(ind_nb, 0));
        SC v_0 = eps_0 * rho_0 * dV_dt;
        SC v_1 = eps_1 * rho_1 * dV_dt * u;
        EXPECT_NEAR(s.Center(0), v_0, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_0));
        EXPECT_NEAR(s.GetRhs(0), v_1, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_1));
        for (auto face : g_x->GetFaces())
            EXPECT_EQ(s.GetValue(face, 0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian2DTest, BuildMomentum_ddt_test) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDefault> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++)
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
    }
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        // LO o_loc = g_x->MapIndexToOrdinalLocalInternal(ind_loc);
        auto s = dare::free_pm_ddt_Cartesian(&pm, dare::ZERO, *g_x, n_loc);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        SC dV_dt = g_x->GetCellVolume() / dt;
        SC u = pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0);
        Index ind_nb(ind);
        ind_nb.i() -= 1;
        SC rho_0 = 0.5 * (rho.GetDataVector(0).At(ind, 0) + rho.GetDataVector(0).At(ind_nb, 0));
        SC rho_1 = 0.5 * (rho.GetDataVector(1).At(ind, 0) + rho.GetDataVector(1).At(ind_nb, 0));
        SC eps_0 = 0.5 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_nb, 0));
        SC eps_1 = 0.5 * (epsilon.GetDataVector(1).At(ind, 0) + epsilon.GetDataVector(1).At(ind_nb, 0));
        SC v_0 = eps_0 * rho_0 * dV_dt;
        SC v_1 = eps_1 * rho_1 * dV_dt * u;
        EXPECT_NEAR(s.Center(0), v_0, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_0));
        EXPECT_NEAR(s.GetRhs(0), v_1, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_1));
        for (auto face : g_x->GetFaces())
            EXPECT_EQ(s.GetValue(face, 0), 0.);
    }
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);
        LO o_loc = g_y->MapIndexToOrdinalLocalInternal(ind_loc);
        auto s = dare::free_pm_ddt_Cartesian(&pm, dare::ONE, *g_y, o_loc);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        SC dV_dt = g_y->GetCellVolume() / dt;
        SC u = pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0);
        Index ind_nb(ind);
        ind_nb.j() -= 1;
        SC rho_0 = 0.5 * (rho.GetDataVector(0).At(ind, 0) + rho.GetDataVector(0).At(ind_nb, 0));
        SC rho_1 = 0.5 * (rho.GetDataVector(1).At(ind, 0) + rho.GetDataVector(1).At(ind_nb, 0));
        SC eps_0 = 0.5 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_nb, 0));
        SC eps_1 = 0.5 * (epsilon.GetDataVector(1).At(ind, 0) + epsilon.GetDataVector(1).At(ind_nb, 0));
        SC v_0 = eps_0 * rho_0 * dV_dt;
        SC v_1 = eps_1 * rho_1 * dV_dt * u;
        EXPECT_NEAR(s.Center(0), v_0, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_0));
        EXPECT_NEAR(s.GetRhs(0), v_1, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_1));
        for (auto face : g_y->GetFaces())
            EXPECT_EQ(s.GetValue(face, 0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian3DTest, BuildMomentum_ddt_test) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDefault> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++)
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
    }
    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector().At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
    }
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    auto g_z = &pm.GetMomentum(2)->GetField()->GetGridRepresentation();

    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        LO o_loc = g_x->MapIndexToOrdinalLocalInternal(ind_loc);
        auto s = dare::free_pm_ddt_Cartesian(&pm, dare::ZERO, *g_x, o_loc);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        SC dV_dt = g_x->GetCellVolume() / dt;
        SC u = pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0);
        Index ind_nb(ind);
        ind_nb.i() -= 1;
        SC rho_0 = 0.5 * (rho.GetDataVector(0).At(ind, 0) + rho.GetDataVector(0).At(ind_nb, 0));
        SC rho_1 = 0.5 * (rho.GetDataVector(1).At(ind, 0) + rho.GetDataVector(1).At(ind_nb, 0));
        SC eps_0 = 0.5 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_nb, 0));
        SC eps_1 = 0.5 * (epsilon.GetDataVector(1).At(ind, 0) + epsilon.GetDataVector(1).At(ind_nb, 0));
        SC v_0 = eps_0 * rho_0 * dV_dt;
        SC v_1 = eps_1 * rho_1 * dV_dt * u;
        EXPECT_NEAR(s.Center(0), v_0, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_0));
        EXPECT_NEAR(s.GetRhs(0), v_1, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_1));
        for (auto face : g_x->GetFaces())
            EXPECT_EQ(s.GetValue(face, 0), 0.);
    }
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);
        LO o_loc = g_y->MapIndexToOrdinalLocalInternal(ind_loc);
        auto s = dare::free_pm_ddt_Cartesian(&pm, dare::ONE, *g_y, o_loc);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        SC dV_dt = g_y->GetCellVolume() / dt;
        SC u = pm.GetMomentum(1)->GetField()->GetDataVector(1).At(ind, 0);
        Index ind_nb(ind);
        ind_nb.j() -= 1;
        SC rho_0 = 0.5 * (rho.GetDataVector(0).At(ind, 0) + rho.GetDataVector(0).At(ind_nb, 0));
        SC rho_1 = 0.5 * (rho.GetDataVector(1).At(ind, 0) + rho.GetDataVector(1).At(ind_nb, 0));
        SC eps_0 = 0.5 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_nb, 0));
        SC eps_1 = 0.5 * (epsilon.GetDataVector(1).At(ind, 0) + epsilon.GetDataVector(1).At(ind_nb, 0));
        SC v_0 = eps_0 * rho_0 * dV_dt;
        SC v_1 = eps_1 * rho_1 * dV_dt * u;
        EXPECT_NEAR(s.Center(0), v_0, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_0));
        EXPECT_NEAR(s.GetRhs(0), v_1, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_1));
        for (auto face : g_y->GetFaces())
            EXPECT_EQ(s.GetValue(face, 0), 0.);
    }
    for (LO n_loc = 0; n_loc < g_z->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z->MapInternalToLocal(ind_loc);
        LO o_loc = g_z->MapIndexToOrdinalLocalInternal(ind_loc);
        auto s = dare::free_pm_ddt_Cartesian(&pm, dare::TWO, *g_z, o_loc);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        SC dV_dt = g_z->GetCellVolume() / dt;
        SC u = pm.GetMomentum(2)->GetField()->GetDataVector(1).At(ind, 0);
        Index ind_nb(ind);
        ind_nb.k() -= 1;
        SC rho_0 = 0.5 * (rho.GetDataVector(0).At(ind, 0) + rho.GetDataVector(0).At(ind_nb, 0));
        SC rho_1 = 0.5 * (rho.GetDataVector(1).At(ind, 0) + rho.GetDataVector(1).At(ind_nb, 0));
        SC eps_0 = 0.5 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_nb, 0));
        SC eps_1 = 0.5 * (epsilon.GetDataVector(1).At(ind, 0) + epsilon.GetDataVector(1).At(ind_nb, 0));
        SC v_0 = eps_0 * rho_0 * dV_dt;
        SC v_1 = eps_1 * rho_1 * dV_dt * u;
        EXPECT_NEAR(s.Center(0), v_0, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_0));
        EXPECT_NEAR(s.GetRhs(0), v_1, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_1));
        for (auto face : g_z->GetFaces())
            EXPECT_EQ(s.GetValue(face, 0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_conv_upwind_test) {
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
    using TVD = dare::TVD<GridType, SC, dare::CDS>;
    // using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        TVD tvd(*g_x, n_loc, velocities);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ZERO, *g_x, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_ww(ind), ind_e(ind);
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        ind_e.i() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind_ww, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind_ww, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_w, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_e, 0) + velocities[0]->At(ind, 0));

        SC v_w{0.}, v_c{0.}, v_e{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
        }
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_EQ(s.GetRhs(0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian2DTest, BuildMomentum_conv_upwind_test) {
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
    using TVD = dare::TVD<GridType, SC, dare::CDS>;
    // using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        TVD tvd(*g_x, n_loc, velocities);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ZERO, *g_x, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_ww(ind), ind_s(ind), ind_sw(ind), ind_n(ind), ind_nw(ind), ind_e(ind);
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_nw.i() -= 1;
        ind_nw.j() += 1;
        ind_e.i() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind_ww, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_s, 0) + rho.GetDataVector().At(ind_sw, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_n, 0) + rho.GetDataVector().At(ind_nw, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind_ww, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_s, 0) + epsilon.GetDataVector().At(ind_sw, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_n, 0) + epsilon.GetDataVector().At(ind_nw, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_w, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_e, 0) + velocities[0]->At(ind, 0));
        SC u_s = 0.5 * (velocities[1]->At(ind_w, 0) + velocities[1]->At(ind, 0));
        SC u_n = 0.5 * (velocities[1]->At(ind_nw, 0) + velocities[1]->At(ind_n, 0));

        SC v_w{0.}, v_s{0.}, v_c{0.}, v_e{0.}, v_n{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
        }

        if (u_s < 0) {
            // downwind south
            v_c -= eps_c * rho_c * u_s * dA[1];
        } else {
            // upwind south
            v_s -= eps_s * rho_s * u_s * dA[1];
        }

        if (u_n < 0) {
            // downwind north
            v_n += eps_n * rho_n * u_n * dA[1];
        } else {
            // upwind north
            v_c += eps_c * rho_c * u_n * dA[1];
        }

        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_EQ(s.GetRhs(0), 0.);
    }

    dA = g_y->GetFaceArea();
    // in Y-momentum
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);

        TVD tvd(*g_y, n_loc, velocities);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ONE, *g_y, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_ss(ind), ind_s(ind), ind_sw(ind), ind_n(ind), ind_se(ind), ind_e(ind);
        ind_s.j() -= 1;
        ind_ss.j() -= 2;
        ind_w.i() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_se.i() += 1;
        ind_se.j() -= 1;
        ind_e.i() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_sw, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_ss, 0) + rho.GetDataVector().At(ind_s, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_s, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_se, 0) + rho.GetDataVector().At(ind_e, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_n, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_sw, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_ss, 0) + epsilon.GetDataVector().At(ind_s, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_s, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_se, 0) + epsilon.GetDataVector().At(ind_e, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_n, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_s, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_se, 0) + velocities[0]->At(ind_e, 0));
        SC u_s = 0.5 * (velocities[1]->At(ind_s, 0) + velocities[1]->At(ind, 0));
        SC u_n = 0.5 * (velocities[1]->At(ind_n, 0) + velocities[1]->At(ind, 0));

        SC v_w{0.}, v_s{0.}, v_c{0.}, v_e{0.}, v_n{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
        }

        if (u_s < 0) {
            // downwind south
            v_c -= eps_c * rho_c * u_s * dA[1];
        } else {
            // upwind south
            v_s -= eps_s * rho_s * u_s * dA[1];
        }

        if (u_n < 0) {
            // downwind north
            v_n += eps_n * rho_n * u_n * dA[1];
        } else {
            // upwind north
            v_c += eps_c * rho_c * u_n * dA[1];
        }

        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_EQ(s.GetRhs(0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian3DTest, BuildMomentum_conv_upwind_test) {
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
    using TVD = dare::TVD<GridType, SC, dare::CDS>;
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    auto g_z = &pm.GetMomentum(2)->GetField()->GetGridRepresentation();
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        TVD tvd(*g_x, n_loc, velocities);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ZERO, *g_x, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_ww(ind), ind_s(ind), ind_sw(ind), ind_n(ind), ind_nw(ind), ind_e(ind),
              ind_b(ind), ind_bw(ind), ind_t(ind), ind_tw(ind);
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_nw.i() -= 1;
        ind_nw.j() += 1;
        ind_e.i() += 1;
        ind_b.k() -= 1;
        ind_bw.i() -= 1;
        ind_bw.k() -= 1;
        ind_tw.i() -= 1;
        ind_tw.k() += 1;
        ind_t.k() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind_ww, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_s, 0) + rho.GetDataVector().At(ind_sw, 0));
        SC rho_b = 0.5 * (rho.GetDataVector().At(ind_b, 0) + rho.GetDataVector().At(ind_bw, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_t = 0.5 * (rho.GetDataVector().At(ind_t, 0) + rho.GetDataVector().At(ind_tw, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_n, 0) + rho.GetDataVector().At(ind_nw, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind_ww, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_s, 0) + epsilon.GetDataVector().At(ind_sw, 0));
        SC eps_b = 0.5 * (epsilon.GetDataVector().At(ind_b, 0) + epsilon.GetDataVector().At(ind_bw, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_t = 0.5 * (epsilon.GetDataVector().At(ind_t, 0) + epsilon.GetDataVector().At(ind_tw, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_n, 0) + epsilon.GetDataVector().At(ind_nw, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_w, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_e, 0) + velocities[0]->At(ind, 0));
        SC u_s = 0.5 * (velocities[1]->At(ind_w, 0) + velocities[1]->At(ind, 0));
        SC u_n = 0.5 * (velocities[1]->At(ind_nw, 0) + velocities[1]->At(ind_n, 0));
        SC u_b = 0.5 * (velocities[2]->At(ind_w, 0) + velocities[2]->At(ind, 0));
        SC u_t = 0.5 * (velocities[2]->At(ind_tw, 0) + velocities[2]->At(ind_t, 0));

        SC v_w{0.}, v_s{0.}, v_c{0.}, v_e{0.}, v_n{0.}, v_t{0.}, v_b{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
        }

        if (u_s < 0) {
            // downwind south
            v_c -= eps_c * rho_c * u_s * dA[1];
        } else {
            // upwind south
            v_s -= eps_s * rho_s * u_s * dA[1];
        }

        if (u_n < 0) {
            // downwind north
            v_n += eps_n * rho_n * u_n * dA[1];
        } else {
            // upwind north
            v_c += eps_c * rho_c * u_n * dA[1];
        }

        if (u_b < 0) {
            // downwind bottom
            v_c -= eps_c * rho_c * u_b * dA[2];
        } else {
            // upwind bottom
            v_b -= eps_b * rho_b * u_b * dA[2];
        }

        if (u_t < 0) {
            // downwind top
            v_t += eps_t * rho_t * u_t * dA[2];
        } else {
            // upwind top
            v_c += eps_c * rho_c * u_t * dA[2];
        }

        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_EQ(s.GetRhs(0), 0.);
    }

    dA = g_y->GetFaceArea();
    // in Y-momentum
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ONE, *g_y, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_ss(ind), ind_s(ind), ind_sw(ind), ind_n(ind), ind_se(ind), ind_e(ind),
              ind_b(ind), ind_bs(ind), ind_t(ind), ind_ts(ind);

        ind_s.j() -= 1;
        ind_ss.j() -= 2;
        ind_w.i() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_se.i() += 1;
        ind_se.j() -= 1;
        ind_e.i() += 1;
        ind_b.k() -= 1;
        ind_bs.j() -= 1;
        ind_bs.k() -= 1;
        ind_ts.j() -= 1;
        ind_ts.k() += 1;
        ind_t.k() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_sw, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_ss, 0) + rho.GetDataVector().At(ind_s, 0));
        SC rho_b = 0.5 * (rho.GetDataVector().At(ind_b, 0) + rho.GetDataVector().At(ind_bs, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_s, 0));
        SC rho_t = 0.5 * (rho.GetDataVector().At(ind_t, 0) + rho.GetDataVector().At(ind_ts, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_se, 0) + rho.GetDataVector().At(ind_e, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_n, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_sw, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_ss, 0) + epsilon.GetDataVector().At(ind_s, 0));
        SC eps_b = 0.5 * (epsilon.GetDataVector().At(ind_b, 0) + epsilon.GetDataVector().At(ind_bs, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_s, 0));
        SC eps_t = 0.5 * (epsilon.GetDataVector().At(ind_t, 0) + epsilon.GetDataVector().At(ind_ts, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_se, 0) + epsilon.GetDataVector().At(ind_e, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_n, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_s, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_se, 0) + velocities[0]->At(ind_e, 0));
        SC u_s = 0.5 * (velocities[1]->At(ind_s, 0) + velocities[1]->At(ind, 0));
        SC u_n = 0.5 * (velocities[1]->At(ind_n, 0) + velocities[1]->At(ind, 0));
        SC u_b = 0.5 * (velocities[2]->At(ind_s, 0) + velocities[2]->At(ind, 0));
        SC u_t = 0.5 * (velocities[2]->At(ind_ts, 0) + velocities[2]->At(ind_t, 0));

        SC v_w{0.}, v_s{0.}, v_c{0.}, v_e{0.}, v_n{0.}, v_b{0.}, v_t{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
        }

        if (u_s < 0) {
            // downwind south
            v_c -= eps_c * rho_c * u_s * dA[1];
        } else {
            // upwind south
            v_s -= eps_s * rho_s * u_s * dA[1];
        }

        if (u_n < 0) {
            // downwind north
            v_n += eps_n * rho_n * u_n * dA[1];
        } else {
            // upwind north
            v_c += eps_c * rho_c * u_n * dA[1];
        }

        if (u_b < 0) {
            // downwind bottom
            v_c -= eps_c * rho_c * u_b * dA[2];
        } else {
            // upwind bottom
            v_b -= eps_b * rho_b * u_b * dA[2];
        }

        if (u_t < 0) {
            // downwind top
            v_t += eps_t * rho_t * u_t * dA[2];
        } else {
            // upwind top
            v_c += eps_c * rho_c * u_t * dA[2];
        }

        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_EQ(s.GetRhs(0), 0.);
    }

    dA = g_z->GetFaceArea();
    // in Z-momentum
    for (LO n_loc = 0; n_loc < g_z->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z->MapInternalToLocal(ind_loc);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::TWO, *g_z, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_s(ind), ind_b(ind), ind_bb(ind),
            ind_bw(ind), ind_be(ind), ind_bs(ind), ind_bn(ind),
            ind_t(ind), ind_n(ind), ind_e(ind);
        ind_w.i() -= 1;
        ind_s.j() -= 1;
        ind_b.k() -= 1;
        ind_bb.k() -= 2;
        ind_bw.i() -= 1;
        ind_bw.k() -= 1;
        ind_be.i() += 1;
        ind_be.k() -= 1;
        ind_bs.j() -= 1;
        ind_bs.k() -= 1;
        ind_bn.j() += 1;
        ind_bn.k() -= 1;
        ind_t.k() += 1;
        ind_n.j() += 1;
        ind_e.i() += 1;

        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_bw, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_bs, 0) + rho.GetDataVector().At(ind_s, 0));
        SC rho_b = 0.5 * (rho.GetDataVector().At(ind_bb, 0) + rho.GetDataVector().At(ind_b, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_b, 0));
        SC rho_t = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_t, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_be, 0) + rho.GetDataVector().At(ind_e, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_bn, 0) + rho.GetDataVector().At(ind_n, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_bw, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_bs, 0) + epsilon.GetDataVector().At(ind_s, 0));
        SC eps_b = 0.5 * (epsilon.GetDataVector().At(ind_bb, 0) + epsilon.GetDataVector().At(ind_b, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_b, 0));
        SC eps_t = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_t, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_be, 0) + epsilon.GetDataVector().At(ind_e, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_bn, 0) + epsilon.GetDataVector().At(ind_n, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_b, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_be, 0) + velocities[0]->At(ind_e, 0));
        SC u_s = 0.5 * (velocities[1]->At(ind_b, 0) + velocities[1]->At(ind, 0));
        SC u_n = 0.5 * (velocities[1]->At(ind_bn, 0) + velocities[1]->At(ind_n, 0));
        SC u_b = 0.5 * (velocities[2]->At(ind_b, 0) + velocities[2]->At(ind, 0));
        SC u_t = 0.5 * (velocities[2]->At(ind_t, 0) + velocities[2]->At(ind, 0));

        SC v_w{0.}, v_s{0.}, v_c{0.}, v_e{0.}, v_n{0.}, v_b{0.}, v_t{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
        }

        if (u_s < 0) {
            // downwind south
            v_c -= eps_c * rho_c * u_s * dA[1];
        } else {
            // upwind south
            v_s -= eps_s * rho_s * u_s * dA[1];
        }

        if (u_n < 0) {
            // downwind north
            v_n += eps_n * rho_n * u_n * dA[1];
        } else {
            // upwind north
            v_c += eps_c * rho_c * u_n * dA[1];
        }

        if (u_b < 0) {
            // downwind bottom
            v_c -= eps_c * rho_c * u_b * dA[2];
        } else {
            // upwind bottom
            v_b -= eps_b * rho_b * u_b * dA[2];
        }

        if (u_t< 0) {
            // downwind top
            v_t += eps_t * rho_t * u_t * dA[2];
        } else {
            // upwind top
            v_c += eps_c * rho_c * u_t * dA[2];
        }

        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_EQ(s.GetRhs(0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_conv_cds_test) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ZERO, *g_x, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        // first check implicit upwind
        Index ind_w(ind), ind_ww(ind), ind_e(ind);
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        ind_e.i() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind_ww, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind_ww, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_w, 0) + velocities[0]->At(ind, 0));
        SC u_W = velocities[0]->At(ind_w, 0);
        SC u_e = 0.5 * (velocities[0]->At(ind_e, 0) + velocities[0]->At(ind, 0));
        SC u_E = velocities[0]->At(ind_e, 0);
        SC u_C = velocities[0]->At(ind, 0);

        SC v_w{0.}, v_c{0.}, v_e{0.}, rhs_e{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
            rhs_e -= 0.5* (eps_c * rho_c * u_C - eps_w * rho_w * u_W) * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_w * rho_w * u_W) * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_e * rho_e * u_E) * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_e * rho_e * u_E) * u_e * dA[0];
        }
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_EQ(s.GetRhs(0), rhs_e);
    }
}

TEST_F(ProjectionMethodCartesian2DTest, BuildMomentum_conv_cds_test) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e3;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ZERO, *g_x, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_ww(ind), ind_s(ind), ind_sw(ind), ind_n(ind), ind_nw(ind), ind_e(ind);
        ind_w.i() -= 1;
        ind_ww.i() -= 2;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_nw.i() -= 1;
        ind_nw.j() += 1;
        ind_e.i() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind_ww, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_s, 0) + rho.GetDataVector().At(ind_sw, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_n, 0) + rho.GetDataVector().At(ind_nw, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind_ww, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_s, 0) + epsilon.GetDataVector().At(ind_sw, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_n, 0) + epsilon.GetDataVector().At(ind_nw, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_w, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_e, 0) + velocities[0]->At(ind, 0));
        SC u_s = 0.5 * (velocities[1]->At(ind_w, 0) + velocities[1]->At(ind, 0));
        SC u_n = 0.5 * (velocities[1]->At(ind_nw, 0) + velocities[1]->At(ind_n, 0));
        SC u_C = velocities[0]->At(ind, 0);
        SC u_W = velocities[0]->At(ind_w, 0);
        SC u_E = velocities[0]->At(ind_e, 0);
        SC u_S = velocities[0]->At(ind_s, 0);
        SC u_N = velocities[0]->At(ind_n, 0);

        SC v_w{0.}, v_s{0.}, v_c{0.}, v_e{0.}, v_n{0.}, rhs_e{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_w * rho_w * u_W) * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_w * rho_w * u_W) * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_e * rho_e * u_E) * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_e * rho_e * u_E) * u_e * dA[0];
        }

        if (u_s < 0) {
            // downwind south
            v_c -= eps_c * rho_c * u_s * dA[1];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_s * rho_s * u_S) * u_s * dA[1];
        } else {
            // upwind south
            v_s -= eps_s * rho_s * u_s * dA[1];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_s * rho_s * u_S) * u_s * dA[1];
        }

        if (u_n < 0) {
            // downwind north
            v_n += eps_n * rho_n * u_n * dA[1];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_n * rho_n * u_N) * u_n * dA[1];
        } else {
            // upwind north
            v_c += eps_c * rho_c * u_n * dA[1];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_n * rho_n * u_N) * u_n * dA[1];
        }

        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }

    dA = g_y->GetFaceArea();
    // in Y-momentum
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ONE, *g_y, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_ss(ind), ind_s(ind), ind_sw(ind), ind_n(ind), ind_se(ind), ind_e(ind);
        ind_s.j() -= 1;
        ind_ss.j() -= 2;
        ind_w.i() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_se.i() += 1;
        ind_se.j() -= 1;
        ind_e.i() += 1;
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_sw, 0) + rho.GetDataVector().At(ind_w, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_ss, 0) + rho.GetDataVector().At(ind_s, 0));
        SC rho_c = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_s, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_se, 0) + rho.GetDataVector().At(ind_e, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind, 0) + rho.GetDataVector().At(ind_n, 0));
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_sw, 0) + epsilon.GetDataVector().At(ind_w, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_ss, 0) + epsilon.GetDataVector().At(ind_s, 0));
        SC eps_c = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_s, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_se, 0) + epsilon.GetDataVector().At(ind_e, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_n, 0));
        SC u_w = 0.5 * (velocities[0]->At(ind_s, 0) + velocities[0]->At(ind, 0));
        SC u_e = 0.5 * (velocities[0]->At(ind_se, 0) + velocities[0]->At(ind_e, 0));
        SC u_s = 0.5 * (velocities[1]->At(ind_s, 0) + velocities[1]->At(ind, 0));
        SC u_n = 0.5 * (velocities[1]->At(ind_n, 0) + velocities[1]->At(ind, 0));
        SC u_C = velocities[1]->At(ind, 0);
        SC u_W = velocities[1]->At(ind_w, 0);
        SC u_E = velocities[1]->At(ind_e, 0);
        SC u_S = velocities[1]->At(ind_s, 0);
        SC u_N = velocities[1]->At(ind_n, 0);

        SC v_w{0.}, v_s{0.}, v_c{0.}, v_e{0.}, v_n{0.}, rhs_e{0.};
        if (u_w < 0) {
            // downwind west
            v_c -= eps_c * rho_c * u_w * dA[0];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_w * rho_w * u_W) * u_w * dA[0];
        } else {
            // upwind west
            v_w -= eps_w * rho_w * u_w * dA[0];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_w * rho_w * u_W) * u_w * dA[0];
        }

        if (u_e < 0) {
            // downwind east
            v_e += eps_e * rho_e * u_e * dA[0];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_e * rho_e * u_E) * u_e * dA[0];
        } else {
            // upwind east
            v_c += eps_c * rho_c * u_e * dA[0];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_e * rho_e * u_E) * u_e * dA[0];
        }

        if (u_s < 0) {
            // downwind south
            v_c -= eps_c * rho_c * u_s * dA[1];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_s * rho_s * u_S) * u_s * dA[1];
        } else {
            // upwind south
            v_s -= eps_s * rho_s * u_s * dA[1];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_s * rho_s * u_S) * u_s * dA[1];
        }

        if (u_n < 0) {
            // downwind north
            v_n += eps_n * rho_n * u_n * dA[1];
            rhs_e -= 0.5 * (eps_c * rho_c * u_C - eps_n * rho_n * u_N) * u_n * dA[1];
        } else {
            // upwind north
            v_c += eps_c * rho_c * u_n * dA[1];
            rhs_e += 0.5 * (eps_c * rho_c * u_C - eps_n * rho_n * u_N) * u_n * dA[1];
        }

        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_stress_standard_test) {
    struct PDict {
        using density = double;
        using viscosity = Field;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
        using viscous_stress = dare::PMDefaultStressTensor;
    };
    using CNB = dare::CartesianNeighbor;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    Field mu("mu", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(&mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < mu.GetDataVector().GetSize(); i++) {
        mu.GetDataVector().At(i) = rd.Generate();
        mu.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_x, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_x, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ZERO, *g_x, ind, eps_mu_f, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        EXPECT_EQ(s.GetValue(CNB::CENTER, 0), 0.);
        for (auto face : g_x->GetFaces()) {
            EXPECT_EQ(s.GetValue(face, 0), 0.);
        }
        EXPECT_EQ(s.GetRhs(0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_stress_dijkhuizen_test) {
    struct PDict {
        using density = double;
        using viscosity = Field;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
        using viscous_stress = dare::PMDijkhuizenStressTensor;
    };
    using CNB = dare::CartesianNeighbor;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    Field mu("mu", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(&mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < mu.GetDataVector().GetSize(); i++) {
        mu.GetDataVector().At(i) = rd.Generate();
        mu.GetDataVector(1).At(i) = rd.Generate();
        epsilon.GetDataVector().At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_x, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_x, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ZERO, *g_x, ind, eps_mu_f, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        EXPECT_EQ(s.GetValue(CNB::CENTER, 0), 0.);
        for (auto face : g_x->GetFaces()) {
            EXPECT_EQ(s.GetValue(face, 0), 0.);
        }
        EXPECT_EQ(s.GetRhs(0), 0.);
    }
}

TEST_F(ProjectionMethodCartesian2DTest, BuildMomentum_stress_standard_test) {
    struct PDict {
        using density = double;
        using viscosity = Field;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
        using viscous_stress = dare::PMDefaultStressTensor;
    };
    using CNB = dare::CartesianNeighbor;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    SC tol_eps = 1e3;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    Field mu("mu", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(&mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < mu.GetDataVector().GetSize(); i++) {
        mu.GetDataVector().At(i) = std::abs(rd.Generate());
        mu.GetDataVector(1).At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector().At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector(1).At(i) = std::abs(rd.Generate());
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    dare::Vector<Dim, SC> dn = g_x->GetDistances();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_x, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_x, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ZERO, *g_x, ind, eps_mu_f, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_e(ind), ind_n(ind), ind_nw(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_nw.i() -= 1;
        ind_nw.j() += 1;
        SC mu_w = mu.GetDataVector(0).At(ind_w, 0);
        SC mu_e = mu.GetDataVector(0).At(ind, 0);
        SC mu_s = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_n = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_nw, 0) + mu.GetDataVector(0).At(ind_n, 0));
        SC eps_w = epsilon.GetDataVector(0).At(ind_w, 0);
        SC eps_e = epsilon.GetDataVector(0).At(ind, 0);
        SC eps_s = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_n = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_nw, 0) + epsilon.GetDataVector(0).At(ind_n, 0));
        SC U_n = velocities[0]->At(ind_n, 0);
        SC U_c = velocities[0]->At(ind, 0);
        SC U_s = velocities[0]->At(ind_s, 0);
        SC dv_s = velocities[1]->At(ind, 0) - velocities[1]->At(ind_w, 0);
        SC dv_n = velocities[1]->At(ind_n, 0) - velocities[1]->At(ind_nw, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, rhs_e{0.};
        v_w = -2. * eps_w * mu_w * dA.x() / dn.x();
        v_e = -2. * eps_e * mu_e * dA.x() / dn.x();
        v_c = -v_w - v_e - v_s - v_n;
        rhs_e += mu_n * eps_n * (U_n - U_c) * dA.y() / dn.y();
        rhs_e -= mu_s * eps_s * (U_c - U_s) * dA.y() / dn.y();
        rhs_e += mu_n * eps_n * (dv_n) * dA.y() / dn.x();
        rhs_e -= mu_s * eps_s * (dv_s) * dA.y() / dn.x();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }

    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_y, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_y, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ONE, *g_y, ind, eps_mu_f, velocities);
        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_e(ind), ind_n(ind), ind_se(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_se.i() += 1;
        ind_se.j() -= 1;
        SC mu_w = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_e = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_e, 0)
                        + mu.GetDataVector(0).At(ind_se, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_s = mu.GetDataVector(0).At(ind_s, 0);
        SC mu_n = mu.GetDataVector(0).At(ind, 0);
        SC eps_w = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_e = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_e, 0)
                         + epsilon.GetDataVector(0).At(ind_se, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_s = epsilon.GetDataVector(0).At(ind_s, 0);
        SC eps_n = epsilon.GetDataVector(0).At(ind, 0);
        SC V_w = velocities[1]->At(ind_w, 0);
        SC V_c = velocities[1]->At(ind, 0);
        SC V_e = velocities[1]->At(ind_e, 0);
        SC du_w = velocities[0]->At(ind, 0) - velocities[0]->At(ind_s, 0);
        SC du_e = velocities[0]->At(ind_e, 0) - velocities[0]->At(ind_se, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, rhs_e{0.};
        v_s = -2. * eps_s * mu_s * dA.y() / dn.y();
        v_n = -2. * eps_n * mu_n * dA.y() / dn.y();
        v_c = -v_w - v_e - v_s - v_n;
        rhs_e += mu_e * eps_e * (V_e - V_c) * dA.x() / dn.x();
        rhs_e -= mu_w * eps_w * (V_c - V_w) * dA.x() / dn.x();
        rhs_e += mu_e * eps_e * (du_e)*dA.x() / dn.y();
        rhs_e -= mu_w * eps_w * (du_w)*dA.x() / dn.y();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }
}


TEST_F(ProjectionMethodCartesian2DTest, BuildMomentum_stress_Dijkhuizen_test) {
    struct PDict {
        using density = double;
        using viscosity = Field;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
        using viscous_stress = dare::PMDijkhuizenStressTensor;
    };
    using CNB = dare::CartesianNeighbor;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    SC tol_eps = 1e4;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    Field mu("mu", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(&mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < mu.GetDataVector().GetSize(); i++) {
        mu.GetDataVector().At(i) = std::abs(rd.Generate());
        mu.GetDataVector(1).At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector().At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector(1).At(i) = std::abs(rd.Generate());
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    dare::Vector<Dim, SC> dn = g_x->GetDistances();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_x, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_x, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ZERO, *g_x, ind, eps_mu_f, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_e(ind), ind_n(ind), ind_nw(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_nw.i() -= 1;
        ind_nw.j() += 1;
        SC mu_w = mu.GetDataVector(0).At(ind_w, 0);
        SC mu_e = mu.GetDataVector(0).At(ind, 0);
        SC mu_s = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_n = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_nw, 0) + mu.GetDataVector(0).At(ind_n, 0));
        SC eps_w = epsilon.GetDataVector(0).At(ind_w, 0);
        SC eps_e = epsilon.GetDataVector(0).At(ind, 0);
        SC eps_s = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_n = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_nw, 0) + epsilon.GetDataVector(0).At(ind_n, 0));
        SC dv_s = velocities[1]->At(ind, 0) - velocities[1]->At(ind_w, 0);
        SC dv_n = velocities[1]->At(ind_n, 0) - velocities[1]->At(ind_nw, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, rhs_e{0.};
        v_w = -2. * eps_w * mu_w * dA.x() / dn.x();
        v_e = -2. * eps_e * mu_e * dA.x() / dn.x();
        v_s = -eps_s * mu_s * dA.y() / dn.y();
        v_n = -eps_n * mu_n * dA.y() / dn.y();
        v_c = -v_w - v_e - v_s - v_n;
        rhs_e += mu_n * eps_n * (dv_n) * dA.y() / dn.x();
        rhs_e -= mu_s * eps_s * (dv_s) * dA.y() / dn.x();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }

    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_y, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_y, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ONE, *g_y, ind, eps_mu_f, velocities);
        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_e(ind), ind_n(ind), ind_se(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_n.j() += 1;
        ind_se.i() += 1;
        ind_se.j() -= 1;
        SC mu_w = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_e = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_e, 0)
                        + mu.GetDataVector(0).At(ind_se, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_s = mu.GetDataVector(0).At(ind_s, 0);
        SC mu_n = mu.GetDataVector(0).At(ind, 0);
        SC eps_w = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_e = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_e, 0)
                         + epsilon.GetDataVector(0).At(ind_se, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_s = epsilon.GetDataVector(0).At(ind_s, 0);
        SC eps_n = epsilon.GetDataVector(0).At(ind, 0);
        SC du_w = velocities[0]->At(ind, 0) - velocities[0]->At(ind_s, 0);
        SC du_e = velocities[0]->At(ind_e, 0) - velocities[0]->At(ind_se, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, rhs_e{0.};
        v_w = -eps_w * mu_w * dA.x() / dn.x();
        v_e = -eps_e * mu_e * dA.x() / dn.x();
        v_s = -2. * eps_s * mu_s * dA.y() / dn.y();
        v_n = -2. * eps_n * mu_n * dA.y() / dn.y();
        v_c = -v_w - v_e - v_s - v_n;
        rhs_e += mu_e * eps_e * (du_e)*dA.x() / dn.y();
        rhs_e -= mu_w * eps_w * (du_w)*dA.x() / dn.y();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, BuildMomentum_stress_standard_test) {
    struct PDict {
        using density = double;
        using viscosity = Field;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
        using viscous_stress = dare::PMDefaultStressTensor;
    };
    using CNB = dare::CartesianNeighbor;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    SC tol_eps = 1e4;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    Field mu("mu", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    auto g_z = &pm.GetMomentum(2)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(&mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < mu.GetDataVector().GetSize(); i++) {
        mu.GetDataVector().At(i) = std::abs(rd.Generate());
        mu.GetDataVector(1).At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector().At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector(1).At(i) = std::abs(rd.Generate());
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    dare::Vector<Dim, SC> dn = g_x->GetDistances();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_x, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_x, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ZERO, *g_x, ind, eps_mu_f, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_b(ind), ind_bw(ind), ind_e(ind),
              ind_n(ind), ind_nw(ind), ind_t(ind), ind_tw(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_b.k() -= 1;
        ind_bw.i() -= 1;
        ind_bw.k() -= 1;
        ind_n.j() += 1;
        ind_nw.i() -= 1;
        ind_nw.j() += 1;
        ind_t.k() += 1;
        ind_tw.i() -= 1;
        ind_tw.k() += 1;
        SC mu_w = mu.GetDataVector(0).At(ind_w, 0);
        SC mu_e = mu.GetDataVector(0).At(ind, 0);
        SC mu_s = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_n = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_nw, 0) + mu.GetDataVector(0).At(ind_n, 0));
        SC mu_b = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_bw, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_t = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_tw, 0) + mu.GetDataVector(0).At(ind_t, 0));
        SC eps_w = epsilon.GetDataVector(0).At(ind_w, 0);
        SC eps_e = epsilon.GetDataVector(0).At(ind, 0);
        SC eps_s = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_n = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_nw, 0) + epsilon.GetDataVector(0).At(ind_n, 0));
        SC eps_b = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_bw, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_t = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_tw, 0) + epsilon.GetDataVector(0).At(ind_t, 0));
        SC U_n = velocities[0]->At(ind_n, 0);
        SC U_t = velocities[0]->At(ind_t, 0);
        SC U_c = velocities[0]->At(ind, 0);
        SC U_b = velocities[0]->At(ind_b, 0);
        SC U_s = velocities[0]->At(ind_s, 0);
        SC dv_s = velocities[1]->At(ind, 0) - velocities[1]->At(ind_w, 0);
        SC dv_n = velocities[1]->At(ind_n, 0) - velocities[1]->At(ind_nw, 0);
        SC dw_b = velocities[2]->At(ind, 0) - velocities[2]->At(ind_w, 0);
        SC dw_t = velocities[2]->At(ind_t, 0) - velocities[2]->At(ind_tw, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, v_b{0.}, v_t{0.}, rhs_e{0.};
        v_w = -2. * eps_w * mu_w * dA.x() / dn.x();
        v_e = -2. * eps_e * mu_e * dA.x() / dn.x();
        v_c = -v_w - v_e - v_s - v_n - v_b - v_t;
        rhs_e += mu_n * eps_n * (U_n - U_c) * dA.y() / dn.y();
        rhs_e -= mu_s * eps_s * (U_c - U_s) * dA.y() / dn.y();
        rhs_e += mu_t * eps_t * (U_t - U_c) * dA.z() / dn.z();
        rhs_e -= mu_b * eps_b * (U_c - U_b) * dA.z() / dn.z();
        rhs_e += mu_n * eps_n * (dv_n) * dA.y() / dn.x();
        rhs_e -= mu_s * eps_s * (dv_s) * dA.y() / dn.x();
        rhs_e += mu_t * eps_t * (dw_t) * dA.z() / dn.x();
        rhs_e -= mu_b * eps_b * (dw_b) * dA.z() / dn.x();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }

    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_y, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_y, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ONE, *g_y, ind, eps_mu_f, velocities);
        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_sb(ind), ind_e(ind), ind_n(ind), ind_se(ind),
              ind_st(ind), ind_b(ind), ind_t(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_sb.j() -= 1;
        ind_sb.k() -= 1;
        ind_n.j() += 1;
        ind_se.i() += 1;
        ind_se.j() -= 1;
        ind_t.k() += 1;
        ind_b.k() -= 1;
        ind_st.j() -= 1;
        ind_st.k() += 1;
        SC mu_w = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_e = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_e, 0)
                        + mu.GetDataVector(0).At(ind_se, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_s = mu.GetDataVector(0).At(ind_s, 0);
        SC mu_n = mu.GetDataVector(0).At(ind, 0);
        SC mu_b = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_b, 0)
                        + mu.GetDataVector(0).At(ind_sb, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_t = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_t, 0)
                        + mu.GetDataVector(0).At(ind_st, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC eps_w = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_e = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_e, 0)
                         + epsilon.GetDataVector(0).At(ind_se, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_s = epsilon.GetDataVector(0).At(ind_s, 0);
        SC eps_n = epsilon.GetDataVector(0).At(ind, 0);
        SC eps_b = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_b, 0)
                         + epsilon.GetDataVector(0).At(ind_sb, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_t = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_t, 0)
                         + epsilon.GetDataVector(0).At(ind_st, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC V_w = velocities[1]->At(ind_w, 0);
        SC V_b = velocities[1]->At(ind_b, 0);
        SC V_c = velocities[1]->At(ind, 0);
        SC V_e = velocities[1]->At(ind_e, 0);
        SC V_t = velocities[1]->At(ind_t, 0);
        SC du_w = velocities[0]->At(ind, 0) - velocities[0]->At(ind_s, 0);
        SC du_e = velocities[0]->At(ind_e, 0) - velocities[0]->At(ind_se, 0);
        SC dw_b = velocities[2]->At(ind, 0) - velocities[2]->At(ind_s, 0);
        SC dw_t = velocities[2]->At(ind_t, 0) - velocities[2]->At(ind_st, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, v_b{0.}, v_t{0.}, rhs_e{0.};
        v_s = -2. * eps_s * mu_s * dA.y() / dn.y();
        v_n = -2. * eps_n * mu_n * dA.y() / dn.y();
        v_c = -v_w - v_e - v_s - v_n - v_b - v_t;
        rhs_e += mu_e * eps_e * (V_e - V_c) * dA.x() / dn.x();
        rhs_e -= mu_w * eps_w * (V_c - V_w) * dA.x() / dn.x();
        rhs_e += mu_t * eps_t * (V_t - V_c) * dA.z() / dn.z();
        rhs_e -= mu_b * eps_b * (V_c - V_b) * dA.z() / dn.z();
        rhs_e += mu_e * eps_e * (du_e)*dA.x() / dn.y();
        rhs_e -= mu_w * eps_w * (du_w)*dA.x() / dn.y();
        rhs_e += mu_t * eps_t * (dw_t)*dA.z() / dn.y();
        rhs_e -= mu_b * eps_b * (dw_b)*dA.z() / dn.y();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }

    for (LO n_loc = 0; n_loc < g_z->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_z, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_z, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::TWO, *g_z, ind, eps_mu_f, velocities);
        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_w(ind), ind_s(ind), ind_bw(ind), ind_bs(ind), ind_e(ind), ind_n(ind), ind_be(ind),
              ind_bn(ind), ind_b(ind), ind_t(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_bw.i() -= 1;
        ind_bw.k() -= 1;
        ind_bs.j() -= 1;
        ind_bs.k() -= 1;
        ind_n.j() += 1;
        ind_be.i() += 1;
        ind_be.k() -= 1;
        ind_t.k() += 1;
        ind_b.k() -= 1;
        ind_bn.j() += 1;
        ind_bn.k() -= 1;
        SC mu_w = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_bw, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_e = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_e, 0)
                        + mu.GetDataVector(0).At(ind_be, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_s = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_s, 0)
                        + mu.GetDataVector(0).At(ind_bs, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_n = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_n, 0)
                        + mu.GetDataVector(0).At(ind_bn, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_b = mu.GetDataVector(0).At(ind_b, 0);
        SC mu_t = mu.GetDataVector(0).At(ind, 0);
        SC eps_w = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_bw, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_e = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_e, 0)
                         + epsilon.GetDataVector(0).At(ind_be, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_s = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_s, 0)
                         + epsilon.GetDataVector(0).At(ind_bs, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_n = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_n, 0)
                         + epsilon.GetDataVector(0).At(ind_bn, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_b = epsilon.GetDataVector(0).At(ind_b, 0);
        SC eps_t = epsilon.GetDataVector(0).At(ind, 0);
        SC V_w = velocities[2]->At(ind_w, 0);
        SC V_s = velocities[2]->At(ind_s, 0);
        SC V_c = velocities[2]->At(ind, 0);
        SC V_e = velocities[2]->At(ind_e, 0);
        SC V_n = velocities[2]->At(ind_n, 0);
        SC du_w = velocities[0]->At(ind, 0) - velocities[0]->At(ind_b, 0);
        SC du_e = velocities[0]->At(ind_e, 0) - velocities[0]->At(ind_be, 0);
        SC dv_s = velocities[1]->At(ind, 0) - velocities[1]->At(ind_b, 0);
        SC dv_n = velocities[1]->At(ind_n, 0) - velocities[1]->At(ind_bn, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, v_b{0.}, v_t{0.}, rhs_e{0.};
        v_b = -2. * eps_b * mu_b * dA.z() / dn.z();
        v_t = -2. * eps_t * mu_t * dA.z() / dn.z();
        v_c = -v_w - v_e - v_s - v_n - v_b - v_t;
        rhs_e += mu_e * eps_e * (V_e - V_c) * dA.x() / dn.x();
        rhs_e -= mu_w * eps_w * (V_c - V_w) * dA.x() / dn.x();
        rhs_e += mu_n * eps_n * (V_n - V_c) * dA.y() / dn.y();
        rhs_e -= mu_s * eps_s * (V_c - V_s) * dA.y() / dn.y();
        rhs_e += mu_e * eps_e * (du_e)*dA.x() / dn.z();
        rhs_e -= mu_w * eps_w * (du_w)*dA.x() / dn.z();
        rhs_e += mu_n * eps_n * (dv_n)*dA.y() / dn.z();
        rhs_e -= mu_s * eps_s * (dv_s)*dA.y() / dn.z();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, BuildMomentum_stress_Dijkhuizen_test) {
    struct PDict {
        using density = double;
        using viscosity = Field;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::CDS;
        using time_scheme_convective = dare::EULER_BACKWARD;
        using viscous_stress = dare::PMDijkhuizenStressTensor;
    };
    using CNB = dare::CartesianNeighbor;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    SC tol_eps = 1e4;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    Field mu("mu", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    auto g_z = &pm.GetMomentum(2)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(&mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < mu.GetDataVector().GetSize(); i++) {
        mu.GetDataVector().At(i) = std::abs(rd.Generate());
        mu.GetDataVector(1).At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector().At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector(1).At(i) = std::abs(rd.Generate());
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    dare::Vector<Dim, SC> dn = g_x->GetDistances();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_x, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_x, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ZERO, *g_x, ind, eps_mu_f, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_b(ind), ind_bw(ind), ind_e(ind),
            ind_n(ind), ind_nw(ind), ind_t(ind), ind_tw(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_b.k() -= 1;
        ind_bw.i() -= 1;
        ind_bw.k() -= 1;
        ind_n.j() += 1;
        ind_nw.i() -= 1;
        ind_nw.j() += 1;
        ind_t.k() += 1;
        ind_tw.i() -= 1;
        ind_tw.k() += 1;
        SC mu_w = mu.GetDataVector(0).At(ind_w, 0);
        SC mu_e = mu.GetDataVector(0).At(ind, 0);
        SC mu_s = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_n = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_nw, 0) + mu.GetDataVector(0).At(ind_n, 0));
        SC mu_b = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_bw, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_t = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_tw, 0) + mu.GetDataVector(0).At(ind_t, 0));
        SC eps_w = epsilon.GetDataVector(0).At(ind_w, 0);
        SC eps_e = epsilon.GetDataVector(0).At(ind, 0);
        SC eps_s = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_n = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_nw, 0) + epsilon.GetDataVector(0).At(ind_n, 0));
        SC eps_b = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_bw, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_t = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_tw, 0) + epsilon.GetDataVector(0).At(ind_t, 0));
        SC dv_s = velocities[1]->At(ind, 0) - velocities[1]->At(ind_w, 0);
        SC dv_n = velocities[1]->At(ind_n, 0) - velocities[1]->At(ind_nw, 0);
        SC dw_b = velocities[2]->At(ind, 0) - velocities[2]->At(ind_w, 0);
        SC dw_t = velocities[2]->At(ind_t, 0) - velocities[2]->At(ind_tw, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, v_b{0.}, v_t{0.}, rhs_e{0.};
        v_w = -2. * eps_w * mu_w * dA.x() / dn.x();
        v_e = -2. * eps_e * mu_e * dA.x() / dn.x();
        v_s = -eps_s * mu_s * dA.y() / dn.y();
        v_n = -eps_n * mu_n * dA.y() / dn.y();
        v_b = -eps_b * mu_b * dA.z() / dn.z();
        v_t = -eps_t * mu_t * dA.z() / dn.z();
        v_c = -v_w - v_e - v_s - v_n - v_b - v_t;
        rhs_e += mu_n * eps_n * (dv_n)*dA.y() / dn.x();
        rhs_e -= mu_s * eps_s * (dv_s)*dA.y() / dn.x();
        rhs_e += mu_t * eps_t * (dw_t)*dA.z() / dn.x();
        rhs_e -= mu_b * eps_b * (dw_b)*dA.z() / dn.x();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }

    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_y, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_y, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::ONE, *g_y, ind, eps_mu_f, velocities);
        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_w(ind), ind_s(ind), ind_sw(ind), ind_sb(ind), ind_e(ind), ind_n(ind), ind_se(ind),
            ind_st(ind), ind_b(ind), ind_t(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_sw.i() -= 1;
        ind_sw.j() -= 1;
        ind_sb.j() -= 1;
        ind_sb.k() -= 1;
        ind_n.j() += 1;
        ind_se.i() += 1;
        ind_se.j() -= 1;
        ind_t.k() += 1;
        ind_b.k() -= 1;
        ind_st.j() -= 1;
        ind_st.k() += 1;
        SC mu_w = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_sw, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_e = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_e, 0)
                        + mu.GetDataVector(0).At(ind_se, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_s = mu.GetDataVector(0).At(ind_s, 0);
        SC mu_n = mu.GetDataVector(0).At(ind, 0);
        SC mu_b = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_b, 0)
                        + mu.GetDataVector(0).At(ind_sb, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC mu_t = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_t, 0)
                        + mu.GetDataVector(0).At(ind_st, 0) + mu.GetDataVector(0).At(ind_s, 0));
        SC eps_w = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_sw, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_e = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_e, 0)
                         + epsilon.GetDataVector(0).At(ind_se, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_s = epsilon.GetDataVector(0).At(ind_s, 0);
        SC eps_n = epsilon.GetDataVector(0).At(ind, 0);
        SC eps_b = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_b, 0)
                         + epsilon.GetDataVector(0).At(ind_sb, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC eps_t = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_t, 0)
                         + epsilon.GetDataVector(0).At(ind_st, 0) + epsilon.GetDataVector(0).At(ind_s, 0));
        SC du_w = velocities[0]->At(ind, 0) - velocities[0]->At(ind_s, 0);
        SC du_e = velocities[0]->At(ind_e, 0) - velocities[0]->At(ind_se, 0);
        SC dw_b = velocities[2]->At(ind, 0) - velocities[2]->At(ind_s, 0);
        SC dw_t = velocities[2]->At(ind_t, 0) - velocities[2]->At(ind_st, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, v_b{0.}, v_t{0.}, rhs_e{0.};
        v_w = -eps_w * mu_w * dA.x() / dn.x();
        v_e = -eps_e * mu_e * dA.x() / dn.x();
        v_s = -2. * eps_s * mu_s * dA.y() / dn.y();
        v_n = -2. * eps_n * mu_n * dA.y() / dn.y();
        v_b = -eps_b * mu_b * dA.z() / dn.z();
        v_t = -eps_t * mu_t * dA.z() / dn.z();
        v_c = -v_w - v_e - v_s - v_n - v_b - v_t;
        rhs_e += mu_e * eps_e * (du_e)*dA.x() / dn.y();
        rhs_e -= mu_w * eps_w * (du_w)*dA.x() / dn.y();
        rhs_e += mu_t * eps_t * (dw_t)*dA.z() / dn.y();
        rhs_e -= mu_b * eps_b * (dw_b)*dA.z() / dn.y();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }

    for (LO n_loc = 0; n_loc < g_z->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z->MapInternalToLocal(ind_loc);

        FVStencil mu_f = dare::InterpolateToFaceStencil(*g_z, ind, mu.GetDataVector(0));
        FVStencil epsilon_f = dare::InterpolateToFaceStencil(*g_z, ind, epsilon.GetDataVector(0));
        FVStencil eps_mu_f = mu_f * epsilon_f;

        auto s = dare::free_pm_viscious_stress_Cartesian(&pm, dare::TWO, *g_z, ind, eps_mu_f, velocities);
        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_w(ind), ind_s(ind), ind_bw(ind), ind_bs(ind), ind_e(ind), ind_n(ind), ind_be(ind),
            ind_bn(ind), ind_b(ind), ind_t(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_bw.i() -= 1;
        ind_bw.k() -= 1;
        ind_bs.j() -= 1;
        ind_bs.k() -= 1;
        ind_n.j() += 1;
        ind_be.i() += 1;
        ind_be.k() -= 1;
        ind_t.k() += 1;
        ind_b.k() -= 1;
        ind_bn.j() += 1;
        ind_bn.k() -= 1;
        SC mu_w = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_w, 0)
                        + mu.GetDataVector(0).At(ind_bw, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_e = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_e, 0)
                        + mu.GetDataVector(0).At(ind_be, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_s = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_s, 0)
                        + mu.GetDataVector(0).At(ind_bs, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_n = 0.25 * (mu.GetDataVector(0).At(ind, 0) + mu.GetDataVector(0).At(ind_n, 0)
                        + mu.GetDataVector(0).At(ind_bn, 0) + mu.GetDataVector(0).At(ind_b, 0));
        SC mu_b = mu.GetDataVector(0).At(ind_b, 0);
        SC mu_t = mu.GetDataVector(0).At(ind, 0);
        SC eps_w = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_w, 0)
                         + epsilon.GetDataVector(0).At(ind_bw, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_e = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_e, 0)
                         + epsilon.GetDataVector(0).At(ind_be, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_s = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_s, 0)
                         + epsilon.GetDataVector(0).At(ind_bs, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_n = 0.25 * (epsilon.GetDataVector(0).At(ind, 0) + epsilon.GetDataVector(0).At(ind_n, 0)
                         + epsilon.GetDataVector(0).At(ind_bn, 0) + epsilon.GetDataVector(0).At(ind_b, 0));
        SC eps_b = epsilon.GetDataVector(0).At(ind_b, 0);
        SC eps_t = epsilon.GetDataVector(0).At(ind, 0);
        SC du_w = velocities[0]->At(ind, 0) - velocities[0]->At(ind_b, 0);
        SC du_e = velocities[0]->At(ind_e, 0) - velocities[0]->At(ind_be, 0);
        SC dv_s = velocities[1]->At(ind, 0) - velocities[1]->At(ind_b, 0);
        SC dv_n = velocities[1]->At(ind_n, 0) - velocities[1]->At(ind_bn, 0);

        SC v_c{0.}, v_w{0.}, v_e{0.}, v_s{0.}, v_n{0.}, v_b{0.}, v_t{0.}, rhs_e{0.};
        v_w = -eps_w * mu_w * dA.x() / dn.x();
        v_e = -eps_e * mu_e * dA.x() / dn.x();
        v_s = -eps_s * mu_s * dA.y() / dn.y();
        v_n = -eps_n * mu_n * dA.y() / dn.y();
        v_b = -2. * eps_b * mu_b * dA.z() / dn.z();
        v_t = -2. * eps_t * mu_t * dA.z() / dn.z();
        v_c = -v_w - v_e - v_s - v_n - v_b - v_t;
        rhs_e += mu_e * eps_e * (du_e)*dA.x() / dn.z();
        rhs_e -= mu_w * eps_w * (du_w)*dA.x() / dn.z();
        rhs_e += mu_n * eps_n * (dv_n)*dA.y() / dn.z();
        rhs_e -= mu_s * eps_s * (dv_s)*dA.y() / dn.z();
        EXPECT_NEAR(s.GetValue(CNB::WEST, 0), v_w, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_w));
        EXPECT_NEAR(s.GetValue(CNB::EAST, 0), v_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_e));
        EXPECT_NEAR(s.GetValue(CNB::SOUTH, 0), v_s, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_s));
        EXPECT_NEAR(s.GetValue(CNB::NORTH, 0), v_n, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_n));
        EXPECT_NEAR(s.GetValue(CNB::BOTTOM, 0), v_b, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_b));
        EXPECT_NEAR(s.GetValue(CNB::TOP, 0), v_t, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_t));
        EXPECT_NEAR(s.GetValue(CNB::CENTER, 0), v_c, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_c));
        EXPECT_NEAR(s.GetRhs(0), rhs_e, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(rhs_e));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_pressure_test) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    using NDict = dare::PMNumericalInfoDefault;


    SC tol_eps = 1e2;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    double mu = 0;
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < epsilon.GetDataVector().GetSize(); i++) {
        epsilon.GetDataVector().At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector(1).At(i) = std::abs(rd.Generate());
        pm.GetPressure()->GetDataVector().At(i) = rd.Generate();
        pm.GetPressure()->GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    SC dV = g_x->GetCellVolume();
    dare::Vector<Dim, SC> dn = g_x->GetDistances();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_pressure_force_Cartesian(&pm, dare::ZERO, ind, &epsilon);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_low(ind);
        ind_low[0] -= 1;
        SC dp = pm.GetPressure()->GetDataVector().At(ind, 0) - pm.GetPressure()->GetDataVector().At(ind_low, 0);
        SC grad_p = dp / dn[0];
        SC eps_f = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                        + epsilon.GetDataVector().At(ind_low, 0));
        SC e_grad_p = -eps_f * grad_p * dV;
        EXPECT_NEAR(s.GetRhs(0), e_grad_p, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(e_grad_p));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, BuildMomentum_pressure_test) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    using NDict = dare::PMNumericalInfoDefault;

    SC tol_eps = 1e2;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    double mu = 0;
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < epsilon.GetDataVector().GetSize(); i++) {
        epsilon.GetDataVector().At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector(1).At(i) = std::abs(rd.Generate());
        pm.GetPressure()->GetDataVector().At(i) = rd.Generate();
        pm.GetPressure()->GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    SC dV = g_x->GetCellVolume();
    dare::Vector<Dim, SC> dn = g_x->GetDistances();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_pressure_force_Cartesian(&pm, dare::ZERO, ind, &epsilon);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_low(ind);
        ind_low[0] -= 1;
        SC dp = pm.GetPressure()->GetDataVector().At(ind, 0) - pm.GetPressure()->GetDataVector().At(ind_low, 0);
        SC grad_p = dp / dn[0];
        SC eps_f = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_low, 0));
        SC e_grad_p = -eps_f * grad_p * dV;
        EXPECT_NEAR(s.GetRhs(0), e_grad_p, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(e_grad_p));
    }

    dV = g_y->GetCellVolume();
    dn = g_y->GetDistances();
    // in Y-momentum
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_pressure_force_Cartesian(&pm, dare::ONE, ind, &epsilon);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_low(ind);
        ind_low[1] -= 1;
        SC dp = pm.GetPressure()->GetDataVector().At(ind, 0) - pm.GetPressure()->GetDataVector().At(ind_low, 0);
        SC grad_p = dp / dn[1];
        SC eps_f = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_low, 0));
        SC e_grad_p = -eps_f * grad_p * dV;
        EXPECT_NEAR(s.GetRhs(0), e_grad_p, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(e_grad_p));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, BuildMomentum_pressure_test) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    using NDict = dare::PMNumericalInfoDefault;

    SC tol_eps = 1e2;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    double mu = 0;
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    auto g_z = &pm.GetMomentum(2)->GetField()->GetGridRepresentation();
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    pm.SetPorosity(&epsilon);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < epsilon.GetDataVector().GetSize(); i++) {
        epsilon.GetDataVector().At(i) = std::abs(rd.Generate());
        epsilon.GetDataVector(1).At(i) = std::abs(rd.Generate());
        pm.GetPressure()->GetDataVector().At(i) = rd.Generate();
        pm.GetPressure()->GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    SC dV = g_x->GetCellVolume();
    dare::Vector<Dim, SC> dn = g_x->GetDistances();
    // in X-momentum
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_pressure_force_Cartesian(&pm, dare::ZERO, ind, &epsilon);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_low(ind);
        ind_low[0] -= 1;
        SC dp = pm.GetPressure()->GetDataVector().At(ind, 0) - pm.GetPressure()->GetDataVector().At(ind_low, 0);
        SC grad_p = dp / dn[0];
        SC eps_f = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_low, 0));
        SC e_grad_p = -eps_f * grad_p * dV;
        EXPECT_NEAR(s.GetRhs(0), e_grad_p, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(e_grad_p));
    }

    dV = g_y->GetCellVolume();
    dn = g_y->GetDistances();
    // in Y-momentum
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_pressure_force_Cartesian(&pm, dare::ONE, ind, &epsilon);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_low(ind);
        ind_low[1] -= 1;
        SC dp = pm.GetPressure()->GetDataVector().At(ind, 0) - pm.GetPressure()->GetDataVector().At(ind_low, 0);
        SC grad_p = dp / dn[1];
        SC eps_f = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_low, 0));
        SC e_grad_p = -eps_f * grad_p * dV;
        EXPECT_NEAR(s.GetRhs(0), e_grad_p, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(e_grad_p));
    }

    dV = g_z->GetCellVolume();
    dn = g_z->GetDistances();
    // in Z-momentum
    for (LO n_loc = 0; n_loc < g_z->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_pressure_force_Cartesian(&pm, dare::TWO, ind, &epsilon);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        Index ind_low(ind);
        ind_low[2] -= 1;
        SC dp = pm.GetPressure()->GetDataVector().At(ind, 0) - pm.GetPressure()->GetDataVector().At(ind_low, 0);
        SC grad_p = dp / dn[2];
        SC eps_f = 0.5 * (epsilon.GetDataVector().At(ind, 0) + epsilon.GetDataVector().At(ind_low, 0));
        SC e_grad_p = -eps_f * grad_p * dV;
        EXPECT_NEAR(s.GetRhs(0), e_grad_p, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(e_grad_p));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_explicit_force_test) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = dare::None;
        using explicit_force = Field;
        using implicit_force = dare::None;
    };
    using NDict = dare::PMNumericalInfoDefault;

    SC tol_eps = 1e2;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    double mu = 0;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    Field beta_ex("beta_ex", *g_x, 1);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    pm.AddExplicitForce(&beta_ex, 0);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            beta_ex.GetDataVector().At(i) = rd.Generate();
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < pm.GetPressure()->GetDataVector().GetSize(); i++) {
        pm.GetPressure()->GetDataVector().At(i) = rd.Generate();
        pm.GetPressure()->GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    // in X-momentum
    SC dV = g_x->GetCellVolume();
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_explicit_force_Cartesian(&pm, dare::ZERO, ind);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        SC beta = beta_ex.GetDataVector().At(ind, 0) * dV;
        EXPECT_NEAR(s.GetRhs(0), beta, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(beta));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, BuildMomentum_explicit_force_test) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = dare::None;
        using explicit_force = Field;
        using implicit_force = dare::None;
    };
    using NDict = dare::PMNumericalInfoDefault;

    SC tol_eps = 1e2;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    double mu = 0;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    Field beta_x("beta_x", *g_x, 1);
    Field beta_y("beta_y", *g_y, 1);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    pm.AddExplicitForce(&beta_x, 0);
    pm.AddExplicitForce(&beta_y, 1);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            if (d == 0)
                beta_x.GetDataVector().At(i) = rd.Generate();
            else if (d == 1)
                beta_y.GetDataVector().At(i) = rd.Generate();
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < pm.GetPressure()->GetDataVector().GetSize(); i++) {
        pm.GetPressure()->GetDataVector().At(i) = rd.Generate();
        pm.GetPressure()->GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    // in X-momentum
    SC dV = g_x->GetCellVolume();
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_explicit_force_Cartesian(&pm, dare::ZERO, ind);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        SC beta = beta_x.GetDataVector().At(ind, 0) * dV;
        EXPECT_NEAR(s.GetRhs(0), beta, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(beta));
    }

    // in Y-momentum
    dV = g_y->GetCellVolume();
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_explicit_force_Cartesian(&pm, dare::ONE, ind);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        SC beta = beta_y.GetDataVector().At(ind, 0) * dV;
        EXPECT_NEAR(s.GetRhs(0), beta, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(beta));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, BuildMomentum_explicit_force_test) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = dare::None;
        using explicit_force = Field;
        using implicit_force = dare::None;
    };
    using NDict = dare::PMNumericalInfoDefault;

    SC tol_eps = 1e2;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    double mu = 0;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    auto g_y = &pm.GetMomentum(1)->GetField()->GetGridRepresentation();
    auto g_z = &pm.GetMomentum(2)->GetField()->GetGridRepresentation();
    Field beta_x("beta_x", *g_x, 1);
    Field beta_y("beta_y", *g_y, 1);
    Field beta_z("beta_z", *g_z, 1);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    pm.AddExplicitForce(&beta_x, 0);
    pm.AddExplicitForce(&beta_y, 1);
    pm.AddExplicitForce(&beta_z, 2);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            if (d == 0)
                beta_x.GetDataVector().At(i) = rd.Generate();
            else if (d == 1)
                beta_y.GetDataVector().At(i) = rd.Generate();
            else if (d == 2)
                beta_z.GetDataVector().At(i) = rd.Generate();
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }
    for (std::size_t i{0}; i < pm.GetPressure()->GetDataVector().GetSize(); i++) {
        pm.GetPressure()->GetDataVector().At(i) = rd.Generate();
        pm.GetPressure()->GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    // in X-momentum
    SC dV = g_x->GetCellVolume();
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_x->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_explicit_force_Cartesian(&pm, dare::ZERO, ind);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        SC beta = beta_x.GetDataVector().At(ind, 0) * dV;
        EXPECT_NEAR(s.GetRhs(0), beta, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(beta));
    }

    // in Y-momentum
    dV = g_y->GetCellVolume();
    for (LO n_loc = 0; n_loc < g_y->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_y->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_y->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_explicit_force_Cartesian(&pm, dare::ONE, ind);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        SC beta = beta_y.GetDataVector().At(ind, 0) * dV;
        EXPECT_NEAR(s.GetRhs(0), beta, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(beta));
    }

    // in Z-momentum
    dV = g_z->GetCellVolume();
    for (LO n_loc = 0; n_loc < g_z->GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_z->MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_z->MapInternalToLocal(ind_loc);
        auto s = dare::free_pm_explicit_force_Cartesian(&pm, dare::TWO, ind);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);
        SC beta = beta_z.GetDataVector().At(ind, 0) * dV;
        EXPECT_NEAR(s.GetRhs(0), beta, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(beta));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_normalization_test) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = dare::None;
        using explicit_force = Field;
        using implicit_force = dare::None;
    };
    struct NDict {
        using momentum_normalizer = double;
        using continuity_normalizer = double;
    };
    using CNB = dare::CartesianNeighbor;
    auto g_s = grid->GetRepresentation(opt_s);
    double rho = 0.;
    double mu = 0;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();
    Field beta_x("beta_x", *g_x, 1);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);

    SC w_base{-1.}, e_base{1.}, c_base{1.2};
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        dare::MatrixBlock<GridType, LO, SC, 1> mblock(g_x, n_loc, dare::Vector<1, SC>{3});
        SC normalizer = rd.Generate();
        mblock.Get(0, 0, CNB::CENTER) = c_base;
        mblock.Get(0, 0, CNB::WEST) = w_base;
        mblock.Get(0, 0, CNB::EAST) = e_base;
        pm.GetMomentum(0)->GetCustomMember()->normalizer = normalizer;
        dare::free_pm_apply_normalizer(&pm, pm.GetMomentum(0)->GetCustomMember()->normalizer, &mblock);
        EXPECT_EQ(mblock.Get(0, 0, CNB::CENTER), c_base * normalizer);
        EXPECT_EQ(mblock.Get(0, 0, CNB::WEST), w_base * normalizer);
        EXPECT_EQ(mblock.Get(0, 0, CNB::EAST), e_base * normalizer);
    }
}
