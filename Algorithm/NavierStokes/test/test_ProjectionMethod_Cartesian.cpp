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

#include <gtest/gtest.h>

#include <array>
#include <type_traits>

#include "Algorithm/NavierStokes/ProjectionMethod.h"
#include "Grid/Cartesian.h"
#include "Algorithm/NavierStokes/ProjectionMethod_Cartesian.h"
#include "Algorithm/ConstantTimeStep.h"
#include "Utilities/RandomNumberGenerator.h"


namespace dare::test {

template <std::size_t Dim>
dare::Vector<Dim, defaults::GlobalOrdinalType> GetResolutionTestPMCartesian() {
    dare::Vector<Dim, defaults::GlobalOrdinalType> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 10 + n;
    return res;
}
template <std::size_t Dim>
dare::Vector<Dim, defaults::ScalarType> GetSizeTestPMCartesian() {
    dare::Vector<Dim, defaults::ScalarType> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

// a dummy for the boundary strategy
struct BStrat {
    template <typename T>
    void operator()(T t) {}
};

}  // namespace dare::test

/*!
 * @brief Fixture for testing the interpolation functions with the Cartesian Grid
 */
template<std::size_t D>
class ProjectionMethodCartesianTest : public testing::Test {
public:
    static const std::size_t Dim{D};
    static const std::size_t N{1};
    using GridType = dare::Cartesian<Dim>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using Index = typename GridType::Index;
    using GridVector = dare::GridVector<GridType, SC, N>;
    using Field = dare::Field<GridType, SC, N>;
    using Options = typename GridType::Options;

    using PDefault = dare::PMPropertyInfoDefault<GridType>;
    using NDefault = dare::PMNumericalInfoDefault;
    using BStrat = dare::test::BStrat;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(&exec_man,
                                          dare::test::GetResolutionTestPMCartesian<Dim>(),
                                          dare::test::GetSizeTestPMCartesian<Dim>(),
                                          num_ghost);
        for (std::size_t d{0}; d < Dim; d++) {
            opt_s[d] = 0;
            for (std::size_t e{0}; e < Dim; e++) {
                opt_m[d][e] = e == d;
            }
        }
    }

    std::unique_ptr<GridType> grid;   //!< the grid
    dare::ExecutionManager exec_man;  //!< the execution manager
    Options opt_s;                    //!< options scalar grid
    std::array<Options, Dim> opt_m;   //!< options momentum grid
};

using ProjectionMethodCartesian1DTest = ProjectionMethodCartesianTest<1>;
using ProjectionMethodCartesian2DTest = ProjectionMethodCartesianTest<2>;
using ProjectionMethodCartesian3DTest = ProjectionMethodCartesianTest<3>;

TEST_F(ProjectionMethodCartesian1DTest, Initialization) {
    dare::ProjectionMethod<GridType, BStrat, PDefault, NDefault> pm;
    dare::test::BStrat bstrat;
    dare::ConstantTimeStep dt(1.);
    EXPECT_FALSE(pm.IsInitialized());
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    const double tol_eps = 1e2;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
    pm.Initialize(grid, &dt, bstrat);
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
