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

TEST_F(ProjectionMethodCartesian1DTest, BuildMomentum_conv_test) {
    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = Field;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::MINMOD;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };
    using TVD = dare::TVD<GridType, SC, dare::MINMOD>;
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    // const double tol_eps = 1e2;
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
    SC v_step{1.};
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            Index ind = g_x->MapOrdinalToIndexLocal(i);
            SC vloc = ind[d] * v_step;
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = vloc;
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

    // dare::Vector<Dim, SC> dA = g_x->GetFaceArea();
    // in x-direction
    for (LO n_loc = 0; n_loc < g_x->GetNumberLocalCellsInternal(); n_loc++) {
        // Index ind_loc = g_x->MapOrdinalToIndexLocalInternal(n_loc);
        // Index ind = g_x->MapInternalToLocal(ind_loc);

        TVD tvd(*g_x, n_loc, velocities);

        auto s = dare::free_pm_convection_Cartesian(&pm, dare::ZERO, *g_x, n_loc, &epsilon, &rho, velocities);

        static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

        // SC u = pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0);
        // Index ind_nb(ind);
        // ind_nb.i() -= 1;
        FVStencil rho_s = tvd.Interpolate(rho.GetDataVector(1));
        FVStencil eps_s = tvd.Interpolate(epsilon.GetDataVector(1));
        FVStencil uloc = tvd.Interpolate(pm.GetMomentum(0)->GetField()->GetDataVector(1));

        FVStencil mom_loc = eps_s * rho_s * uloc;

        // SC rho_w = rho.GetDataVector(0).At(ind_nb, 0);
        // SC eps_e = epsilon.GetDataVector(0).At(ind, 0);
        // SC eps_w = epsilon.GetDataVector(0).At(ind_nb, 0);
        // SC v_0 = eps_0 * rho_0 * dV_dt;
        // SC v_1 = eps_1 * rho_1 * dV_dt * u;
        // EXPECT_NEAR(s.Center(0), v_0, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_0));
        // EXPECT_NEAR(s.GetRhs(0), v_1, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(v_1));
        // for (auto face : g_x->GetFaces())
        //     EXPECT_EQ(s.GetValue(face, 0), 0.);
    }
}
