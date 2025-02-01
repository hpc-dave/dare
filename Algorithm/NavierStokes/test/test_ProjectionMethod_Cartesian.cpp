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
        using density = double;
        using viscosity = double;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    double rho = 1.;
    double mu = 0.;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);
    dare::ProjectionMethod<GridType, BStrat, PDict, NDefault> pm;
    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++)
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
    }
    pm.GetContinuity()->GetPressure()->SetValues(0.);
    auto g_x = &pm.GetMomentum(0)->GetField()->GetGridRepresentation();

    Index ind_loc(1);
    Index ind = g_x->MapInternalToLocal(ind_loc);
    LO o_loc = g_x->MapIndexToOrdinalLocalInternal(ind_loc);
    auto s = dare::free_pm_ddt_Cartesian(&pm, dare::ZERO, *g_x, o_loc);

    static_assert(std::is_same_v<decltype(s), dare::CenterMatrixStencil<GridType, SC, 1>>);

    SC dV_dt = g_x->GetCellVolume()/dt;
    SC u = pm.GetMomentum(0)->GetField()->GetDataVector(1).At(ind, 0);
    EXPECT_EQ(s.Center(0), rho * dV_dt);
    EXPECT_EQ(s.GetRhs(0), rho * dV_dt * u);
}
