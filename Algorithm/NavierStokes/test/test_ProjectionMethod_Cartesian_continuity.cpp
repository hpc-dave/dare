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

TEST_F(ProjectionMethodCartesian1DTest, ContinuityDefectIncompressible) {
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_incompressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                               + epsilon.GetDataVector().At(ind_low, 0)) * velocities[d]->At(ind, 0);
            SC flux_hi = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                              + epsilon.GetDataVector().At(ind_hi, 0)) * velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, ContinuityDefectIncompressible_no_porosity) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = dare::None;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    double rho = 0.;
    const double tol_eps = 1e2;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_incompressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, ContinuityDefectIncompressible) {
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
    const double tol_eps = 1e3;
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_incompressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                               + epsilon.GetDataVector().At(ind_low, 0)) * velocities[d]->At(ind, 0);
            SC flux_hi = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                              + epsilon.GetDataVector().At(ind_hi, 0)) * velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, ContinuityDefectIncompressible_no_porosity) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = dare::None;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    double rho = 0.;
    const double tol_eps = 1e3;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_incompressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, ContinuityDefectIncompressible) {
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
    const double tol_eps = 1e4;
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_incompressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                               + epsilon.GetDataVector().At(ind_low, 0)) * velocities[d]->At(ind, 0);
            SC flux_hi = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                              + epsilon.GetDataVector().At(ind_hi, 0)) * velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, ContinuityDefectIncompressible_no_porosity) {
    struct PDict {
        using density = double;
        using viscosity = double;
        using porosity = dare::None;
        using explicit_force = dare::None;
        using implicit_force = dare::None;
    };
    struct NDict {
        using tvd = dare::UPWIND;
        using time_scheme_convective = dare::EULER_BACKWARD;
    };

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    double rho = 0.;
    const double tol_eps = 1e4;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(!pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(rho);
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_incompressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, Continuity_compute_defect_incompressible) {
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    free_pm_compute_defect(&pm);

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = pm.GetContinuity()->GetDefect()->GetDataVector().At(ind, 0);

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                               + epsilon.GetDataVector().At(ind_low, 0)) * velocities[d]->At(ind, 0);
            SC flux_hi = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                              + epsilon.GetDataVector().At(ind_hi, 0)) * velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, Continuity_compute_defect_incompressible) {
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
    const double tol_eps = 1e3;
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(0);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    free_pm_compute_defect(&pm);

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = pm.GetContinuity()->GetDefect()->GetDataVector().At(ind, 0);

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                               + epsilon.GetDataVector().At(ind_low, 0)) * velocities[d]->At(ind, 0);
            SC flux_hi = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                              + epsilon.GetDataVector().At(ind_hi, 0)) * velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, Continuity_compute_defect_incompressible) {
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
    const double tol_eps = 1e4;
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
            pm.GetMomentum(d)->GetField()->GetDataVector(0).At(i) = rd.Generate();
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

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    free_pm_compute_defect(&pm);

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = pm.GetContinuity()->GetDefect()->GetDataVector().At(ind, 0);

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                               + epsilon.GetDataVector().At(ind_low, 0)) * velocities[d]->At(ind, 0);
            SC flux_hi = 0.5 * (epsilon.GetDataVector().At(ind, 0)
                              + epsilon.GetDataVector().At(ind_hi, 0)) * velocities[d]->At(ind_hi, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, ContinuityDefectCompressible_no_porosity) {
    SC dd_val{1.12};
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.): v(dd) {}
        SC operator()(Index ind) { return 1.12; }
        SC v;
    };

    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = dare::None;
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
    Field rho = Field("rho", g_s, 2);
    const double tol_eps = 1e2;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(0).At(i) = rd.Generate();
        }
    }

    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector(0).At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(0);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.)
                flux_low *= rho.GetDataVector().At(ind, 0);
            else
                flux_low *= rho.GetDataVector().At(ind_low, 0);
            if (velocities[d]->At(ind_hi, 0) < 0.)
                flux_hi *= rho.GetDataVector().At(ind_hi, 0);
            else
                flux_hi *= rho.GetDataVector().At(ind, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        SC rho_c = rho.GetDataVector().At(ind, 0);
        SC rho_c_old = rho.GetDataVector(1).At(ind, 0);
        d_ex += (rho_c - rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, ContinuityDefectCompressible) {
    SC dd_val{0.};
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return 1.12; }
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
    double rho = 3.14;
    Field rho_f = Field("rho", g_s, 2);
    Field epsilon = Field("epsilon", g_s, 2);
    const double tol_eps = 1e2;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(&rho_f);
    pm.SetPorosity(&epsilon);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    for (std::size_t i{0}; i < epsilon.GetDataVector().GetSize(); i++) {
        epsilon.GetDataVector(0).At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        rho_f.GetDataVector(0).At(i) = rho;
        rho_f.GetDataVector(1).At(i) = rho;
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.)
                flux_low *= epsilon.GetDataVector().At(ind, 0);
            else
                flux_low *= epsilon.GetDataVector().At(ind_low, 0);
            if (velocities[d]->At(ind_hi, 0) < 0.)
                flux_hi *= epsilon.GetDataVector().At(ind_hi, 0);
            else
                flux_hi *= epsilon.GetDataVector().At(ind, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        d_ex *= rho;
        SC eps_rho_c = epsilon.GetDataVector().At(ind, 0) * rho;
        SC eps_rho_c_old = epsilon.GetDataVector(1).At(ind, 0) * rho;
        d_ex += (eps_rho_c - eps_rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, ContinuityDefectCompressible_no_porosity) {
    SC dd_val{1.12};
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return 1.12; }
        SC v;
    };

    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = dare::None;
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
    Field rho = Field("rho", g_s, 2);
    const double tol_eps = 1e2;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector(0).At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.)
                flux_low *= rho.GetDataVector().At(ind, 0);
            else
                flux_low *= rho.GetDataVector().At(ind_low, 0);
            if (velocities[d]->At(ind_hi, 0) < 0.)
                flux_hi *= rho.GetDataVector().At(ind_hi, 0);
            else
                flux_hi *= rho.GetDataVector().At(ind, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        SC rho_c = rho.GetDataVector().At(ind, 0);
        SC rho_c_old = rho.GetDataVector(1).At(ind, 0);
        d_ex += (rho_c - rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, ContinuityDefectCompressible) {
    SC dd_val{0.};
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return 1.12; }
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
    double rho = 3.14;
    Field rho_f = Field("rho", g_s, 2);
    Field epsilon = Field("epsilon", g_s, 2);
    const double tol_eps = 1e2;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho_f);
    pm.SetPorosity(&epsilon);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    for (std::size_t i{0}; i < epsilon.GetDataVector().GetSize(); i++) {
        epsilon.GetDataVector(0).At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        rho_f.GetDataVector(0).At(i) = rho;
        rho_f.GetDataVector(1).At(i) = rho;
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.)
                flux_low *= epsilon.GetDataVector().At(ind, 0);
            else
                flux_low *= epsilon.GetDataVector().At(ind_low, 0);
            if (velocities[d]->At(ind_hi, 0) < 0.)
                flux_hi *= epsilon.GetDataVector().At(ind_hi, 0);
            else
                flux_hi *= epsilon.GetDataVector().At(ind, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        d_ex *= rho;
        SC eps_rho_c = epsilon.GetDataVector().At(ind, 0) * rho;
        SC eps_rho_c_old = epsilon.GetDataVector(1).At(ind, 0) * rho;
        d_ex += (eps_rho_c - eps_rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, ContinuityDefectCompressible_no_porosity) {
    SC dd_val{1.12};
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return 1.12; }
        SC v;
    };

    struct PDict {
        using density = Field;
        using viscosity = double;
        using porosity = dare::None;
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
    Field rho = Field("rho", g_s, 2);
    const double tol_eps = 1e4;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    for (std::size_t i{0}; i < rho.GetDataVector().GetSize(); i++) {
        rho.GetDataVector(0).At(i) = rd.Generate();
        rho.GetDataVector(1).At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.)
                flux_low *= rho.GetDataVector().At(ind, 0);
            else
                flux_low *= rho.GetDataVector().At(ind_low, 0);
            if (velocities[d]->At(ind_hi, 0) < 0.)
                flux_hi *= rho.GetDataVector().At(ind_hi, 0);
            else
                flux_hi *= rho.GetDataVector().At(ind, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        SC rho_c = rho.GetDataVector().At(ind, 0);
        SC rho_c_old = rho.GetDataVector(1).At(ind, 0);
        d_ex += (rho_c - rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, ContinuityDefectCompressible) {
    SC dd_val{0.};
    struct density_derivative_functor {
        explicit density_derivative_functor(SC dd = 0.) : v(dd) {}
        SC operator()(Index ind) { return 1.12; }
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
    double rho = 3.14;
    Field rho_f = Field("rho", g_s, 2);
    Field epsilon = Field("epsilon", g_s, 2);
    const double tol_eps = 1e4;
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho_f);
    pm.SetPorosity(&epsilon);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
    pm.SetViscosity(mu);
    for (std::size_t d{0}; d < Dim; d++) {
        for (std::size_t i{0}; i < pm.GetMomentum(d)->GetField()->GetDataVector().GetSize(); i++) {
            pm.GetMomentum(d)->GetField()->GetDataVector(1).At(i) = rd.Generate();
        }
    }

    for (std::size_t i{0}; i < epsilon.GetDataVector().GetSize(); i++) {
        epsilon.GetDataVector(0).At(i) = rd.Generate();
        epsilon.GetDataVector(1).At(i) = rd.Generate();
        rho_f.GetDataVector(0).At(i) = rho;
        rho_f.GetDataVector(1).At(i) = rho;
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.)
                flux_low *= epsilon.GetDataVector().At(ind, 0);
            else
                flux_low *= epsilon.GetDataVector().At(ind_low, 0);
            if (velocities[d]->At(ind_hi, 0) < 0.)
                flux_hi *= epsilon.GetDataVector().At(ind_hi, 0);
            else
                flux_hi *= epsilon.GetDataVector().At(ind, 0);
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        d_ex *= rho;
        SC eps_rho_c = epsilon.GetDataVector().At(ind, 0) * rho;
        SC eps_rho_c_old = epsilon.GetDataVector(1).At(ind, 0) * rho;
        d_ex += (eps_rho_c - eps_rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, Continuity_compute_defect_compressible) {
    SC dd_val{0.};
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
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
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

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    free_pm_compute_defect(&pm);

    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.) {
                flux_low *= epsilon.GetDataVector().At(ind, 0);
                flux_low *= rho.GetDataVector().At(ind, 0);
            } else {
                flux_low *= epsilon.GetDataVector().At(ind_low, 0);
                flux_low *= rho.GetDataVector().At(ind_low, 0);
            }
            if (velocities[d]->At(ind_hi, 0) < 0.) {
                flux_hi *= epsilon.GetDataVector().At(ind_hi, 0);
                flux_hi *= rho.GetDataVector().At(ind_hi, 0);
            } else {
                flux_hi *= epsilon.GetDataVector().At(ind, 0);
                flux_hi *= rho.GetDataVector().At(ind, 0);
            }
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        SC eps_rho_c = epsilon.GetDataVector().At(ind, 0);
        eps_rho_c *= rho.GetDataVector().At(ind, 0);
        SC eps_rho_c_old = epsilon.GetDataVector(1).At(ind, 0);
        eps_rho_c_old *= rho.GetDataVector(1).At(ind, 0);
        d_ex += (eps_rho_c - eps_rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, Continuity_compute_defect_compressible) {
    SC dd_val{0.};
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
    const double tol_eps = 1e3;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
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

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    free_pm_compute_defect(&pm);

    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.) {
                flux_low *= epsilon.GetDataVector().At(ind, 0);
                flux_low *= rho.GetDataVector().At(ind, 0);
            } else {
                flux_low *= epsilon.GetDataVector().At(ind_low, 0);
                flux_low *= rho.GetDataVector().At(ind_low, 0);
            }
            if (velocities[d]->At(ind_hi, 0) < 0.) {
                flux_hi *= epsilon.GetDataVector().At(ind_hi, 0);
                flux_hi *= rho.GetDataVector().At(ind_hi, 0);
            } else {
                flux_hi *= epsilon.GetDataVector().At(ind, 0);
                flux_hi *= rho.GetDataVector().At(ind, 0);
            }
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        SC eps_rho_c = epsilon.GetDataVector().At(ind, 0);
        eps_rho_c *= rho.GetDataVector().At(ind, 0);
        SC eps_rho_c_old = epsilon.GetDataVector(1).At(ind, 0);
        eps_rho_c_old *= rho.GetDataVector(1).At(ind, 0);
        d_ex += (eps_rho_c - eps_rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, Continuity_compute_defect_compressible) {
    SC dd_val{0.};
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
    const double tol_eps = 1e4;
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
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

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();

    free_pm_compute_defect(&pm);

    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        SC d = free_pm_defect_compressible_Cartesian(&pm, n_loc, ind, velocities, pm.GetDensity(), pm.GetPorosity());

        SC d_ex{0.};
        for (std::size_t d{0}; d < Dim; d++) {
            Index ind_low(ind), ind_hi(ind);
            ind_low[d] -= 1;
            ind_hi[d] += 1;
            SC flux_low = velocities[d]->At(ind, 0);
            SC flux_hi = velocities[d]->At(ind_hi, 0);
            if (velocities[d]->At(ind, 0) < 0.) {
                flux_low *= epsilon.GetDataVector().At(ind, 0);
                flux_low *= rho.GetDataVector().At(ind, 0);
            } else {
                flux_low *= epsilon.GetDataVector().At(ind_low, 0);
                flux_low *= rho.GetDataVector().At(ind_low, 0);
            }
            if (velocities[d]->At(ind_hi, 0) < 0.) {
                flux_hi *= epsilon.GetDataVector().At(ind_hi, 0);
                flux_hi *= rho.GetDataVector().At(ind_hi, 0);
            } else {
                flux_hi *= epsilon.GetDataVector().At(ind, 0);
                flux_hi *= rho.GetDataVector().At(ind, 0);
            }
            d_ex += dA[d] * (flux_hi - flux_low);
        }
        SC eps_rho_c = epsilon.GetDataVector().At(ind, 0);
        eps_rho_c *= rho.GetDataVector().At(ind, 0);
        SC eps_rho_c_old = epsilon.GetDataVector(1).At(ind, 0);
        eps_rho_c_old *= rho.GetDataVector(1).At(ind, 0);
        d_ex += (eps_rho_c - eps_rho_c_old) * dV / dt;
        EXPECT_NEAR(d, d_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(d_ex));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, ContinuityJacobianIncompressible) {
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
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    dare::Vector<Dim, SC> dn_r = 1./g_s.GetDistances();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        auto J = free_pm_continuity_Jacobian_Cartesian(&pm, n_loc, ind, pm.GetDensity(), pm.GetPorosity());
        static_assert(std::is_same_v<decltype(J), dare::CenterMatrixStencil<GridType, SC, 1>>,
                      "Should be a center matrix stencil");

        Index ind_w(ind), ind_e(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        FVStencil coef_faces;
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC beta_w = 0.5 * (beta.GetDataVector().At(ind_w, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_e = 0.5 * (beta.GetDataVector().At(ind_e, 0) + beta.GetDataVector().At(ind, 0));

        coef_faces(CNB::WEST, 0) = -dt * eps_w / (rho_w + beta_w * dt / eps_w);
        coef_faces(CNB::EAST, 0) = -dt * eps_e / (rho_e + beta_e * dt / eps_e);

        SC c_ex{0.};
        for (auto face : g_s.GetFaces()) {
            coef_faces(face, 0) *= dA[dare::MapCartesianFaceToDim(face)] * dn_r[dare::MapCartesianFaceToDim(face)];
            c_ex += -coef_faces(face, 0);
        }

        EXPECT_NEAR(J.GetValue(CNB::WEST, 0), coef_faces(CNB::WEST, 0),
                      tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::WEST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::EAST, 0), coef_faces(CNB::EAST, 0),
                      tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::EAST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::CENTER, 0), c_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(c_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, ContinuityJacobianIncompressible) {
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
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        auto J = free_pm_continuity_Jacobian_Cartesian(&pm, n_loc, ind, pm.GetDensity(), pm.GetPorosity());
        static_assert(std::is_same_v<decltype(J), dare::CenterMatrixStencil<GridType, SC, 1>>,
                      "Should be a center matrix stencil");

        Index ind_w(ind), ind_e(ind), ind_s(ind), ind_n(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_n.j() += 1;
        FVStencil coef_faces;
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_s, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_n, 0) + epsilon.GetDataVector().At(ind, 0));
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_s, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_n, 0) + rho.GetDataVector().At(ind, 0));
        SC beta_w = 0.5 * (beta.GetDataVector().At(ind_w, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_e = 0.5 * (beta.GetDataVector().At(ind_e, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_s = 0.5 * (beta.GetDataVector().At(ind_s, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_n = 0.5 * (beta.GetDataVector().At(ind_n, 0) + beta.GetDataVector().At(ind, 0));

        coef_faces(CNB::WEST, 0) = -dt * eps_w / (rho_w + beta_w * dt / eps_w);
        coef_faces(CNB::EAST, 0) = -dt * eps_e / (rho_e + beta_e * dt / eps_e);
        coef_faces(CNB::SOUTH, 0) = -dt * eps_s / (rho_s + beta_s * dt / eps_s);
        coef_faces(CNB::NORTH, 0) = -dt * eps_n / (rho_n + beta_n * dt / eps_n);

        SC c_ex{0.};
        for (auto face : g_s.GetFaces()) {
            coef_faces(face, 0) *= dA[dare::MapCartesianFaceToDim(face)] * dn_r[dare::MapCartesianFaceToDim(face)];
            c_ex += -coef_faces(face, 0);
        }

        EXPECT_NEAR(J.GetValue(CNB::WEST, 0), coef_faces(CNB::WEST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::WEST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::EAST, 0), coef_faces(CNB::EAST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::EAST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::SOUTH, 0), coef_faces(CNB::SOUTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::SOUTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::NORTH, 0), coef_faces(CNB::NORTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::NORTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::CENTER, 0), c_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(c_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, ContinuityJacobianIncompressible) {
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
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e4;
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
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();

    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        auto J = free_pm_continuity_Jacobian_Cartesian(&pm, n_loc, ind, pm.GetDensity(), pm.GetPorosity());
        static_assert(std::is_same_v<decltype(J), dare::CenterMatrixStencil<GridType, SC, 1>>,
                      "Should be a center matrix stencil");

        Index ind_w(ind), ind_e(ind), ind_s(ind), ind_n(ind), ind_b(ind), ind_t(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_n.j() += 1;
        ind_b.k() -= 1;
        ind_t.k() += 1;
        FVStencil coef_faces;
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_s, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_n, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_b = 0.5 * (epsilon.GetDataVector().At(ind_b, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_t = 0.5 * (epsilon.GetDataVector().At(ind_t, 0) + epsilon.GetDataVector().At(ind, 0));
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_s, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_n, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_b = 0.5 * (rho.GetDataVector().At(ind_b, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_t = 0.5 * (rho.GetDataVector().At(ind_t, 0) + rho.GetDataVector().At(ind, 0));
        SC beta_w = 0.5 * (beta.GetDataVector().At(ind_w, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_e = 0.5 * (beta.GetDataVector().At(ind_e, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_s = 0.5 * (beta.GetDataVector().At(ind_s, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_n = 0.5 * (beta.GetDataVector().At(ind_n, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_b = 0.5 * (beta.GetDataVector().At(ind_b, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_t = 0.5 * (beta.GetDataVector().At(ind_t, 0) + beta.GetDataVector().At(ind, 0));

        coef_faces(CNB::WEST, 0) = -dt * eps_w / (rho_w + beta_w * dt / eps_w);
        coef_faces(CNB::EAST, 0) = -dt * eps_e / (rho_e + beta_e * dt / eps_e);
        coef_faces(CNB::SOUTH, 0) = -dt * eps_s / (rho_s + beta_s * dt / eps_s);
        coef_faces(CNB::NORTH, 0) = -dt * eps_n / (rho_n + beta_n * dt / eps_n);
        coef_faces(CNB::BOTTOM, 0) = -dt * eps_b / (rho_b + beta_b * dt / eps_b);
        coef_faces(CNB::TOP, 0) = -dt * eps_t / (rho_t + beta_t * dt / eps_t);

        SC c_ex{0.};
        for (auto face : g_s.GetFaces()) {
            coef_faces(face, 0) *= dA[dare::MapCartesianFaceToDim(face)] * dn_r[dare::MapCartesianFaceToDim(face)];
            c_ex += -coef_faces(face, 0);
        }

        EXPECT_NEAR(J.GetValue(CNB::WEST, 0), coef_faces(CNB::WEST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::WEST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::EAST, 0), coef_faces(CNB::EAST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::EAST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::SOUTH, 0), coef_faces(CNB::SOUTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::SOUTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::NORTH, 0), coef_faces(CNB::NORTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::NORTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::BOTTOM, 0), coef_faces(CNB::BOTTOM, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::BOTTOM, 0)));
        EXPECT_NEAR(J.GetValue(CNB::TOP, 0), coef_faces(CNB::TOP, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::TOP, 0)));
        EXPECT_NEAR(J.GetValue(CNB::CENTER, 0), c_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(c_ex));
    }
}

TEST_F(ProjectionMethodCartesian1DTest, ContinuityJacobianCompressible) {
    struct density_derivative_functor {
        density_derivative_functor() : field(nullptr) {}
        explicit density_derivative_functor(const Field& dd) : field(&dd) {}
        SC operator()(Index ind) const { return field->GetDataVector().At(ind, 0); }
        const Field* field;
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
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e5;
    Field dd_val("dd", g_s, 1);    // value of density derivative
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
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
        dd_val.GetDataVector().At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        auto J = free_pm_continuity_Jacobian_Cartesian(&pm, n_loc, ind, pm.GetDensity(), pm.GetPorosity());
        static_assert(std::is_same_v<decltype(J), dare::CenterMatrixStencil<GridType, SC, 1>>,
                      "Should be a center matrix stencil");

        Index ind_w(ind), ind_e(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        FVStencil coef_faces;
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC beta_w = 0.5 * (beta.GetDataVector().At(ind_w, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_e = 0.5 * (beta.GetDataVector().At(ind_e, 0) + beta.GetDataVector().At(ind, 0));

        coef_faces(CNB::WEST, 0) = -dt * eps_w / (1. + beta_w * dt / (eps_w * rho_w));
        coef_faces(CNB::EAST, 0) = -dt * eps_e / (1. +  beta_e * dt / (eps_e * rho_e));

        SC c_ex{0.};
        for (auto face : g_s.GetFaces()) {
            coef_faces(face, 0) *= dA[dare::MapCartesianFaceToDim(face)] * dn_r[dare::MapCartesianFaceToDim(face)];
            c_ex += -coef_faces(face, 0);
        }
        SC dd = dd_val.GetDataVector().At(ind, 0);
        c_ex += dd * epsilon.GetDataVector().At(ind, 0) * dV / dt;

        EXPECT_NEAR(J.GetValue(CNB::WEST, 0), coef_faces(CNB::WEST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::WEST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::EAST, 0), coef_faces(CNB::EAST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::EAST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::CENTER, 0), c_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(c_ex));
    }
}

TEST_F(ProjectionMethodCartesian2DTest, ContinuityJacobianCompressible) {
    struct density_derivative_functor {
        density_derivative_functor() : field(nullptr) {}
        explicit density_derivative_functor(const Field& dd) : field(&dd) {}
        SC operator()(Index ind) const { return field->GetDataVector().At(ind, 0); }
        const Field* field;
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
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e4;
    Field dd_val("dd", g_s, 1);  // value of density derivative
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
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
        dd_val.GetDataVector().At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        auto J = free_pm_continuity_Jacobian_Cartesian(&pm, n_loc, ind, pm.GetDensity(), pm.GetPorosity());
        static_assert(std::is_same_v<decltype(J), dare::CenterMatrixStencil<GridType, SC, 1>>,
                      "Should be a center matrix stencil");

        Index ind_w(ind), ind_e(ind), ind_s(ind), ind_n(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_n.j() += 1;
        FVStencil coef_faces;
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_s, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_n, 0) + epsilon.GetDataVector().At(ind, 0));
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_s, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_n, 0) + rho.GetDataVector().At(ind, 0));
        SC beta_w = 0.5 * (beta.GetDataVector().At(ind_w, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_e = 0.5 * (beta.GetDataVector().At(ind_e, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_s = 0.5 * (beta.GetDataVector().At(ind_s, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_n = 0.5 * (beta.GetDataVector().At(ind_n, 0) + beta.GetDataVector().At(ind, 0));

        coef_faces(CNB::WEST, 0) = -dt * eps_w / (1. + beta_w * dt / (eps_w * rho_w));
        coef_faces(CNB::EAST, 0) = -dt * eps_e / (1. + beta_e * dt / (eps_e * rho_e));
        coef_faces(CNB::SOUTH, 0) = -dt * eps_s / (1. + beta_s * dt / (eps_s * rho_s));
        coef_faces(CNB::NORTH, 0) = -dt * eps_n / (1. + beta_n * dt / (eps_n * rho_n));

        SC c_ex{0.};
        for (auto face : g_s.GetFaces()) {
            coef_faces(face, 0) *= dA[dare::MapCartesianFaceToDim(face)] * dn_r[dare::MapCartesianFaceToDim(face)];
            c_ex += -coef_faces(face, 0);
        }
        SC dd = dd_val.GetDataVector().At(ind, 0);
        c_ex += dd * epsilon.GetDataVector().At(ind, 0) * dV / dt;

        EXPECT_NEAR(J.GetValue(CNB::WEST, 0), coef_faces(CNB::WEST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::WEST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::EAST, 0), coef_faces(CNB::EAST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::EAST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::SOUTH, 0), coef_faces(CNB::SOUTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::SOUTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::NORTH, 0), coef_faces(CNB::NORTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::NORTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::CENTER, 0), c_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(c_ex));
    }
}

TEST_F(ProjectionMethodCartesian3DTest, ContinuityJacobianCompressible) {
    struct density_derivative_functor {
        density_derivative_functor() : field(nullptr) {}
        explicit density_derivative_functor(const Field& dd) : field(&dd) {}
        SC operator()(Index ind) const { return field->GetDataVector().At(ind, 0); }
        const Field* field;
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
    using FVStencil = dare::FaceValueStencil<GridType, SC, 1>;
    using CNB = dare::CartesianNeighbor;

    auto g_s = grid->GetRepresentation(opt_s);
    double mu = 0.;
    const double tol_eps = 1e5;
    Field dd_val("dd", g_s, 1);  // value of density derivative
    Field rho("rho", g_s, 2);
    Field epsilon("epsilon", g_s, 2);
    Field beta("beta", g_s, 1);
    dare::PseudoRandomTGenerator<SC> rd(-1000, 1000);
    rd.SetPreFactor(1e-4);

    dare::ProjectionMethod<GridType, BStrat, PDict, NDict> pm;
    static_assert(pm.IsCompressible(), "This is for the incompressible approach!");

    dare::ConstantTimeStep dt(1.);
    dare::test::BStrat bstrat;
    pm.Initialize(grid, &dt, bstrat, bstrat, bstrat, bstrat);
    pm.SetDensity(&rho);
    pm.SetDensityDerivative(density_derivative_functor{dd_val});
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
        dd_val.GetDataVector().At(i) = rd.Generate();
    }

    dare::Vector<Dim, const GridVector*> velocities;
    for (std::size_t d{0}; d < Dim; d++) {
        velocities[d] = &pm.GetMomentum(d)->GetField()->GetDataVector(1);
    }

    dare::Vector<Dim, SC> dA = g_s.GetFaceArea();
    dare::Vector<Dim, SC> dn_r = 1. / g_s.GetDistances();
    SC dV = g_s.GetCellVolume();
    for (LO n_loc = 0; n_loc < g_s.GetNumberLocalCellsInternal(); n_loc++) {
        Index ind_loc = g_s.MapOrdinalToIndexLocalInternal(n_loc);
        Index ind = g_s.MapInternalToLocal(ind_loc);
        auto J = free_pm_continuity_Jacobian_Cartesian(&pm, n_loc, ind, pm.GetDensity(), pm.GetPorosity());
        static_assert(std::is_same_v<decltype(J), dare::CenterMatrixStencil<GridType, SC, 1>>,
                      "Should be a center matrix stencil");

        Index ind_w(ind), ind_e(ind), ind_s(ind), ind_n(ind), ind_b(ind), ind_t(ind);
        ind_w.i() -= 1;
        ind_e.i() += 1;
        ind_s.j() -= 1;
        ind_n.j() += 1;
        ind_b.k() -= 1;
        ind_t.k() += 1;

        FVStencil coef_faces;
        SC eps_w = 0.5 * (epsilon.GetDataVector().At(ind_w, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_e = 0.5 * (epsilon.GetDataVector().At(ind_e, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_s = 0.5 * (epsilon.GetDataVector().At(ind_s, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_n = 0.5 * (epsilon.GetDataVector().At(ind_n, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_b = 0.5 * (epsilon.GetDataVector().At(ind_b, 0) + epsilon.GetDataVector().At(ind, 0));
        SC eps_t = 0.5 * (epsilon.GetDataVector().At(ind_t, 0) + epsilon.GetDataVector().At(ind, 0));
        SC rho_w = 0.5 * (rho.GetDataVector().At(ind_w, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_e = 0.5 * (rho.GetDataVector().At(ind_e, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_s = 0.5 * (rho.GetDataVector().At(ind_s, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_n = 0.5 * (rho.GetDataVector().At(ind_n, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_b = 0.5 * (rho.GetDataVector().At(ind_b, 0) + rho.GetDataVector().At(ind, 0));
        SC rho_t = 0.5 * (rho.GetDataVector().At(ind_t, 0) + rho.GetDataVector().At(ind, 0));
        SC beta_w = 0.5 * (beta.GetDataVector().At(ind_w, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_e = 0.5 * (beta.GetDataVector().At(ind_e, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_s = 0.5 * (beta.GetDataVector().At(ind_s, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_n = 0.5 * (beta.GetDataVector().At(ind_n, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_b = 0.5 * (beta.GetDataVector().At(ind_b, 0) + beta.GetDataVector().At(ind, 0));
        SC beta_t = 0.5 * (beta.GetDataVector().At(ind_t, 0) + beta.GetDataVector().At(ind, 0));

        coef_faces(CNB::WEST, 0) = -dt * eps_w / (1. + beta_w * dt / (eps_w * rho_w));
        coef_faces(CNB::EAST, 0) = -dt * eps_e / (1. + beta_e * dt / (eps_e * rho_e));
        coef_faces(CNB::SOUTH, 0) = -dt * eps_s / (1. + beta_s * dt / (eps_s * rho_s));
        coef_faces(CNB::NORTH, 0) = -dt * eps_n / (1. + beta_n * dt / (eps_n * rho_n));
        coef_faces(CNB::BOTTOM, 0) = -dt * eps_b / (1. + beta_b * dt / (eps_b * rho_b));
        coef_faces(CNB::TOP, 0) = -dt * eps_t / (1. + beta_t * dt / (eps_t * rho_t));

        SC c_ex{0.};
        for (auto face : g_s.GetFaces()) {
            coef_faces(face, 0) *= dA[dare::MapCartesianFaceToDim(face)] * dn_r[dare::MapCartesianFaceToDim(face)];
            c_ex += -coef_faces(face, 0);
        }
        SC dd = dd_val.GetDataVector().At(ind, 0);
        c_ex += dd * epsilon.GetDataVector().At(ind, 0) * dV / dt;

        EXPECT_NEAR(J.GetValue(CNB::WEST, 0), coef_faces(CNB::WEST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::WEST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::EAST, 0), coef_faces(CNB::EAST, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::EAST, 0)));
        EXPECT_NEAR(J.GetValue(CNB::SOUTH, 0), coef_faces(CNB::SOUTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::SOUTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::NORTH, 0), coef_faces(CNB::NORTH, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::NORTH, 0)));
        EXPECT_NEAR(J.GetValue(CNB::BOTTOM, 0), coef_faces(CNB::BOTTOM, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::BOTTOM, 0)));
        EXPECT_NEAR(J.GetValue(CNB::TOP, 0), coef_faces(CNB::TOP, 0),
                    tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(coef_faces(CNB::TOP, 0)));
        EXPECT_NEAR(J.GetValue(CNB::CENTER, 0), c_ex, tol_eps * std::numeric_limits<SC>::epsilon() * std::abs(c_ex));
    }
}
