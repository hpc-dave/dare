/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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

#include <memory>

#include "../Data/Field.h"
#include "../Grid/Cartesian.h"

namespace dare::test {

template <std::size_t Dim, typename GO>
dare::utils::Vector<Dim, GO> GetResolutionIntegrationTestCartesianField() {
    dare::utils::Vector<Dim, GO> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 20 + n;
    return res;
}
template <std::size_t Dim, typename SC>
dare::utils::Vector<Dim, SC> GetSizeIntegrationTestCartesianField() {
    dare::utils::Vector<Dim, SC> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

}  // namespace dare::test

template <std::size_t Dimension>
class IntegrationCartesianField
    : public testing::TestWithParam<typename dare::Grid::Cartesian<Dimension>::Options> {
public:
    static const std::size_t N{3};
    static const std::size_t num_time_steps{3};
    using GridType = dare::Grid::Cartesian<Dimension>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using SC = double;
    using Field = dare::Data::Field<GridType, SC, N>;
    using GridVector = typename Field::VectorType;

    void SetUp() {
        const LO num_ghost{2};
        grid = std::make_unique<GridType>(
            &exec_man,
            dare::test::GetResolutionIntegrationTestCartesianField<Dimension, GO>(),
            dare::test::GetSizeIntegrationTestCartesianField<Dimension, SC>(),
            num_ghost);
    }

    std::unique_ptr<GridType> grid;
    dare::mpi::ExecutionManager exec_man;
};

using IntegrationCartesianField1D = IntegrationCartesianField<1>;
using IntegrationCartesianField2D = IntegrationCartesianField<2>;
using IntegrationCartesianField3D = IntegrationCartesianField<3>;

TEST_F(IntegrationCartesianField1D, Initialize) {
    typename GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    Field field("test_1D", g_rep, num_time_steps);
    GridVector* data = nullptr;
    for (std::size_t n{0}; n < num_time_steps; n++) {
        data = &field.GetDataVector(n);
        EXPECT_EQ(data->GetSize(), g_rep.GetNumberLocalCells() * N);
        EXPECT_EQ(data->GetIdentifier(), "test_1D_" + std::to_string(n));
    }

    typename GridType::Options opt_s(1);  // staggered
    auto g_rep_staggered = grid->GetRepresentation(opt_s);
    Field field_s("test_1D_staggered", g_rep_staggered, num_time_steps);
    for (std::size_t n{0}; n < num_time_steps; n++) {
        data = &field_s.GetDataVector(n);
        EXPECT_EQ(data->GetSize(), g_rep_staggered.GetNumberLocalCells() * N);
        EXPECT_EQ(data->GetIdentifier(), "test_1D_staggered_" + std::to_string(n));
    }
}

TEST_F(IntegrationCartesianField2D, Initialize) {
    typename GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    Field field("test_2D", g_rep, num_time_steps);
    GridVector* data = nullptr;
    for (std::size_t n{0}; n < num_time_steps; n++) {
        data = &field.GetDataVector(n);
        EXPECT_EQ(data->GetSize(), g_rep.GetNumberLocalCells() * N);
        EXPECT_EQ(data->GetIdentifier(), "test_2D_" + std::to_string(n));
    }

    typename GridType::Options opt_s(1);  // staggered
    auto g_rep_staggered = grid->GetRepresentation(opt_s);
    Field field_s("test_2D_staggered", g_rep_staggered, num_time_steps);
    for (std::size_t n{0}; n < num_time_steps; n++) {
        data = &field_s.GetDataVector(n);
        EXPECT_EQ(data->GetSize(), g_rep_staggered.GetNumberLocalCells() * N);
        EXPECT_EQ(data->GetIdentifier(), "test_2D_staggered_" + std::to_string(n));
    }
}

TEST_F(IntegrationCartesianField3D, Initialize) {
    typename GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    Field field("test_3D", g_rep, num_time_steps);
    GridVector* data = nullptr;
    for (std::size_t n{0}; n < num_time_steps; n++) {
        data = &field.GetDataVector(n);
        EXPECT_EQ(data->GetSize(), g_rep.GetNumberLocalCells() * N);
        EXPECT_EQ(data->GetIdentifier(), "test_3D_" + std::to_string(n));
    }

    typename GridType::Options opt_s(1);  // staggered
    auto g_rep_staggered = grid->GetRepresentation(opt_s);
    Field field_s("test_3D_staggered", g_rep_staggered, num_time_steps);
    for (std::size_t n{0}; n < num_time_steps; n++) {
        data = &field_s.GetDataVector(n);
        EXPECT_EQ(data->GetSize(), g_rep_staggered.GetNumberLocalCells() * N);
        EXPECT_EQ(data->GetIdentifier(), "test_3D_staggered_" + std::to_string(n));
    }
}

TEST_F(IntegrationCartesianField1D, CopyToOld) {
    typename GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    Field field("test_1D", g_rep, num_time_steps);
    GridVector* data = nullptr;
    for (std::size_t tstep{1}; tstep < num_time_steps; tstep++) {
        data = &field.GetDataVector(tstep);
        for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
            Index ind = g_rep.MapOrdinalToIndexLocal(node);
            for (std::size_t n{0}; n < N; n++) {
                data->At(ind, n) = 0.;
            }
        }
    }
    data = &field.GetDataVector();
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            data->At(ind, n) = node * N + n;
        }
    }
    field.CopyDataVectorsToOldTimeStep();
    data = &field.GetDataVector(1);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(ind, n), node * N + n);
        }
    }
    data = &field.GetDataVector(2);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(node, n), 0.);
        }
    }
    field.CopyDataVectorsToOldTimeStep();
    data = &field.GetDataVector(2);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(ind, n), node * N + n);
        }
    }
}

TEST_F(IntegrationCartesianField2D, CopyToOld) {
    typename GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    Field field("test_2D", g_rep, num_time_steps);
    GridVector* data = nullptr;
    for (std::size_t tstep{1}; tstep < num_time_steps; tstep++) {
        data = &field.GetDataVector(tstep);
        for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
            Index ind = g_rep.MapOrdinalToIndexLocal(node);
            for (std::size_t n{0}; n < N; n++) {
                data->At(ind, n) = 0.;
            }
        }
    }
    data = &field.GetDataVector();
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            data->At(ind, n) = node * N + n;
        }
    }
    field.CopyDataVectorsToOldTimeStep();
    data = &field.GetDataVector(1);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(ind, n), node * N + n);
        }
    }
    data = &field.GetDataVector(2);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(node, n), 0.);
        }
    }
    field.CopyDataVectorsToOldTimeStep();
    data = &field.GetDataVector(2);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(ind, n), node * N + n);
        }
    }
}

TEST_F(IntegrationCartesianField3D, CopyToOld) {
    typename GridType::Options opt(0);  // not staggered
    auto g_rep = grid->GetRepresentation(opt);
    Field field("test_3D", g_rep, num_time_steps);
    GridVector* data = nullptr;
    for (std::size_t tstep{1}; tstep < num_time_steps; tstep++) {
        data = &field.GetDataVector(tstep);
        for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
            Index ind = g_rep.MapOrdinalToIndexLocal(node);
            for (std::size_t n{0}; n < N; n++) {
                data->At(ind, n) = 0.;
            }
        }
    }
    data = &field.GetDataVector();
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            data->At(ind, n) = node * N + n;
        }
    }
    field.CopyDataVectorsToOldTimeStep();
    data = &field.GetDataVector(1);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(ind, n), node * N + n);
        }
    }
    data = &field.GetDataVector(2);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(node, n), 0.);
        }
    }
    field.CopyDataVectorsToOldTimeStep();
    data = &field.GetDataVector(2);
    for (std::size_t node{0}; node < g_rep.GetNumberLocalCells(); node++) {
        Index ind = g_rep.MapOrdinalToIndexLocal(node);
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(data->At(ind, n), node * N + n);
        }
    }
}

TEST_F(IntegrationCartesianField1D, Exchange) {
}