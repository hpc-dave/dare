/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder
 *
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
#include "../Trilinos.h"
#include "../../Grid/DefaultTypes.h"

namespace dare::Matrix::test {
class TrilinosTestGrid {
public:
    using GlobalOrdinalType = dare::Grid::details::GlobalOrdinalType;
    using LocalOrdinalType = dare::Grid::details::LocalOrdinalType;
    using ScalarType = double;
    using Index = dare::utils::Vector<1, LocalOrdinalType>;
    using IndexGlobal = dare::utils::Vector<1, GlobalOrdinalType>;
    class TestRepresentation {
    public:
        using GlobalOrdinalType = TrilinosTestGrid::GlobalOrdinalType;
        using LocalOrdinalType = TrilinosTestGrid::LocalOrdinalType;

        TestRepresentation() : TestRepresentation(nullptr) {}

        explicit TestRepresentation(dare::Matrix::test::TrilinosTestGrid* _grid) {
            grid = _grid;
        }

        GlobalOrdinalType MapLocalToGlobalInternal(LocalOrdinalType node) const {
            return node + grid->offset;
        }

        LocalOrdinalType MapGlobalToLocalInternal(GlobalOrdinalType node) const {
            return node - grid->offset;
        }

        LocalOrdinalType MapInternalToLocal(LocalOrdinalType n_internal) const {
            return n_internal;
        }

        bool IsLocalInternal(GlobalOrdinalType id_glob) const {
            return id_glob >= grid->offset && id_glob < (grid->offset * grid->local_size);
        }

        LocalOrdinalType GetNumberLocalCellsInternal() const {
            return grid->local_size;
        }

        LocalOrdinalType GetNumberLocalCells() const {
            return GetNumberLocalCellsInternal();
        }

        GlobalOrdinalType GetNumberGlobalCellsInternal() const {
            return grid->size_global;
        }

        dare::Matrix::test::TrilinosTestGrid* grid;
    };
    using Representation = TestRepresentation;

    Representation GetRepresentation() { return Representation(this); }

    void Initialize(dare::mpi::ExecutionManager* _exman) {
        exman = _exman;
        local_size = size_global / exman->GetNumberProcesses();
        if (exman->GetRank() == (exman->GetNumberProcesses() - 1))
            local_size += size_global - _exman->GetNumberProcesses() * local_size;
        offset = local_size * exman->GetRank();
    }

    GlobalOrdinalType size_global{10000};
    LocalOrdinalType local_size{0};
    GlobalOrdinalType offset{0};
    dare::mpi::ExecutionManager* exman;
};

static const std::size_t N{4};

}  // end namespace dare::Matrix::test

/*!
 * @brief Fixture for testing the trilinos implementation
 */
class TrilinosTest : public testing::Test {
public:
    using GridType = dare::Matrix::test::TrilinosTestGrid;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using GridRepresentation = typename GridType::Representation;
    static const std::size_t N = dare::Matrix::test::N;
    using FieldType = dare::Data::GridVector<GridType, SC, N>;
    GridType grid;
    FieldType field;
    dare::mpi::ExecutionManager exec_man;

    void SetUp() {
        grid.Initialize(&exec_man);
        field = FieldType("test field", grid.GetRepresentation());
    }
};

TEST_F(TrilinosTest, Initialize) {
    dare::Matrix::Trilinos<SC, LO, GO> trilinos(&exec_man);
    EXPECT_TRUE(trilinos.IsInitialized());
}

TEST_F(TrilinosTest, Build) {
    SC offset_rhs{1.5};
    SC offset_init{0.75};
    SC offset_coef{0.3};
    GridRepresentation g_rep{grid.GetRepresentation()};
    for (LO n{0}; n < grid.local_size; n++) {
        LO row = n * N;
        for (LO i{0}; i < N; i++) {
            field.At(n, i) = row + i;
        }
    }

    auto functor = KOKKOS_LAMBDA(auto mblock) {
        const std::size_t num_rows = grid.size_global * N;
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_TRUE(mblock->GetInitialGuess(n) == field.At(mblock->GetNode(), n));
            if (mblock->GetRow(n) > 0)
                mblock->SetCoefficient(n, mblock->GetRow(n) - 1, mblock->GetInitialGuess(n));
            mblock->SetCoefficient(n, mblock->GetRow(n), mblock->GetInitialGuess(n) + offset_coef);
            if (mblock->GetRow(n) < (num_rows - 1))
                mblock->SetCoefficient(n, mblock->GetRow(n) + 1, mblock->GetInitialGuess(n) + 2. * offset_coef);
            mblock->SetRhs(n, mblock->GetRow(n) + offset_rhs);
            mblock->SetInitialGuess(n, mblock->GetRow(n) + offset_init);
        }
    };

    dare::Matrix::Trilinos<SC, LO, GO> trilinos(&exec_man);
    trilinos.Build(g_rep, field, functor, false);

}
