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
            if (grid) {
                offset = _grid->offset;
                local_size = _grid->local_size;
                size_global = _grid->size_global;
            }
        }

        GlobalOrdinalType MapLocalToGlobalInternal(LocalOrdinalType node) const {
            return node + offset;
        }

        LocalOrdinalType MapGlobalToLocalInternal(GlobalOrdinalType node) const {
            return node - offset;
        }

        LocalOrdinalType MapInternalToLocal(LocalOrdinalType n_internal) const {
            return n_internal;
        }

        GlobalOrdinalType MapInternalToLocal(GlobalOrdinalType n_internal) const {
            return n_internal;
        }

        bool IsLocalInternal(GlobalOrdinalType id_glob) const {
            return id_glob >= grid->offset && id_glob < (grid->offset * grid->local_size);
            // return id_glob >= offset && id_glob < (offset * local_size);
        }

        LocalOrdinalType GetNumberLocalCellsInternal() const {
            // return local_size;
            return grid->local_size;
        }

        LocalOrdinalType GetNumberLocalCells() const {
            return GetNumberLocalCellsInternal();
        }

        GlobalOrdinalType GetNumberGlobalCellsInternal() const {
            return grid->size_global;
            // return size_global;
        }

        dare::Matrix::test::TrilinosTestGrid* grid;
        GlobalOrdinalType size_global{0};
        LocalOrdinalType local_size{0};
        GlobalOrdinalType offset{0};
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
    using GOViewType = typename dare::Matrix::Trilinos<SC, LO, GO>::GOViewType;
    using LOViewType = typename dare::Matrix::Trilinos<SC, LO, GO>::LOViewType;
    using SViewType = typename dare::Matrix::Trilinos<SC, LO, GO>::SViewType;
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
    for (LO node = 0; node < grid.local_size; node++) {
        LO row = node * N;
        for (LO i{0}; i < N; i++) {
            field.At(node, i) = row + i;
        }
    }

    auto functor = [&](auto mblock) {
        const std::size_t num_rows = grid.size_global * N;
        GO node_g = mblock->GetNode();
        LO node_l = g_rep.MapInternalToLocal(g_rep.MapGlobalToLocalInternal(node_g));
        for (std::size_t n{0}; n < N; n++) {
            EXPECT_EQ(mblock->GetInitialGuess(n), field.At(node_l, n));
            EXPECT_TRUE(mblock->IsGlobal());
            if (mblock->GetRow(n) == 0 || mblock->GetRow(n) == (num_rows - 1)) {
                mblock->Resize(n, 2);
            } else {
                mblock->Resize(n, 3);
            }
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

    LO num_rows = grid.local_size * N;
    GO num_rows_g = grid.size_global * N;
    ASSERT_EQ(trilinos.GetMap()->getLocalNumElements(), num_rows);
    ASSERT_EQ(trilinos.GetX()->getLocalLength(), num_rows);
    ASSERT_EQ(trilinos.GetB()->getLocalLength(), num_rows);

    for (LO node{0}; node < grid.local_size; node++) {
        GO node_g = g_rep.MapLocalToGlobalInternal(node);
        for (LO i{0}; i < N; i++) {
            LO row = node * N + i;
            GO row_g = node_g * N + i;
            EXPECT_EQ(trilinos.GetX()->getData()[row], row_g + offset_init);
            EXPECT_EQ(trilinos.GetB()->getData()[row], row_g + offset_rhs);
            if (row == 0 || row == (grid.local_size * N - 1)) {
                std::size_t num_entries = trilinos.GetA()->getNumEntriesInGlobalRow(row_g);
                GOViewType indices("ind", num_entries);
                SViewType values("coef", num_entries);
                trilinos.GetA()->getGlobalRowCopy(row_g, indices, values, num_entries);
                std::vector<GO> ind_sorted(num_entries);
                std::vector<SC> val_sorted(num_entries);
                for (std::size_t n{0}; n < num_entries; n++) {
                    ind_sorted[n] = indices[n];
                    val_sorted[n] = values[n];
                }
                for (std::size_t n{0}; n < num_entries; n++) {
                    std::size_t lowest_id_pos = n;
                    for (std::size_t pos{n + 1}; pos < num_entries; pos++) {
                        if (ind_sorted[pos] < ind_sorted[lowest_id_pos])
                            lowest_id_pos = pos;
                    }
                    if (lowest_id_pos != n) {
                        std::swap(ind_sorted[n], ind_sorted[lowest_id_pos]);
                        std::swap(val_sorted[n], val_sorted[lowest_id_pos]);
                    }
                }
                if (node_g == 0) {
                    EXPECT_EQ(num_entries, 2);
                    EXPECT_EQ(ind_sorted.size(), 2);
                    EXPECT_EQ(val_sorted.size(), 2);
                    EXPECT_EQ(ind_sorted[0], 0);
                    EXPECT_EQ(val_sorted[0], field.At(0, 0) + offset_coef);
                    EXPECT_EQ(ind_sorted[1], 1);
                    EXPECT_EQ(val_sorted[1], field.At(0, 0) + 2. * offset_coef);
                } else if (node_g == (grid.size_global - 1)) {
                    EXPECT_EQ(num_entries, 2);
                    EXPECT_EQ(ind_sorted.size(), 2);
                    EXPECT_EQ(val_sorted.size(), 2);
                    EXPECT_EQ(ind_sorted[0], num_rows_g - 2);
                    EXPECT_EQ(val_sorted[0], field.At(grid.local_size - 1, N - 1));
                    EXPECT_EQ(ind_sorted[1], num_rows_g - 1);
                    EXPECT_EQ(val_sorted[1], field.At(grid.local_size - 1, N - 1) + offset_coef);
                } else {
                    EXPECT_EQ(num_entries, 3);
                    EXPECT_EQ(ind_sorted.size(), 3);
                    EXPECT_EQ(val_sorted.size(), 3);
                    EXPECT_EQ(ind_sorted[0], row_g - 1);
                    EXPECT_EQ(val_sorted[0], field.At(node, i));
                    EXPECT_EQ(ind_sorted[1], row_g);
                    EXPECT_EQ(val_sorted[1], field.At(node, i) + offset_coef);
                    EXPECT_EQ(ind_sorted[2], row_g + 1);
                    EXPECT_EQ(val_sorted[2], field.At(node, i) + 2. * offset_coef);
                }
            } else {
                std::size_t num_entries = trilinos.GetA()->getNumEntriesInLocalRow(row);
                LOViewType indices("ind", num_entries);
                SViewType values("coef", num_entries);
                trilinos.GetA()->getLocalRowCopy(row, indices, values, num_entries);
                std::vector<LO> ind_sorted(num_entries);
                std::vector<SC> val_sorted(num_entries);
                for (std::size_t n{0}; n < num_entries; n++) {
                    ind_sorted[n] = indices[n];
                    val_sorted[n] = values[n];
                }
                for (std::size_t n{0}; n < num_entries; n++) {
                    std::size_t lowest_id_pos = n;
                    for (std::size_t pos{n + 1}; pos < num_entries; pos++) {
                        if (ind_sorted[pos] < ind_sorted[lowest_id_pos])
                            lowest_id_pos = pos;
                    }
                    if (lowest_id_pos != n) {
                        std::swap(ind_sorted[n], ind_sorted[lowest_id_pos]);
                        std::swap(val_sorted[n], val_sorted[lowest_id_pos]);
                    }
                }
                EXPECT_EQ(num_entries, 3);
                EXPECT_EQ(ind_sorted.size(), 3);
                EXPECT_EQ(val_sorted.size(), 3);
                EXPECT_EQ(ind_sorted[0], row - 1);
                EXPECT_EQ(val_sorted[0], field.At(node, i));
                EXPECT_EQ(ind_sorted[1], row);
                EXPECT_EQ(val_sorted[1], field.At(node, i) + offset_coef);
                EXPECT_EQ(ind_sorted[2], row + 1);
                EXPECT_EQ(val_sorted[2], field.At(node, i) + 2. * offset_coef);
            }
        }
    }
}
