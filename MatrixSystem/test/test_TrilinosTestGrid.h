/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder
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

#ifndef MATRIXSYSTEM_TEST_TEST_TRILINOSTESTGRID_H_
#define MATRIXSYSTEM_TEST_TEST_TRILINOSTESTGRID_H_

#include "../../Grid/DefaultTypes.h"
#include "../../Utilities/Vector.h"
#include "../../MPI/ExecutionManager.h"

namespace dare::Matrix::test {
class TrilinosTestGrid {
public:
    using GlobalOrdinalType = dare::defaults::GlobalOrdinalType;
    using LocalOrdinalType = dare::defaults::LocalOrdinalType;
    using ScalarType = dare::defaults::ScalarType;
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
        size_global = 21 + exman->GetNumberProcesses() * 2;
        local_size = size_global / exman->GetNumberProcesses();
        LocalOrdinalType subdomain_size = local_size;
        if (exman->GetRank() == (exman->GetNumberProcesses() - 1))
            local_size += size_global - _exman->GetNumberProcesses() * local_size;
        offset = subdomain_size * exman->GetRank();
    }

    GlobalOrdinalType size_global{0};
    LocalOrdinalType local_size{0};
    GlobalOrdinalType offset{0};
    dare::mpi::ExecutionManager* exman;
};

static const std::size_t N{4};

}  // end namespace dare::Matrix::test

#endif  // MATRIXSYSTEM_TEST_TEST_TRILINOSTESTGRID_H_
