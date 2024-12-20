/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder

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
#include <float.h>
#include <unordered_map>

#include "MPI/ExecutionManager.h"
#include "MPI/HaloBuffer.h"
#include "TestField.h"

TEST(HaloBufferTest, Exchange) {
    using SC = double;
    using LO = typename dare::mpi::HaloBuffer<SC>::LO;
    using GO = typename dare::mpi::HaloBuffer<SC>::GO;

    dare::mpi::ExecutionManager exman;
    dare::mpi::HaloBuffer<SC> buffer;
    dare::test::TestField field;
    std::size_t cell_per_domain = 2;
    LO num_ghost = 1;
    // std::size_t global_field_size = exman.GetNumberProcesses() * cell_per_domain + 2 * num_ghost;
    field.ResizeByGridSize(cell_per_domain + 2*num_ghost);
    std::vector<GO> list_required_IDs(2 * num_ghost);
    std::unordered_map<GO, GO> map_periodic;
    for (int n{0}; n < num_ghost; n++) {
        // cells on lower end
        list_required_IDs[n] = exman.GetRank() * cell_per_domain + n;
        // cells on upper end
        list_required_IDs[n + num_ghost] = (exman.GetRank() + 1) * cell_per_domain + num_ghost + n;
        // map periodic IDs
        if (exman.GetRank() == 0 || exman.GetRank() == (exman.GetNumberProcesses() - 1)) {
            map_periodic[n] = exman.GetNumberProcesses() * cell_per_domain + n;
            map_periodic[exman.GetNumberProcesses() * cell_per_domain + n] = n;
            map_periodic[exman.GetNumberProcesses() * cell_per_domain + num_ghost + n] = num_ghost + n;
            map_periodic[num_ghost + n] = exman.GetNumberProcesses() * cell_per_domain + num_ghost + n;
        }
    }
    GO ID_low = exman.GetRank() * cell_per_domain;
    GO ID_high = ID_low + cell_per_domain + 2 * num_ghost;
    auto is_local = [&](GO id) {
        return id >= ID_low && id < ID_high;
    };
    auto is_internal = [&](GO id) {
        return id >= (ID_low+num_ghost) && id < (ID_high - num_ghost);
    };
    auto map_global_local = [&](GO id) {
        return id - ID_low;
    };
    buffer.Initialize(&exman, list_required_IDs, map_periodic, is_local, is_internal, map_global_local);

    // Initializing field with default value
    for (std::size_t n{0}; n < (cell_per_domain + 2 * num_ghost); n++) {
        for (std::size_t e{0}; e < field.GetNumEquations(); e++) {
            field.At(n * field.GetNumEquations() + e) = -1.;
        }
    }
    // Filling internal field with values;
    for (std::size_t n{static_cast<std::size_t>(num_ghost)}; n < (num_ghost + cell_per_domain); n++) {
        GO v_start = static_cast<GO>((exman.GetRank() * cell_per_domain + n) * field.GetNumEquations());
        for (std::size_t e{0}; e < field.GetNumEquations(); e++)
            field.At(n * field.GetNumEquations() + e) = v_start + e;
    }

    // Exchanging the values
    buffer.Exchange(&field);

    for (LO n{0}; n < num_ghost; n++) {
        LO id_low{static_cast<LO>(n * field.GetNumEquations())};
        LO id_high{ static_cast<LO>((num_ghost + cell_per_domain) * field.GetNumEquations())};
        GO v_start_low = (exman.GetRank() * cell_per_domain + n) * field.GetNumEquations();
        GO v_start_high = ((exman.GetRank() + 1) * cell_per_domain + num_ghost + n) * field.GetNumEquations();
        if (exman.GetRank() == 0) {
            v_start_low = (exman.GetNumberProcesses() * cell_per_domain + n) * field.GetNumEquations();
        }
        if (exman.GetRank() == (exman.GetNumberProcesses() - 1)) {
            v_start_high = num_ghost * field.GetNumEquations();
        }
        for (std::size_t e{0}; e < field.GetNumEquations(); e++) {
            double value_low_f = field.At(id_low + e);
            double value_high_f = field.At(id_high + e);
            double value_low_e = v_start_low + e;
            double value_high_e = v_start_high + e;
            double err_low = std::abs(value_low_f - value_low_e);
            double err_high = std::abs(value_high_f - value_high_e);
            EXPECT_TRUE(err_low < DBL_EPSILON)
                << "Rank " << exman.GetRank() << ": lower value wrong! "
                << "Expected " << value_low_e << " but found " << value_low_f << "!";
            EXPECT_TRUE(err_high < DBL_EPSILON)
                << "Rank " << exman.GetRank() << ": upper value wrong! "
                << "Expected " << value_high_e << " but found " << value_high_f << "!";
        }
    }
}
