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

#include <float.h>
#include <gtest/gtest.h>
#include "TestField.h"
#include "MPI/SingleHaloBuffer.h"

TEST(SingleHaloBufferTest, Initialize) {
    using SC = double;
    using BufferType = dare::mpi::SingleHaloBuffer<SC>;
    int test_rank{0};
    dare::mpi::ExecutionManager exman;
    BufferType buffer;
    buffer.Initialize(&exman, 0);
    ASSERT_EQ(buffer.GetPartnerRank(), test_rank);
}

TEST(SingleHaloBufferTest, CommunicateAmountHaloCellIDs) {
    using SC = double;
    using BufferType = dare::mpi::SingleHaloBuffer<SC>;
    dare::mpi::ExecutionManager exman;
    using BufferType = dare::mpi::SingleHaloBuffer<SC>;
    std::vector<BufferType> list_buffers(exman.GetNumberProcesses());
    for (std::size_t n{0}; n < list_buffers.size(); n++)
        list_buffers[n].Initialize(&exman, n);

    std::size_t send_buffer = exman.GetRank();
    std::vector<std::size_t> recv_buffer(exman.GetNumberProcesses());
    std::vector<MPI_Request> request(exman.GetNumberProcesses()*2);
    for (std::size_t n{0}; n < list_buffers.size(); n++)
        list_buffers[n].CommunicateNumCellIDs(send_buffer, &recv_buffer[n],
                                              &request[2 * n], &request[2 * n + 1]);
    MPI_Waitall(request.size(), request.data(), MPI_STATUSES_IGNORE);
    for (int n{0}; n < exman.GetNumberProcesses(); n++) {
        ASSERT_EQ(n, recv_buffer[n]);
    }
}

TEST(SingleHaloBufferTest, ExchangeAndFill) {
    using SC = double;
    using BufferType = dare::mpi::SingleHaloBuffer<SC>;
    using LO = typename BufferType::LO;
    std::size_t num_halo_IDs = 10;
    dare::mpi::ExecutionManager exman;
    dare::test::TestField field;
    // field needs for each process a halo cell region and local data to send
    field.ResizeByGridSize(num_halo_IDs * (exman.GetNumberProcesses() + 1));
    std::size_t buffer_size = num_halo_IDs;

    std::vector<BufferType> list_buffers(exman.GetNumberProcesses());
    for (std::size_t n{0}; n < list_buffers.size(); n++)
        list_buffers[n].Initialize(&exman, n);

    // we send the last block of IDs in the field to the other processes
    std::vector<LO> send_IDs(buffer_size);
    for (std::size_t n{0}; n < buffer_size; n++) {
        send_IDs[n] = exman.GetNumberProcesses() * num_halo_IDs + n;
    }
    for (std::size_t m{0}; m < list_buffers.size(); m++) {
        std::vector<LO> recv_IDs(buffer_size);
        for (std::size_t n{0}; n < buffer_size; n++) {
            recv_IDs[n] = m * num_halo_IDs + n;
        }
        list_buffers[m].FinalizeInitialization(send_IDs, recv_IDs);
    }
    // Initialize the field values to a distinguishable value
    std::size_t size_field = (num_halo_IDs * (exman.GetNumberProcesses() + 1) * field.GetNumEquations());
    std::size_t size_send = (num_halo_IDs * exman.GetNumberProcesses() * field.GetNumEquations());
    for (std::size_t n{0}; n < size_send; n++)
        field.At(n) = -1.;
    for (std::size_t n{size_send}; n < size_field; n++)
        field.At(n) = exman.GetRank() + n - size_send;

    std::vector<MPI_Request> requests(exman.GetNumberProcesses() * 2);
    for (std::size_t m{0}; m < list_buffers.size(); m++) {
        list_buffers[m].Exchange(field, &requests[m * 2], &requests[m * 2 + 1]);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    for (std::size_t m{0}; m < list_buffers.size(); m++) {
        list_buffers[m].FillValues(&field);
    }
    for (std::size_t n{0}; n < list_buffers.size(); n++) {
        for (std::size_t m{0}; m < buffer_size * field.GetNumEquations(); m++) {
            double value_field = field.At(n * buffer_size * field.GetNumEquations() + m);
            double value_expec = n + m;
            bool success = std::abs(value_field - value_expec) < DBL_EPSILON;
            EXPECT_TRUE(success) << "communication with field was not successful, expect "
                                 << value_expec << " but found " << value_field;
        }
    }
}
