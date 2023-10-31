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

namespace dare::mpi {

template <typename LO, typename GO, typename SC>
HaloBuffer<LO, GO, SC>::HaloBuffer() : exec_man(nullptr) {}

template <typename LO, typename GO, typename SC>
template <typename Lambda1, typename Lambda2>
void HaloBuffer<LO, GO, SC>::Initialize(ExecutionManager* execution_manager,
                                    const std::vector<GO>& list_required_IDs,
                                    Lambda1 is_local,
                                    Lambda2 map_global_to_local) {
    exec_man = execution_manager;
    for (int n{0}; n < exec_man->GetNumberProcesses(); n++) {
        buffers[n] = SingleHaloBuffer<LO, GO, SC>(exec_man, n);
    }

    std::vector<MPI_Request> requests(exec_man->GetNumberProcesses() * 2);
    std::vector<std::size_t> num_send(exec_man->GetNumberProcesses(), list_required_IDs.size());
    std::vector<std::size_t> num_recv(exec_man->GetNumberProcesses(), 0);
    for (int n{0}; n < exec_man->GetNumberProcesses(); n++) {
        buffers[n].CommunicateAmountRequiredHaloCellIDs(num_send[n], &num_recv[n],
                                                        &requests[2 * n], &requests[2 * n + 1]);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    std::vector<std::vector<GO>> list_required_IDs_partner(exec_man->GetNumberProcesses());
    for (int n{0}; n < exec_man->GetNumberProcesses(); n++) {
        list_required_IDs_partner[n].resize(num_recv[n]);
        buffers[n].CommunicateRequiredHaloCellIDs(list_required_IDs, &list_required_IDs_partner[n],
                                                  &requests[2 * n], &requests[2 * n + 1]);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    std::vector<std::vector<GO>> list_filtered_IDs(exec_man->GetNumberProcesses());
    for (int n{0}; n < exec_man->GetNumberProcesses(); n++) {
        for (std::size_t m{0}; m < list_required_IDs_partner[n].size(); n++) {
            GO id = list_required_IDs_partner[n][m];
            if (is_local(id)) {
                list_filtered_IDs[n].push_back(id);
            }
        }
    }

    for (int n{0}; n < exec_man->GetNumberProcesses(); n++)
        num_send[n] = list_filtered_IDs[n].size();
    for (int n{0}; n < exec_man->GetNumberProcesses(); n++) {
        buffers[n].CommunicateAmountRequiredHaloCellIDs(num_send[n], &num_recv[n],
                                                        &requests[2 * n], &requests[2 * n + 1]);
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    for (int n{exec_man->GetNumberProcesses() - 1}; n >= 0; n--) {
        if (num_send[n] == 0) {
#ifndef DARE_NDEBUG
            if (num_recv[n] != 0) {
                std::ostringstream os;
                os << "Rank " << exec_man->GetRank() << ": partner process " << n << " requires communcation "
                   << "but this process does not have anything to send. The buffers cannot deal with that situation. "
                   << "Either debug or revise buffers!";
                std::cout << os.str() << std::endl;
            }
#endif
            buffers.erase(n);
        }
    }

    std::vector<std::vector<GO>> list_filtered_IDs_partner(exec_man->GetNumberProcesses());
    requests.resize(buffers.size() * 2);
    int count{0};
    for (auto& entry : buffers) {
        int proc{entry.first};
        list_filtered_IDs_partner[proc].resize(num_recv[proc]);
        entry.second.CommunicateFilteredHaloCellIDs(list_filtered_IDs[proc], &list_filtered_IDs_partner[proc],
                                                  &requests[2 * count], &requests[2 * count + 1]);
        count++;
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    std::vector<LO> list_IDs_recv, list_IDs_send;
    for (auto& entry : buffers) {
        int proc{entry.first};
        list_IDs_recv.resize(list_filtered_IDs_partner[proc].size());
        list_IDs_send.resize(list_filtered_IDs[proc].size());
        for (std::size_t n{0}; n < list_IDs_recv.size(); n++) {
            GO id_glob = list_filtered_IDs_partner[proc][n];
            LO id_loc = map_global_to_local(id_glob);
            list_IDs_recv[n] = id_loc;
        }
        for (std::size_t n{0}; n < list_IDs_send.size(); n++) {
            GO id_glob = list_filtered_IDs[proc][n];
            LO id_loc = map_global_to_local(id_glob);
            list_IDs_send[n] = id_loc;
        }
        entry.second.FinalizeInitialization(list_IDs_send, list_IDs_recv);
    }

    this->utils::InitializationTracker::Initialize();
}

template <typename LO, typename GO, typename SC>
template <typename Field>
void HaloBuffer<LO, GO, SC>::Exchange(Field* field) {
    std::vector<MPI_Request> requests(buffers.size() * 2);
    std::size_t count{0};
    for (auto& entry : buffers) {
        entry.second.Exchange(*field, &requests[count * 2], &requests[count * 2 + 1]);
        count++;
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    for (auto& entry : buffers) {
        entry.second.FillValues(field);
    }
}

}  // namespace dare::mpi
