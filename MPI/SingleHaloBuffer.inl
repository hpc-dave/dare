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

namespace dare::mpi {

template<typename SC>
SingleHaloBuffer<SC>::SingleHaloBuffer(ExecutionManager* execution_manager,
                                   int rank_partner)
                                   : exec_man(execution_manager),
                                     rank_partner_proc(rank_partner) {}

template<typename SC>
void SingleHaloBuffer<SC>::Initialize(ExecutionManager* execution_manager,
                                          int rank_partner) {
    exec_man = execution_manager;
    rank_partner_proc = rank_partner;
}

template<typename SC>
void SingleHaloBuffer<SC>::CommunicateNumCellIDs(std::size_t num_send,
                                                         std::size_t* num_recv,
                                                         MPI_Request* request_send,
                                                         MPI_Request* request_recv) {
    const int tag_communicated{2001};

    exec_man->Iexchange(&num_send, 1, num_recv, 1, rank_partner_proc, tag_communicated, request_send, request_recv);
}

template<typename SC>
void SingleHaloBuffer<SC>::CommunicateRequiredHaloCellIDs(const std::vector<GO>& list_global_ID,
                                                              std::vector<GO>* list_global_IDs_partner,
                                                              MPI_Request* request_send,
                                                              MPI_Request* request_receive) {
    const int tag{2002};

    exec_man->Iexchange(list_global_ID.data(), list_global_ID.size(), list_global_IDs_partner->data(),
                        list_global_IDs_partner->size(), rank_partner_proc, tag,
                        request_send, request_receive);
}

template<typename SC>
void SingleHaloBuffer<SC>::CommunicateFilteredHaloCellIDs(const std::vector<GO>& list_filtered_IDs,
                                                              std::vector<GO>* list_filtered_IDs_partner,
                                                              MPI_Request* request_send, MPI_Request* request_receive) {
    const int tag{2003};

#ifndef DARE_NDEBUG
    if (list_filtered_IDs.size() == 0 || list_filtered_IDs_partner->size() == 0) {
        std::ostringstream os;
        os << std::string("Rank ");
        os   << exec_man->GetRank()
        << ": One of the buffers of sending and receiving IDs has the size 0!";
        exec_man->Terminate(__func__, os.str());
    }
#endif

    exec_man->Iexchange(list_filtered_IDs.data(), list_filtered_IDs.size(),
                        list_filtered_IDs_partner->data(), list_filtered_IDs_partner->size(),
                        rank_partner_proc, tag, request_send, request_receive);
}

template<typename SC>
void SingleHaloBuffer<SC>::FinalizeInitialization(const std::vector<LO>& list_IDs_send,
                                                          const std::vector<LO>& list_IDs_recv) {
    // the send and receive buffers
    list_local_IDs_send = list_IDs_send;
    list_local_IDs_recv = list_IDs_recv;

    this->utils::InitializationTracker::Initialize();
}

template<typename SC>
template <typename Field>
void SingleHaloBuffer<SC>::Exchange(const Field& field,
                                           MPI_Request* request_send, MPI_Request* request_recv) {
    std::size_t num_eq = field.GetNumEquations();
    const int tag_exchange{2100};
    int num_send = list_local_IDs_send.size() * num_eq;
    int num_recv = list_local_IDs_recv.size() * num_eq;

    if (num_send > buffer_send.size())
        buffer_send.resize(num_send);
    if (num_recv > buffer_recv.size())
        buffer_recv.resize(num_send);

    for (std::size_t n{0}; n < list_local_IDs_send.size(); n++) {
        LO id = list_local_IDs_send[n] * num_eq;
        for (std::size_t e{0}; e < num_eq; e++)
            buffer_send[n * num_eq + e] = field.At(id + e);
    }

    exec_man->Iexchange(buffer_send.data(), num_send, buffer_recv.data(), num_recv,
                               rank_partner_proc, tag_exchange, request_send, request_recv);
}

template<typename SC>
template <typename Field>
void SingleHaloBuffer<SC>::FillValues(Field* field) {
    std::size_t num_eq = field->GetNumEquations();

    for (std::size_t n{0}; n < list_local_IDs_recv.size(); n++) {
        LO id = list_local_IDs_recv[n] * num_eq;
        for (std::size_t e{0}; e < num_eq; e++)
            field->At(id + e) = buffer_recv[n * num_eq + e];
    }
}

template<typename SC>
int SingleHaloBuffer<SC>::GetPartnerRank() const {
    return rank_partner_proc;
}
}  // namespace dare::mpi
