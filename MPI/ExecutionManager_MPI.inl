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

#ifndef MPI_EXECUTIONMANAGER_MPI_INL_
#define MPI_EXECUTIONMANAGER_MPI_INL_
namespace dare::mpi {

template <typename T>
void ExecutionManager::Exchange(const T* data, int count_send, T* recv,
                                int count_recv, int rank_other, int tag,
                                MPI_Status* status_send, MPI_Status* status_recv) {
    MPI_Request request;
    Isend(data, count_send, rank_other, tag, &request);
    Recv(recv, count_recv, rank_other, tag, status_recv);
    MPI_Wait(&request, status_send);
}

template <typename T>
void ExecutionManager::Exchange(const T* data, int count_send,
              std::vector<T>* recv,
              int rank_other, int tag,
              MPI_Status* status_send, MPI_Status* status_recv) {
    MPI_Request request;
    Isend(data, count_send, rank_other, tag, &request);
    Probe(rank_other, tag, status_recv);
    int count_recv{0};
    MPI_Get_count(status_recv, GetMPIType<T>(), count_recv);
    recv->resize(count_recv);
    Recv(recv->data(), count_recv, rank_other, tag, status_recv);
    MPI_Wait(&request, status_send);
}

template <typename T>
void ExecutionManager::Exchange(const std::vector<T>& send,
              std::vector<T>* recv,
              int rank_other, int tag,
              MPI_Status* status_send, MPI_Status* status_recv) {
    Exchange(send.data(), send.size(), recv, rank_other, tag, status_send, status_recv);
}

template <typename T>
void ExecutionManager::Iexchange(const T* data, int count_send,
               T* recv, int count_recv,
               int rank_other, int tag,
               MPI_Request* request_send, MPI_Request* request_recv) {
    Isend(data, count_send, rank_other, tag, request_send);
    Irecv(recv, count_recv, rank_other, tag, request_recv);
}

template <typename T>
T ExecutionManager::Allreduce(const T data, MPI_Op op) {
    T buf{static_cast<T>(0)};
    Allreduce(&data, &buf, 1, op);
    return buf;
}

template <typename T>
void ExecutionManager::Allreduce(const T* data, T* recv, int count, MPI_Op op) {
    MPI_Allreduce(data, recv, count, GetMPIType<T>(), op, communicator);
}

template <typename T>
T ExecutionManager::Allsum(const T data) {
    T buf{static_cast<T>(0)};
    Allsum(&data, &buf, 1);
    return buf;
}

template <typename T>
void ExecutionManager::Allsum(const T* data, T* recv, int count) {
    MPI_Allreduce(data, recv, count, GetMPIType<T>(), MPI_SUM, communicator);
}

template <typename T>
T ExecutionManager::Allmax(const T data) {
    T buf{static_cast<T>(0)};
    MPI_Allmax(&data, &buf, 1);
    return buf;
}

template <typename T>
void ExecutionManager::Allmax(const T* data, T* recv, int count) {
    MPI_Allreduce(data, recv, count, GetMPIType<T>(), MPI_MAX);
}

template <typename T>
T ExecutionManager::Allmin(T data) {
    T buf{static_cast<T>(0)};
    MPI_Allmmin(&data, &buf, 1);
    return buf;
}

template <typename T>
void ExecutionManager::Allmin(const T* data, T* recv, int count) {
    MPI_Allreduce(data, recv, count, GetMPIType<T>(), MPI_MIN);
}

inline bool ExecutionManager::AllLogicAnd(const bool data) {
    return Allreduce(data, MPI_LAND);
}

inline void ExecutionManager::AllLogicAnd(const bool* data, bool* recv, int count) {
    Allreduce(data, recv, count, MPI_LAND);
}

template <typename T>
inline int ExecutionManager::Allgather(const T* data, int count_send, T* recv, int count_recv) {
    return MPI_Allgather(data, count_send, GetMPIType<T>(),
                         recv, count_recv, GetMPIType<T>(),
                         communicator);
}

template <typename T>
inline int ExecutionManager::Allgather(const std::vector<T>& data,
                                       std::vector<T>* recv,
                                       bool is_equally_distributed) {
    if (is_equally_distributed) {
        recv->resize(data.size() * GetNumberProcesses());
    } else {
        if (recv->size() == 0) {
            ERROR << "Rank " << std::to_string(GetRank()) << ": receive buffer has size 0! "
                  << "Expect segmentation fault!" << ERROR_CLOSE;
        }
    }
    return Allgather(data.data(), static_cast<int>(data.size()), recv->data(), static_cast<int>(data.size()));
}

template <typename T>
int ExecutionManager::Send(const T* data, int count, int receiver, int tag) {
    return MPI_Send(data, count, GetMPIType<T>(), receiver, tag, communicator);
}

template <typename T>
int ExecutionManager::Isend(const T* data, int count, int receiver, int tag, MPI_Request* request) {
    return MPI_Isend(data, count, GetMPIType<T>(), receiver, tag, communicator, request);
}

template <typename T>
int ExecutionManager::Recv(T* buffer, int count, int sender, int tag, MPI_Status* status) {
    return MPI_Recv(buffer, count, GetMPIType<T>(), sender, tag, communicator, status);
}

template <typename T>
int ExecutionManager::Irecv(T* buffer, int count, int sender, int tag, MPI_Request* request) {
    return MPI_Irecv(buffer, count, GetMPIType<T>(), sender, tag, communicator, request);
}

int ExecutionManager::Probe(int source, int tag, MPI_Status* status) {
    return MPI_Probe(source, tag, communicator, status);
}
}  // namespace dare::mpi

#endif  // MPI_EXECUTIONMANAGER_MPI_INL_
