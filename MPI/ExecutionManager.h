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

#ifndef MPI_EXECUTIONMANAGER_H_
#define MPI_EXECUTIONMANAGER_H_

#include <mpi.h>
#include <omp.h>

#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "MPITypeConverter.h"
#include "OBlackHoleStream.h"
#include "Utilities/Errors.h"

namespace dare::mpi {

/*! \class Verbosity
 * \brief identifiers for verbosity
 * Defines various level of output, which can be
 * set via the ExecutionManager.
 * \note None is meant for swallowing all output!
 */
enum class Verbosity : uint8_t {
    None,    //!< All output is swallowed
    Low,     //!< minimum production output
    Medium,  //!< a bit more information for rough debugging
    High     //!< huge amount of ouput, only for detailed debuggin
};

/*! \class ExecutionManager
 * \brief takes care of different ways of execution (parallel, multithreaded, accelerated)
 */
class ExecutionManager {
public:
    /*!
     * \brief constructor
     * @param communicator MPI communicator
     * @param output_level controls the output levels
     * \note to developers: Initiates communication, therefore not threadsafe and costly
     */
    explicit ExecutionManager(MPI_Comm communicator = MPI_COMM_WORLD, Verbosity output_level = Verbosity::Low);

    /*!
     * \brief default desctructor
     */
    ~ExecutionManager();

    /*!
     * \brief allows printing with output control
     * @param level level below which the output will be swallowed
     * Messages provided via this function will only be printed
     * by the root processor and if the \p level is lower than
     * the internally specified output verbosity
     * \note not threadsafe, as all output to terminal
     */
    std::ostream& operator()(Verbosity level);

    /*!
     * \brief allows printing with output control
     * @param level level below which the output will be swallowed
     * Messages provided via this function will only be printed
     * by the root processor and if the \p level is lower than
     * the internally specified output verbosity
     * \note not threadsafe, as all output to terminal
     */
    std::ostream& Print(Verbosity level);

    /*!
     * \brief prints output of all processes
     * @param level level below which the output will be swallowed
     * Messages provided via this function will only be printed
     * if the \p level is lower than the internally specified output verbosity,
     * but it will print for each processor
     * \note not threadsafe, as all output to terminal
     */
    std::ostream& PrintAll(Verbosity level);

    /*!
     * \brief returns number of processes within communicator
     * \note threadsafe, if not accesssed via RCP or similar constructs
     */
    inline int GetNumberProcesses() const;

    /*!
     * \brief returns rank within communicator
     * \note threadsafe, if not accesssed via Teuchos::RCP or similar constructs
     */
    inline int GetRank() const;

    /*!
     * \brief returns true, if process is root in communicator
     * \note threadsafe, if not accessed via Teuchos::RCP or similar constructs
     */
    inline bool AmIRoot() const;

    /*!
     * \brief returns true, if execution is serial within communicator
     * \note threadsafe, if not accessed via Teuchos::RCP or similar constructs
     */
    inline bool IsSerial() const;

    /*!
     * \brief sets the output verbosity level
     * @param level verbosity level
     */
    inline void SetVerbosity(Verbosity level);

    /*!
     * \brief returns the number of threads on local machine
     * \note threadsafe, if not accessed via Teuchos::RCP or similar constructs
     *
     * Calls the OpenMP routine, which returns the same value for the whole machine
     */
    int GetNumberThreadsLocal();

    /*!
     * \brief returns the average number of threads
     * \note invokes communication and returns the average of threads per machine
     */
    int GetNumberThreadsAverage();

    /*!
     * \brief terminates execution
     * @param function name of the function (e.g. via __func__)
     * @param message will be printed to terminal
     */
    void Terminate(std::string function, std::string message) const;

    /*!
     * \brief wrapper around multiple MPI routines to exchange buffers
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to array with values to reduce
     * @param[in] count_send number of values to send to other process
     * @param[out] recv pointer to array which will store the reduced values
     * @param[in] count_recv number of values to expect form other process
     * @param[in] rank_other rank of the receiver
     * @param[in] tag tag to use for the communication
     * @param[out] status_send MPI-status for the sending
     * @param[out] status_recv MPI-status for receiving
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    void Exchange(const T* data, int count_send,
                  T* recv, int count_recv,
                  int rank_other, int tag,
                  MPI_Status* status_send = MPI_STATUS_IGNORE, MPI_Status* status_recv = MPI_STATUS_IGNORE);

    /*!
     * \brief Overload to deal with unknown incoming data size
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to array with values to reduce
     * @param[in] count_send number of values to send to other process
     * @param[out] recv pointer to resizable array which will store the reduced values
     * @param[in] rank_other rank of the receiver
     * @param[in] tag tag to use for the communication
     * @param[out] status_send MPI-status for the sending
     * @param[out] status_recv MPI-status for receiving
     * \note invokes communication, not threadsafe and more costly than previous overload
     */
    template <typename T>
    void Exchange(const T* data, int count_send,
                  std::vector<T>* recv,
                  int rank_other, int tag,
                  MPI_Status* status_send = MPI_STATUS_IGNORE, MPI_Status* status_recv = MPI_STATUS_IGNORE);

    /*!
     * \brief Overload to deal with unknown incoming data size and vector input
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to array with values to reduce
     * @param[in] count_send number of values to send to other process
     * @param[out] recv pointer to resizable array which will store the reduced values
     * @param[in] rank_other rank of the receiver
     * @param[in] tag tag to use for the communication
     * @param[out] status_send MPI-status for the sending
     * @param[out] status_recv MPI-status for receiving
     * \note invokes communication, not threadsafe and more costly than previous overload
     */
    template <typename T>
    void Exchange(const std::vector<T>& send,
                  std::vector<T>* recv,
                  int rank_other, int tag,
                  MPI_Status* status_send = MPI_STATUS_IGNORE, MPI_Status* status_recv = MPI_STATUS_IGNORE);

    /*!
     * \brief wrapper around multiple non-blocking MPI routines to exchange buffers
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to array with values to reduce
     * @param[in] count_send number of values to send to other process
     * @param[out] recv pointer to array which will store the reduced values
     * @param[in] count_recv number of values to expect form other process
     * @param[in] rank_other rank of the receiver
     * @param[in] tag tag to use for the communication
     * @param[out] request_send MPI-status for the sending
     * @param[out] request_recv MPI-status for receiving
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    void Iexchange(const T* data, int count_send,
                  T* recv, int count_recv,
                  int rank_other, int tag,
                  MPI_Request* request_send, MPI_Request* request_recv);

    /*!
     * \brief invokes barrier if required
     * \note requires communication, not threadsafe
     */
    void Barrier() const;

    /*!
     * \brief wrapper around the MPI-routine of the same name for a single value
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param data value to reduce
     * @param op specifier for the reduce operation
     * @return reduced value
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    T Allreduce(const T data, MPI_Op op);

    /*!
     * \brief wrapper around the MPI-routine of the same name for an array of values
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to array with values to reduce
     * @param[out] recv pointer to array which will store the reduced values
     * @param[in] count number of values to reduce
     * @param[in] op specifier for the reduce operation
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    void Allreduce(const T* data, T* recv, int count, MPI_Op op);

    /*!
     * \brief sums up values over all processes and spreads information
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param data value to reduce
     * @return reduced value
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    T Allsum(const T data);

    /*!
     * \brief sums up multiple values over all processes and spreads information
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to the array to sum up
     * @param[out] recv pointer to the array whcih will store the summed up values
     * @param[in] count number of values in the array
     *
     * For clarity, each of the values is summed up independently over all processes,
     * so they won't be added to each other within a process!
     *
     * @return reduced value
     *
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    void Allsum(const T* data, T* recv, int count);

    /*!
     * \brief finds maximum on all processes and spreads information
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param data value to reduce
     * @return reduced value
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    T Allmax(const T data);

    /*!
     * \brief find maximum of multiple values over all processes and spreads information
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to the array to sum up
     * @param[out] recv pointer to the array whcih will store the summed up values
     * @param[in] count number of values in the array
     *
     * For clarity, maximum of each of the values is found independently over all processes,
     * so they won't be compared withtin a process!
     *
     * @return reduced value
     *
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    void Allmax(const T* data, T* recv, int count);

    /*!
     * \brief finds minimum on all processes and spreads information
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param data value to reduce
     * @return reduced value
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    T Allmin(const T data);

    /*!
     * \brief find minimum of multiple values over all processes and spreads information
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to the array to sum up
     * @param[out] recv pointer to the array whcih will store the summed up values
     * @param[in] count number of values in the array
     *
     * For clarity, minimum of each of the values is found independently over all processes,
     * so they won't be compared withtin a process!
     *
     * @return reduced value
     *
     * \note invokes communication, not threadsafe
     */
    template <typename T>
    void Allmin(const T* data, T* recv, int count);

    /*!
     * \brief performs logic and on all processes and spreads information
     * @param data value to reduce
     * @return reduced value
     * \note invokes communication, not threadsafe
     */
    inline bool AllLogicAnd(const bool data);

    /*!
     * \brief performs logic and of multiple values over all processes and spreads information
     * @param[in] data pointer to the array to sum up
     * @param[out] recv pointer to the array whcih will store the summed up values
     * @param[in] count number of values in the array
     *
     * For clarity, minimum of each of the values is found independently over all processes,
     * so they won't be compared withtin a process!
     *
     * @return reduced value
     *
     * \note invokes communication, not threadsafe
     */
    inline void AllLogicAnd(const bool* data, bool* recv, int count);

    /*!
     * \brief blocking MPI send operation
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to the array of data to send
     * @param[in] count number of elements to send
     * @param[in] receiver process to send the data to
     * @param[in] tag identifier of the send operation
     * @return status of success
     */
    template <typename T>
    int Send(const T* data, int count, int receiver, int tag);

    /*!
     * \brief non-blocking MPI send operation
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[in] data pointer to the array of data to send
     * @param[in] count number of elements to send
     * @param[in] receiver process to send the data to
     * @param[in] tag identifier of the send operation
     * @param[out] request identifier for inquiry of the sending status
     * @return status of success
     */
    template <typename T>
    int Isend(const T* data, int count, int receiver, int tag, MPI_Request* request);

    /*!
     * \brief blocking MPI receive operation
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[out] buffer pointer to the array in which the data will be stored
     * @param[in] count number of elements to receive
     * @param[in] sender process which sends the message
     * @param[in] tag identifier of the send operation
     * @param[out] status status identifer, if receive was successful
     * @return status of success
     */
    template <typename T>
    int Recv(T* buffer, int count, int sender, int tag, MPI_Status* status = MPI_STATUS_IGNORE);

    /*!
     * \brief non-blocking MPI receive operation
     * \tparam T type of data, needs to be convertible via GetMPIType()
     * @param[out] buffer pointer to the array in which the data will be stored
     * @param[in] count number of elements to receive
     * @param[in] sender process which sends the message
     * @param[in] tag identifier of the send operation
     * @param[out] request identifer for inquire of receive status
     * @return status of success
     */
    template <typename T>
    int Irecv(T* buffer, int count, int sender, int tag, MPI_Request* request);

    /*!
     * @brief Probes a tag
     * @param source rank of sending process
     * @param tag tag of the message
     * @param status poitner to message
     * @return MPI error
     */
    inline int Probe(int source, int tag, MPI_Status* status);

    /*!
     * @brief Stops execution until all requests are satisfied
     * @param requests list of requests
     * @param status list of statuses
     */
    int Waitall(std::vector<MPI_Request>& requests, std::vector<MPI_Status>& status); // NOLINT

    /*!
     * @brief stops execution until all requests are satisfied
     * @param requests list of all requests
     * @param status pointer to array with statuses
     * @return 
     */
    int Waitall(std::vector<MPI_Request>& requests, MPI_Status* status = MPI_STATUSES_IGNORE);  // NOLINT

    /*!
     * @brief stops execution until all requests are satisfied
     * @param count number of entries
     * @param requests array of requests
     * @param status array of statuses
     */
    int Waitall(int count, MPI_Request* requests, MPI_Status* status = MPI_STATUSES_IGNORE);

    /*!
     * \brief returns internal MPI-communicator
     */
    inline MPI_Comm GetCommunicator() const;

    /*!
     * @brief returns rank of root process
     * @return rank of root process
     */
    inline int GetRankRoot() const;

private:
    MPI_Comm communicator;  //!< mpi communicator of this execution manager
    bool is_root;           //!< identifier, if associated process is root
    bool is_serial;         //!< identifier, if code is executed in serial
    int rank_root;          //!< rank of root processor
    int rank;               //!< rank of process within communicator
    int num_proc;           //!< number of processes within communicator

    Verbosity output_level;  //!< all output at levels above is swallowed

    BlackHoleOStream black_hole_osteam;  //!< can swallow output if required
};

}  // namespace dare::mpi

#include "ExecutionManager.inl"
#include "ExecutionManager_MPI.inl"

#endif  // MPI_EXECUTIONMANAGER_H_
