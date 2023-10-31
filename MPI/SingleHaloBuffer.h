/*
 * MIT License
 *
 * Copyright (c) 2022 David Rieder

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

#ifndef MPI_SINGLEHALOBUFFER_H_
#define MPI_SINGLEHALOBUFFER_H_

#include <mpi.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../Utilities/InitializationTracker.h"
#include "ExecutionManager.h"

namespace dare::mpi {

/*! \class SingleHaloBuffer
 * @brief takes care of communication with one processor to exchange halo data
 * This class encapsulates the necessary communication to exchange data in the
 * halo cell region. It stores the rank of a partner process and an array of local
 * grid cell IDs, which it will exchange and subsequently store in the field
 */
template<typename LO, typename GO, typename SC>
class SingleHaloBuffer : public dare::utils::InitializationTracker{
public:
    /*!
     * @brief constructor
     * @param execution_manager pointer to execution manager
     * @param rank_partner rank of partner process
     */
    SingleHaloBuffer(ExecutionManager* execution_manager = nullptr,
               int rank_partner = -1);

    void Initialize(ExecutionManager* execution_manager,
                    int rank_partner);

    void CommunicateAmountRequiredHaloCellIDs(std::size_t num_send, std::size_t* num_recv,
                                              MPI_Request* request_send, MPI_Request* request_receive);

    void CommunicateRequiredHaloCellIDs(const std::vector<GO>& list_global_ID,
                                        std::vector<GO>* list_global_IDs_partner,
                                        MPI_Request* request_send, MPI_Request* request_receive);

    void CommunicateFilteredHaloCellIDs(const std::vector<GO>& list_filtered_ID,
                                        std::vector<GO>* list_filtered_IDs_partner,
                                        MPI_Request* request_send, MPI_Request* request_receive);

    void FinalizeInitialization(const std::vector<LO>& list_IDs_send, const std::vector<LO>& list_IDs_recv);

    /*!
     * @brief 
     * @tparam Field 
     * @param field 
     * @param buffer_send 
     * @param buffer_recv 
     * @param status_send 
     * @param status_recv 
     * @return 
     */
    template <typename Field>
    int Exchange(const Field& field,
                 MPI_Status* status_send, MPI_Status* status_recv);

    /*!
     * @brief 
     * @tparam Field 
     * @param field 
     * @param buffer 
     */
    template <typename Field>
    void FillValues(Field* field);

private:
    ExecutionManager* exec_man;             //!< reference to execution manager
    int rank_partner_proc;                  //!< rank of partner process
    std::vector<LO> list_local_IDs_send;    //!< list of local grid cell IDs to send data from
    std::vector<LO> list_local_IDs_recv;    //!< list of local grid cell IDs to receive data for
    std::vector<SC> buffer_send;            //!< buffer for the scalar values that will be communicated
    std::vector<SC> buffer_recv;            //!< buffer for scalar values that will be received
};

}  // namespace dare::mpi

#include "SingleHaloBuffer.inl"

#endif  //  MPI_SINGLEHALOBUFFER_H_
