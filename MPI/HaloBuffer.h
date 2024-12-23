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

#ifndef MPI_HALOBUFFER_H_
#define MPI_HALOBUFFER_H_

#include <map>
#include <unordered_map>
#include <sstream>
#include <vector>

#include "Utilities/InitializationTracker.h"
#include "Data/DefaultTypes.h"
#include "SingleHaloBuffer.h"
namespace dare::mpi {

/*! \class HaloBuffer
 * @brief takes care of communication with all processes to exchange halo data
 * This class encapsulates the necessary communication to exchange data in the
 * halo cell region.
 */
template <typename SC>
class HaloBuffer : public utils::InitializationTracker {
public:
    using LO = dare::defaults::LocalOrdinalType;
    using GO = dare::defaults::GlobalOrdinalType;

    /*!
     * @brief default constructor
     */
    HaloBuffer();

    /*!
     * @brief prepares all required data
     * @tparam Lambda1 lambda of form lambda(GO):bool
     * @tparam Lambda2 lambda of form lambda(LO):bool
     * @tparam Lambda3 lambda of form lambda(GO):LO
     * @param exec_man pointer to execution manager
     * @param list_required_IDs list of IDs which are located in halo cells
     * @param map_periodic map for correlation of periodic cells
     * @param is_internal_glob lambda to determine if an ID is located on the internal part of the subdomain using global ordinals
     * @param is_internal_loc lambda to determine if an ID is located on the internal part of the subdomain using local ordinals
     * @param map_global_to_local lambda to convert global to local ID
     *
     * This function prepares a buffer for communication with each process within the
     * communicator, including itself. <b> IMPORTANT: <\b> the <i> list_of_required_IDs <\i>
     * refers to the global IDs in the halo/ghost cell area! DO NOT CORRECT THEM for
     * periodic arrangements. This will be conducted internally with the help of the
     * <i> map_periodic </i> object.
     *
     * Communication with the local process is a convenient solution for copying
     * the data within local periodic cells, so the user does not have to care
     * about that.
     */
    template <typename Lambda1, typename Lambda2, typename Lambda3>
    void Initialize(ExecutionManager* exec_man,
                    const std::vector<GO>& list_required_IDs,
                    const std::unordered_map<GO, GO> map_periodic,
                    Lambda1 is_internal_glob,
                    Lambda2 is_internal_loc,
                    Lambda3 map_global_to_local);

    /*!
     * @brief conducts communication to exchange halo cells
     * @tparam Field field type
     * @param field pointer to field
     * \note the Field type requires an interface with two
     * functions:
     *  - GetNumEquations(): std::size_t - number of equations per grid cell
     *  - at(local_id: LO): SC& - access element of field at local_id
     *  - at(local_id: LO): SC const - constant version
     */
    template <typename Field>
    void Exchange(Field* field);

    /*!
     * @brief returns execution manager
     */
    ExecutionManager* GetExecutionManager();

private:
    ExecutionManager* exec_man;                           //!< execution manager
    std::map<int, SingleHaloBuffer<SC>> buffers;          //!< list of buffers
};
}  // namespace dare::mpi

#include "HaloBuffer.inl"

#endif  // MPI_HALOBUFFER_H_
