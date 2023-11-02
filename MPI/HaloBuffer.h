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

#ifndef MPI_HALOBUFFER_H_
#define MPI_HALOBUFFER_H_

#include <map>
#include <unordered_map>
#include <sstream>
#include <vector>

#include "../Utilities/InitializationTracker.h"
#include "SingleHaloBuffer.h"
namespace dare::mpi {

/*! \class HaloBuffer
 * @brief takes care of communication with all processes to exchange halo data
 * This class encapsulates the necessary communication to exchange data in the
 * halo cell region.
 */
template <typename LO, typename GO, typename SC>
class HaloBuffer : public utils::InitializationTracker {
public:
    HaloBuffer();

    template <typename Lambda1, typename Lambda2>
    void Initialize(ExecutionManager* exec_man,
                    const std::vector<GO>& list_required_IDs,
                    const std::unordered_map<GO, GO> map_periodic,
                    Lambda1 is_local,
                    Lambda2 map_global_to_local);

    template <typename Field>
    void Exchange(Field* field);

private:
    ExecutionManager* exec_man;
    std::map<int, SingleHaloBuffer<LO, GO, SC>> buffers;
};
}  // namespace dare::mpi

#include "HaloBuffer.inl"

#endif  // MPI_HALOBUFFER_H_
