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

namespace dare::Data {
template <typename Grid, typename SC, std::size_t N>
Field<Grid, SC, N>::Field(std::string _identifier, GridRepresentation grid_rep, std::size_t num_time_levels)
    : identifier(_identifier), data(num_time_levels) {
    if (num_time_levels == 0) {
        GetExecutionManager()->Terminate(__func__, "Number of time levels may not be 0!");
    }
    for (std::size_t n{0}; n < num_time_levels; n++) {
        data[n].Initialize(identifier + "_" + std::to_string(n), grid_rep);
    }
}

template <typename Grid, typename SC, std::size_t N>
Field<Grid, SC, N>::~Field() {}

template <typename Grid, typename SC, std::size_t N>
typename Field<Grid, SC, N>::VectorType& Field<Grid, SC, N>::GetDataVector(std::size_t time_level) {
#ifndef DARE_NDEBUG
    if (time_level >= data.size()) {
        GetExecutionManager()->Terminate(
            __func__, "time level is too high!");
    }
#endif
    return data[time_level];
}

template <typename Grid, typename SC, std::size_t N>
const typename Field<Grid, SC, N>::VectorType& Field<Grid, SC, N>::GetDataVector(std::size_t time_level) const {
#ifndef DARE_NDEBUG
    if (time_level >= data.size()) {
        GetExecutionManager()->Terminate(
            __func__, "time level is too high!");
    }
#endif
    return data[time_level];
}

template <typename Grid, typename SC, std::size_t N>
std::string Field<Grid, SC, N>::GetFieldName() const {
    return identifier;
}

template <typename Grid, typename SC, std::size_t N>
typename Field<Grid, SC, N>::GridRepresentation& Field<Grid, SC, N>::GetGridRepresentation() {
#ifndef DARE_NDEBUG
    if (data.size() == 0) {
        std::cerr << "ERROR here: " << __func__ << " - Field was not initialized first!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
#endif
    return data[0].GetGridRepresentation();
}

template <typename Grid, typename SC, std::size_t N>
const typename Field<Grid, SC, N>::GridRepresentation& Field<Grid, SC, N>::GetGridRepresentation() const {
#ifndef DARE_NDEBUG
    if (data.size() == 0) {
        std::cerr << "ERROR here: " << __func__ << " - Field was not initialized first!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
#endif
    return data[0].GetGridRepresentation();
}

template <typename Grid, typename SC, std::size_t N>
void Field<Grid, SC, N>::CopyDataVectorsToOldTimeStep() {
    for (std::size_t n{data.size() - 1}; n > 0; n--) {
        data[n - 1].GetDeepCopy(&data[n]);
    }
}

template <typename Grid, typename SC, std::size_t N>
dare::mpi::ExecutionManager* Field<Grid, SC, N>::GetExecutionManager() {
    return GetGridRepresentation().GetHaloBuffer().GetExecutionManager();
}

template <typename Grid, typename SC, std::size_t N>
void Field<Grid, SC, N>::ExchangeHaloCells() {
    GetGridRepresentation().GetHaloBuffer().Exchange(&data[0]);
}
}  // end namespace dare::Data
