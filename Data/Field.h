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

#ifndef DATA_FIELD_H_
#define DATA_FIELD_H_

#include <vector>
#include <string>

#include "GridVector.h"
#include "../MPI/HaloBuffer.h"

namespace dare::Data {

template <typename Grid, typename SC, std::size_t N>
class Field {
public:
    using GridType = Grid;
    using LocalOrdinalType = typename GridType::LocalOrdinalType;
    using IndexType = typename GridType::Index;
    using GridRepresentation = typename Grid::Representation;
    using ScalarType = SC;
    using VectorType = GridVector<Grid, SC, N>;

    static const unsigned int version = 0;  //!< used for serialization

    /*!
     * @brief initializing constructor
     * @param identifier string to identify the field
     * @param grid_rep representation of the grid with spatial information
     * @param num_time_fields number of time fields
     */
    Field(std::string identifier, GridRepresentation grid_rep, std::size_t num_time_levels);

    virtual ~Field();

    VectorType& GetDataVector(std::size_t time_level = 0);
    const VectorType& GetDataVector(std::size_t time_level = 0) const;

    std::string GetFieldName() const;

    GridRepresentation& GetGridRepresentation();
    const GridRepresentation& GetGridRepresentation() const;

    dare::mpi::ExecutionManager* GetExecutionManager();

    void CopyDataVectorsToOldTimeStep();

    void ExchangeHaloCells();

private:
    std::string identifier;         //!< name of the field
    std::vector<VectorType> data;   //!< data of the field
};

}  // end namespace dare::Data

#include "Field.inl"

#endif  // DATA_FIELD_H_
