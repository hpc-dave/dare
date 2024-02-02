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
#include <limits>

#include "GridVector.h"
#include "MPI/HaloBuffer.h"
#include "Utilities/Errors.h"

namespace dare::Data {

/*!
 * @brief manages field related data and stores time steps
 * @tparam Grid type of grid
 * @tparam SC type of scalar
 * @tparam N number of components
 * This Field can store arbitrary amounts of timesteps. Recommended use is to
 * utilize the grid at level 0 as current timestep and subsequently store
 * timesteps at levels 1, 2, ...
 */
template <typename Grid, typename SC, std::size_t N>
class Field {
public:
    using GridType = Grid;
    using LocalOrdinalType = typename GridType::LocalOrdinalType;
    using IndexType = typename GridType::Index;
    using GridRepresentation = typename Grid::Representation;
    using ScalarType = SC;
    using VectorType = GridVector<Grid, SC, N>;

    static const std::size_t NUM_COMPONENTS{N};  //!< number of components
    static const unsigned int version = 0;       //!< used for serialization

    /*!
     * @brief initializing constructor
     * @param identifier string to identify the field
     * @param grid_rep representation of the grid with spatial information
     * @param num_time_fields number of time fields
     */
    Field(std::string identifier, GridRepresentation grid_rep, std::size_t num_time_levels);

    /*!
     * @brief default destructor
     */
    virtual ~Field();

    /*!
     * @brief Getter for fields
     * @param time_level time level to access
     */
    VectorType& GetDataVector(std::size_t time_level = 0);

    /*!
     * @brief const getter for fields
     * @param time_level time level to access
     */
    const VectorType& GetDataVector(std::size_t time_level = 0) const;

    /*!
     * @brief provides field name for identification
     */
    std::string GetFieldName() const;

    /*!
     * @brief getter for grid representation
     */
    GridRepresentation& GetGridRepresentation();

    /*!
     * @brief const getter for GridRepresentation
     */
    const GridRepresentation& GetGridRepresentation() const;

    /*!
     * @brief getter for execution manager
     */
    dare::mpi::ExecutionManager* GetExecutionManager();

    /*!
     * @brief initializes cascade of copy steps
     */
    void CopyDataVectorsToOldTimeStep();

    /*!
     * @brief Exchanges halo cells at level 0
     */
    void ExchangeHaloCells();

    /*!
     * @brief provides number of allocated timesteps
     */
    std::size_t GetNumberTimesteps() const;

    void SetValues(SC v, std::size_t time_step = std::numeric_limits<std::size_t>::max());

private:
    std::string identifier;         //!< name of the field
    std::vector<VectorType> data;   //!< data of the field
};

}  // end namespace dare::Data

namespace dare {

/*!
 * @brief type trait to determine if type is a Field
 * @tparam T type to test
 * This is the SFINAE option for false
 */
template <typename T>
struct is_field : std::false_type{
};

/*!
 * @brief type trait to determine if type is a Field
 * @tparam T type to test
 * This is the SFINAE option for true
 */
template<typename Grid, typename SC, std::size_t N>
struct is_field<Data::Field<Grid, SC, N>> : std::true_type{
};

template <typename Grid, typename SC, std::size_t N>
struct is_field<const Data::Field<Grid, SC, N>> : std::true_type {
};

template <typename T>
static const bool is_field_v = is_field<T>::value;

}  // end namespace dare

#include "Field.inl"

#endif  // DATA_FIELD_H_
