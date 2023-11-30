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

#ifndef DATA_GRIDVECTOR_H_
#define DATA_GRIDVECTOR_H_

#include <string>
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

namespace dare::Data {

/*!
 * @brief general data storage object
 * @tparam Grid type of grid used for this object
 * @tparam T type of data that is stored
 * @tparam N number of components per grid cell
 */
template<typename Grid, typename T, std::size_t N = 1>
class GridVector {
public:
    using GridType = Grid;
    using GridRepresentation = typename Grid::Representation;
    using LO = typename Grid::LocalOrdinalType;
    using GO = typename Grid::GlobalOrdinalType;
    using Index = typename Grid::Index;
    using IndexGlobal = typename Grid::IndexGlobal;
    using DualViewType = Kokkos::DualView<T*>;
    using DeviceViewType = typename DualViewType::t_dev;
    using HostViewType = typename DualViewType::t_dev;
    using HostSpace = typename DualViewType::host_mirror_space;
    using ExecutionSpace = typename DualViewType::execution_space;

    /*!
     * @brief default constructor
     */
    GridVector();

    /*!
     * @brief initialization constructor
     * @param identifier 
     * @param grid 
     */
    GridVector(std::string identifier, GridRepresentation grid);

    /*!
     * @brief default destructor
     */
    virtual ~GridVector();

    /*!
     * @brief resize to dedicated number
     * @param n size after reallocation
     * \warning After this function call, you are responsible for
     *  legal accessing of the data!
     */
    void Resize(LO n);

    /*!
     * @brief resizes according to number of cells
     * @param n number of cells
     * Multiplies the number with the number of equations
     */
    void ResizeByGrid(LO n);

    /*!
     * @brief access based on cell and component
     * @param n grid ordinal
     * @param c component number
     * @return reference to specified value
     */
    T& At(LO n, std::size_t c);

    /*!
     * @brief access based on cell and component
     * @param ind grid index
     * @param c component number
     * @return reference to specified value
     */
    T& At(const Index& ind, std::size_t c);

    /*!
     * @brief const access based on cell and component
     * @param n grid ordinal
     * @param c component number
     * @return specified value
     */
    T At(LO n, std::size_t c) const;

    /*!
     * @brief const access based on cell and component
     * @param ind grid index
     * @param c component number
     * @return specified value
     */
    T At(const Index& ind, std::size_t c) const;

    /*!
     * @brief access to element in the data
     * @param n position of element
     * @return reference to element
     * \note grid graph is ignored here
     */
    T& operator[](LO n);

    /*!
     * @brief const access to element in the data
     * @param n position of element
     * @return copy of element
     * \note grid graph is ignored here
     */
    T operator[](LO n) const;

    /*!
     * @brief returns number of stored elements
     */
    std::size_t GetSize() const;

    /*!
     * @brief returns number of equations
     * @return 
     */
    constexpr std::size_t GetNumEquations() const { return N; }

    /*!
     * @brief synchronizes host and execution space
     * @tparam TargetSpace identifer, to which space a copy will be send
     */
    template<typename TargetSpace>
    void Synchronize();

    /*!
     * @brief provides deep copy of this instance
     * @return 
     */
    GridVector<Grid, T, N> GetDeepCopy() const;

    /*!
     * @brief provides deep copy to other instance
     * @param other instance to copy to
     */
    void GetDeepCopy(GridVector<Grid, T, N>* other) const;

    /*!
     * @brief returns data on device
     */
    DeviceViewType& GetDeviceView();

    /*!
     * @brief returns data on device
     */
    const DeviceViewType& GetDeviceView() const;

protected:
    GridVector(std::string identifier, LO num_cells, GridRepresentation grid);

private:
    std::string ident_string;       //!< identification string
    GridRepresentation grid;        //!< representation and reference to grid
    DualViewType data;      //!< array with data
};

}  // namespace dare::Data

#include "GridVector.inl"

#endif  // DATA_GRIDVECTOR_H_
