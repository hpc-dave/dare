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
#ifndef MATRIXSYSTEM_OPERATORS_CARTESIAN_H_
#define MATRIXSYSTEM_OPERATORS_CARTESIAN_H_

#include "Operators.h"
#include "MatrixBlock_Cartesian.h"
#include "Stencils_Cartesian.h"
#include "../Data/GridVector.h"

namespace dare::Matrix {

/*!
 * @brief gradient operator
 * @tparam Grid type of grid
 * The gradient operator works with a stencil belonging to a certain cell
 */
template <std::size_t Dim, typename LO, typename GO, typename SC>
class Gradient<dare::Grid::Cartesian<Dim, LO, GO, SC>> {
public:
    static const std::size_t NUM_FACES{Dim * 2};                   //!< number of faces
    static const std::size_t NUM_ENTRIES{NUM_FACES + 1};           //!< stencil size
    using GridType = dare::Grid::Cartesian<Dim, LO, GO, SC>;       //!< type of grid
    using GridRepresentation = typename GridType::Representation;  //!< representation of grid
    using Positions = CartesianNeighbor;                           // convenient aliasing
    using Index = typename GridType::Index;                        // convenient aliasing

    /*!
     * @brief constructor
     * @param grid a representation of the grid
     * @param ordinal_internal the internal ordinal of the center cell in the stencil
     */
    Gradient(const GridRepresentation& grid, LO ordinal_internal);

    /*!
     * @brief provide matrix stencil
     * @tparam O ordinal type used in the MatrixBlock instance
     * @tparam N number of components
     * @param mb dummy parameter
     * By using a matrix block as an argument, we can determine the number of components in a
     * convenient way without too much code smell.
     */
    template <typename O, std::size_t N>
    dare::Data::FaceMatrixStencil<GridType, N>
    operator()(const MatrixBlock<GridType, O, SC, N>& mb) const;

    /*!
     * @brief provide gradient values from a field
     * @tparam N number of components
     * @param field instance of field
     */
    template <std::size_t N>
    dare::Data::FaceValueStencil<GridType, N>
    operator()(const dare::Data::GridVector<GridType, SC, N>& field) const;

    /*!
     * @brief provide gradient values from a stencil
     * @tparam N number of components
     * @param s stencil to use
     * \note to users: This is intended to be used for applications like numerical differentiation
     */
    template <std::size_t N>
    dare::Data::FaceValueStencil<GridType, N>
    operator()(const dare::Data::CenterValueStencil<GridType, N>& s) const;

    /*!
     * @brief provide gradient values for dedicated component from a field
     * @tparam N number of components
     * @param f field to determine the stencil from
     * @param n component ID
     */
    template <std::size_t N>
    dare::Data::FaceValueStencil<GridType, 1>
    operator()(const dare::Data::GridVector<GridType, SC, N>& field, std::size_t n) const;

    /*!
     * @brief provide gradient values for dedicated component from a stencil
     * @tparam N number of components
     * @param s stencil to determine the stencil from
     * @param n component ID
     */
    template <std::size_t N>
    dare::Data::FaceValueStencil<GridType, 1>
    operator()(const dare::Data::CenterValueStencil<GridType, N>& s, std::size_t n) const;

private:
    const GridRepresentation* grid;         //!< reference to grid representation
    typename GridType::VecSC dn_r;          //!< recursive distances
    LO ordinal_internal;                    //!< local internal ordinal
};

/*!
 * @brief divergence operator
 * @tparam Grid type of grid
 * \note to users: This operator deals on a cell-wise manner with the stencils
 * By using the Cartesian grid instance, we can optimize the procedure significantly,
 * since we don't need to look up any values besides the face area in each direction.
 * All other required information is encoded by the position of the values.
 */
template <std::size_t Dim, typename LO, typename GO, typename SC>
class Divergence<dare::Grid::Cartesian<Dim, LO, GO, SC>> {
public:
    using GridType = dare::Grid::Cartesian<Dim, LO, GO, SC>;       // convenient alias
    using GridRepresentation = typename GridType::Representation;  // convenient alias
    using Positions = dare::Data::CartesianNeighbor;               // convenient alias

    /*!
     * @brief constructor
     * @param grid grid representation
     * @param ordinal_internal internal ordinal of the center cell of the stencil
     * \note the 'ordinal_internal' value is only a dummy to be consistent with other
     * specializations. No data based on its value is requested. 
     */
    explicit Divergence(const GridRepresentation& grid, LO ordinal_internal = 0);

    /*!
     * @brief evaluates divergence
     * @tparam N number of components
     * @param s value stencil
     */
    template <std::size_t N>
    dare::utils::Vector<N, SC>
    operator()(const dare::Data::FaceValueStencil<GridType, N>& s) const;

    /*!
     * @brief evaluates divergence
     * @tparam N number of components
     * @param s matrix stencil
     */
    template <std::size_t N>
    dare::Data::CenterMatrixStencil<GridType, N>
    operator()(const dare::Data::FaceMatrixStencil<GridType, N>& s) const;

private:
    typename GridType::VecSC A;     //!< face area for each dimension
};

}  // end namespace dare::Matrix

#include "Operators_Cartesian.inl"

#endif  // MATRIXSYSTEM_OPERATORS_CARTESIAN_H_
