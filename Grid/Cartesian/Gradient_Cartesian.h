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

#ifndef GRID_CARTESIAN_GRADIENT_CARTESIAN_H_
#define GRID_CARTESIAN_GRADIENT_CARTESIAN_H_

#include "Data/GridVector.h"
#include "Equations/Operators.h"
#include "Grid/Cartesian/Interpolation_Cartesian.h"
#include "Grid/Cartesian/MatrixBlock_Cartesian.h"
#include "Grid/Cartesian/Stencils_Cartesian.h"

namespace dare::Matrix {

/*!
 * @brief gradient operator
 * @tparam Grid type of grid
 * The gradient operator works with a stencil belonging to a certain cell
 */
template <std::size_t Dim>
class Gradient<dare::Grid::Cartesian<Dim>> {
public:
    static const std::size_t NUM_FACES{Dim * 2};                   //!< number of faces
    static const std::size_t NUM_ENTRIES{NUM_FACES + 1};           //!< stencil size
    using GridType = dare::Grid::Cartesian<Dim>;                   //!< type of grid
    using GridRepresentation = typename GridType::Representation;  //!< representation of grid
    using LO = typename GridType::LocalOrdinalType;                //!< convenient aliasing
    using GO = typename GridType::GlobalOrdinalType;               //!< convenient aliasing
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
    template <typename SC, typename O, std::size_t N>
    dare::Data::FaceMatrixStencil<GridType, SC, N>
    operator()(const MatrixBlock<GridType, O, SC, N>& mb) const;

    /*!
     * @brief provide gradient values from a field
     * @tparam N number of components
     * @param field instance of field
     */
    template <typename SC, std::size_t N>
    dare::Data::FaceValueStencil<GridType, SC, N>
    operator()(const dare::Data::GridVector<GridType, SC, N>& field) const;

    /*!
     * @brief provide gradient values from a stencil
     * @tparam N number of components
     * @param s stencil to use
     * \note to users: This is intended to be used for applications like numerical differentiation
     */
    template <typename SC, std::size_t N>
    dare::Data::FaceValueStencil<GridType, SC, N>
    operator()(const dare::Data::CenterValueStencil<GridType, SC, N>& s) const;

    /*!
     * @brief provide gradient values for dedicated component from a field
     * @tparam N number of components
     * @param f field to determine the stencil from
     * @param n component ID
     */
    template <typename SC, std::size_t N>
    dare::Data::FaceValueStencil<GridType, SC, 1>
    operator()(const dare::Data::GridVector<GridType, SC, N>& field, std::size_t n) const;

    /*!
     * @brief provide gradient values for dedicated component from a stencil
     * @tparam N number of components
     * @param s stencil to determine the stencil from
     * @param n component ID
     */
    template <typename SC, std::size_t N>
    dare::Data::FaceValueStencil<GridType, SC, 1>
    operator()(const dare::Data::CenterValueStencil<GridType, SC, N>& s, std::size_t n) const;

private:
    const GridRepresentation* grid;  //!< reference to grid representation
    typename GridType::VecSC dn_r;   //!< recursive distances
    LO ordinal_internal;             //!< local internal ordinal
};

}  // end namespace dare::Matrix

#include "Gradient_Cartesian.inl"

#endif  // GRID_CARTESIAN_GRADIENT_CARTESIAN_H_
