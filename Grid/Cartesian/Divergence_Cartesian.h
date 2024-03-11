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

#ifndef GRID_CARTESIAN_DIVERGENCE_CARTESIAN_H_
#define GRID_CARTESIAN_DIVERGENCE_CARTESIAN_H_

#include "Data/GridVector.h"
#include "Grid/Cartesian/Stencils_Cartesian.h"
#include "Grid/Cartesian/Interpolation_Cartesian.h"
#include "Grid/Cartesian/MatrixBlock_Cartesian.h"
#include "Equations/Operators.h"

namespace dare::Matrix {

/*!
 * @brief divergence operator
 * @tparam Grid type of grid
 * \note to users: This operator deals on a cell-wise manner with the stencils
 * By using the Cartesian grid instance, we can optimize the procedure significantly,
 * since we don't need to look up any values besides the face area in each direction.
 * All other required information is encoded by the position of the values.
 *
 * The divergence integrates the flux over the faces. In case of explicit contributions,
 * they are also scaled, with the same sign as the implicit contributions. Therefore,
 * the user has to make sure that the explicit values are properly treated afterwards.
 *
 * \f[
 *      \int_V div(\phi_{face})\,\mathrm{d}V \approx \sum_{nb} \vec{n}_{nb} A_{nb}
 * \f]
 */
template <std::size_t Dim>
class Divergence<dare::Grid::Cartesian<Dim>> {
public:
    using GridType = dare::Grid::Cartesian<Dim>;                   // convenient alias
    using LO = typename GridType::LocalOrdinalType;                // convenient alias
    using GridRepresentation = typename GridType::Representation;  // convenient alias
    using Positions = typename GridType::NeighborID;               // convenient alias

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
    template <typename SC, std::size_t N>
    dare::utils::Vector<N, SC>
    operator()(const dare::Data::FaceValueStencil<GridType, SC, N>& s) const;

    /*!
     * @brief evaluates divergence
     * @tparam N number of components
     * @param s matrix stencil
     */
    template <typename SC, std::size_t N>
    dare::Data::CenterMatrixStencil<GridType, SC, N>
    operator()(const dare::Data::FaceMatrixStencil<GridType, SC, N>& s) const;

private:
    typename GridType::VecSC A;  //!< face area for each dimension
};

}  // end namespace dare::Matrix

#include "Divergence_Cartesian.inl"

#endif  // GRID_CARTESIAN_DIVERGENCE_CARTESIAN_H_
