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
template <std::size_t Dim>
class Gradient<dare::Grid::Cartesian<Dim>> {
public:
    static const std::size_t NUM_FACES{Dim * 2};                   //!< number of faces
    static const std::size_t NUM_ENTRIES{NUM_FACES + 1};           //!< stencil size
    using GridType = dare::Grid::Cartesian<Dim>;       //!< type of grid
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
    typename GridType::VecSC A;     //!< face area for each dimension
};

template <std::size_t Dim, typename SC, typename FluxLimiter>
class TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter> {
public:
    static const std::size_t NUM_FACES{2 * Dim};
    using GridType = dare::Grid::Cartesian<Dim>;
    using GridRepresentation = typename GridType::GridRepresentation;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using Options = typename GridType::Options;
    using Positions = typename GridType::NeighborID;

    /*!
     * @brief initialization with references to velocity fields
     * @param grid representation of the target grid
     * @param ordinal_internal internal ordinal
     * @param v vector with velocity fields
     */
    TVD(const GridRepresentation& grid,
        LO ordinal_internal,
        dare::utils::Vector<Dim, const dare::Data::GridVector<GridType, SC, 1>&> v);

    /*!
     * @brief initialization with constant velocities
     * @param grid 
     * @param ordinal_internal 
     * @param v 
     */
    TVD(const GridRepresentation& grid,
        LO ordinal_internal,
        const dare::utils::Vector<Dim, SC>& v);

    /*!
     * @brief destructor
     */
    ~TVD();

    /*!
     * @brief When used for interpolation
     * @tparam N number of components
     * @param s_close CenterMatrixStencil for close neighbors
     * @param s_far CenterValueStencil for far neighbors
     */
    template <std::size_t N>
    dare::Data::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::Data::CenterValueStencil<GridType, SC, N>& s_close,
        const dare::Data::CenterValueStencil<GridType, SC, N>& s_far,
        Options opt) const;

    /*!
     * @brief When used for interpolation from a field
     * @tparam N number of components
     * @param field field
     */
    template <std::size_t N>
    dare::Data::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::Data::GridVector<GridType, SC, N>& field) const;

    /*!
     * \brief when used for matrix assembly
     * @tparam N number of components
     * @param field reference to relevant field
     */
    template <std::size_t N>
    dare::Data::FaceMatrixStencil<GridType, SC, N> operator()(
        const dare::Data::GridVector<GridType, SC, N>& field) const;

private:
    Index ind;                                               //!< triplet of indices
    dare::Data::FaceValueStencil<GridType, SC, 1> velocity;  //!< stencil with velocity
    dare::utils::Vector<NUM_FACES, bool> upwind;             //!< identifier for upwind at each face
    dare::utils::Vector<Dim, bool> self_convection;          //!< identifier, if any direction is self convection
};

}  // end namespace dare::Matrix

#include "Operators_Cartesian.inl"

#endif  // MATRIXSYSTEM_OPERATORS_CARTESIAN_H_
