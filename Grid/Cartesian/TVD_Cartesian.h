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
#ifndef EQUATIONS_TVD_CARTESIAN_H_
#define EQUATIONS_TVD_CARTESIAN_H_

#include "Data/GridVector.h"
#include "Data/Stencils_Cartesian.h"
#include "Grid/Cartesian/Interpolation_Cartesian.h"
#include "MatrixSystem/MatrixBlock_Cartesian.h"
#include "Operators.h"

namespace dare::Matrix {

template <std::size_t Dim, typename SC, typename FluxLimiter>
class TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter> {
public:
    static const std::size_t NUM_FACES{2 * Dim};
    using GridType = dare::Grid::Cartesian<Dim>;
    using GridRepresentation = typename GridType::Representation;
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
        dare::utils::Vector<Dim, const dare::Data::GridVector<GridType, SC, 1>*> v);

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
     * Here, the face values are computed according to the TVD scheme,
     * excluding the velocity component
     */
    template <std::size_t N>
    [[nodiscard]] dare::Data::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::Data::CenterValueStencil<GridType, SC, N>& s_close,
        const dare::Data::CenterValueStencil<GridType, SC, N>& s_far) const;

    /*!
     * @brief When used for interpolation from a field
     * @tparam N number of components
     * @param field field
     * Here, the face values are computed according to the TVD scheme,
     * excluding the velocity component
     */
    template <std::size_t N>
    [[nodiscard]] dare::Data::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::Data::GridVector<GridType, SC, N>& field) const;

    /*!
     * \brief when used for matrix assembly
     * @tparam N number of components
     * @param field reference to relevant field
     * Here, the flux is computed by the TVD scheme, including the velocity component.
     * \f[
     *    J = u \cdot \phi
     * \f]
     * An example code can look like this:
     * @code{.cpp}
     *   GridVector<...> phi;
     *   TVD<...> u(...);
     *   FaceValueStencil<...> rho;
     *
     *   FaceMatrixStenci<...> J = rho * u * phi
     *
     * @endcode
     */
    template <std::size_t N>
    [[nodiscard]] dare::Data::FaceMatrixStencil<GridType, SC, N> operator*(
        const dare::Data::GridVector<GridType, SC, N>& field) const;

    /*!
     * @brief returns the velocity values
     */
    [[nodiscard]] const dare::Data::FaceValueStencil<GridType, SC, 1>& GetVelocities() const;

private:
    Index ind;                                               //!< triplet of indices
    dare::Data::FaceValueStencil<GridType, SC, 1> velocity;  //!< stencil with velocity
    dare::utils::Vector<NUM_FACES, bool> upwind;             //!< identifier for upwind at each face
};

}  // end namespace dare::Matrix

#include "TVD_Cartesian.inl"

#endif  // EQUATIONS_TVD_CARTESIAN_H_
