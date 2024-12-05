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

#include <tuple>

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
template <std::size_t Dim, typename TimeDiscretization>
class Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization> {
public:
    static const std::size_t NUM_TIMESTEPS = TimeDiscretization::NUM_TIMESTEPS;
    static const std::size_t NUM_TFIELDS{NUM_TIMESTEPS + 1};
    using GridType = dare::Grid::Cartesian<Dim>;                   // convenient alias
    using LO = typename GridType::LocalOrdinalType;                // convenient alias
    using Index = typename GridType::Index;                        // convenient alias
    using GridRepresentation = typename GridType::Representation;  // convenient alias
    using Positions = typename GridType::NeighborID;               // convenient alias
    template <typename SC, std::size_t N>
    using TFaceMatrixStencil = dare::utils::Vector<NUM_TFIELDS, dare::Data::FaceMatrixStencil<GridType, SC, N>>;
    template <typename SC, std::size_t N>
    using TFaceValueStencil = dare::utils::Vector<NUM_TFIELDS, dare::Data::FaceValueStencil<GridType, SC, N>>;

    /*!
     * @brief constructor
     * @param grid grid representation
     * @param ordinal_internal internal ordinal of the center cell of the stencil
     * \note the 'ordinal_internal' value is only a dummy to be consistent with other
     * specializations. No data based on its value is requested.
     */
    explicit Divergence(const GridRepresentation& grid, LO ordinal_internal);

    template <typename... Args>
    auto operator()(const Args&... args);

private:
    /*!
     * @brief central differencing interpolation of a face value stencil with values from the specified field
     * @tparam SC scalar type
     * @tparam N number of components
     * @param f grid vector
     * @return FaceValueStencil with interpolated data
     */
    template<typename SC, std::size_t N>
    dare::Data::FaceValueStencil<GridType, SC, N>
    PopulateFaceValueFromField(const dare::Data::GridVector<GridType, SC, N>& f) const;

    /*!
     * @brief central differencing interpolation of a face value stencil with values from the specified field
     * @tparam SC scalar type
     * @tparam N number of components
     * @param f field
     * @return TFaceValueStencil with interpolated data
     */
    template <typename SC, std::size_t N>
    TFaceValueStencil<SC, N>
    PopulateFaceValueFromField(const dare::Data::Field<GridType, SC, N>& f) const;

    /*!
     * @brief multiplies with scalar value
     * @tparam SC type of scalar
     * @tparam N number of components
     * @param value value to scale with
     * @param s FaceMatrixStencils
     */
    template <typename SC, std::size_t N>
    void Multiply(SC value, TFaceMatrixStencil<SC, N>* s) const;

    /*!
     * @brief non-conservative multiplication without temporal information
     * @tparam SC type of scalar
     * @tparam N number of components
     * @param f face value stencil
     * @param s FaceMatrixStencils
     */
    template <typename SC, std::size_t N>
    void Multiply(const dare::Data::FaceValueStencil<GridType, SC, N>& f, TFaceMatrixStencil<SC, N>* s) const;

    /*!
     * @brief conservative multiplication with temporal information
     * @tparam SC scalar type
     * @tparam N number of components
     * @param f multiple face value stencils with temporal information
     * @param s FaceMatrixStencils
     */
    template <typename SC, std::size_t N>
    void Multiply(const TFaceValueStencil<SC, N>& f, TFaceMatrixStencil<SC, N>* s) const;

    /*!
     * @brief non-conservative multiplication without temporal information
     * @tparam SC scalar type
     * @tparam N number of components
     * @param f grid vector
     * @param s FaceMatrixStencils
     */
    template <typename SC, std::size_t N>
    void Multiply(const dare::Data::GridVector<GridType, SC, N>& f, TFaceMatrixStencil<SC, N>* s) const;

    /*!
     * @brief conservative multiplication with temporal information
     * @tparam SC scalar type
     * @tparam N number of components
     * @param f field
     * @param s FaceMatrixStencils
     */
    template <typename SC, std::size_t N>
    void Multiply(const dare::Data::Field<GridType, SC, N>& f, TFaceMatrixStencil<SC, N>* s) const;

    /*!
     * @brief loop through all arguments in parameter pack
     * @tparam Stencil stencil which is multiplied
     * @tparam Arg object type that the stencil is multiplied with
     * @tparam ...Args parameter pack with remaining types
     * @param f stencil which is multiplicated with
     * @param arg object with which the stencil is multiplicated
     * @param ...args remaining arguments
     */
    template <typename Stencil, typename Arg, typename... Args>
    void MultiplyAll(Stencil* f, const Arg& arg, const Args&... args);

    template <typename SC, std::size_t N>
    dare::Data::CenterMatrixStencil<GridType, SC, N>
    ApplyDivergence(const TFaceMatrixStencil<SC, N>& s) const;

    template <typename SC, std::size_t N>
    dare::utils::Vector<N, SC>
    ApplyDivergence(const dare::Data::FaceValueStencil<GridType, SC, N>& s) const;

    template<typename SC, std::size_t N>
    dare::Data::FaceMatrixStencil<GridType, SC, N>
    GetFaceMatrixStencil(const dare::Data::Field<GridType, SC, N>& field) const;

    template <typename SC, std::size_t N>
    dare::Data::FaceMatrixStencil<GridType, SC, N>
    GetFaceMatrixStencil(const dare::Data::FaceMatrixStencil<GridType, SC, N>& s) const;

    typename GridType::VecSC A;  //!< face area for each dimension
    Index ind;                   //!< triplet of indices
    const typename GridType::Representation* grep;
};

}  // end namespace dare::Matrix

#include "Divergence_Cartesian.inl"

#endif  // GRID_CARTESIAN_DIVERGENCE_CARTESIAN_H_
