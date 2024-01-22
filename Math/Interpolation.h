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

#ifndef MATH_INTERPOLATION_H_
#define MATH_INTERPOLATION_H_

#include <limits>
#include <typeinfo>

#include "Data/GridVector.h"
#include "Pow.h"
#include "Utilities/Errors.h"

namespace dare::math {

/*!
 * @brief Interpolates to a specified face of the target grid
 * @tparam GridType type of grid
 * @tparam SC type of scalar
 * @tparam N number of components
 * @param target target grid
 * @param ind_target indices of the target cell
 * @param field field with values
 * @param n component ID
 * For some grids, optimization can be gained, if a certain face ID is known, e.g. WEST or EAST
 * on a Cartesian grid. In those cases, the weight determination is trivial.
 * \note This is merely a placeholder, for which an overload should be defined by the grid,
 * otherwise a NaN is returned
 */
template <typename GridType, typename SC, std::size_t N>
[[nodiscard]] SC InterpolateToFace(const typename GridType::Representation& target,
                                   const typename GridType::Index& ind_target,
                                   const typename GridType::NeighborID face,
                                   const Data::GridVector<GridType, SC, N>& field,
                                   std::size_t n) {
    ERROR << "This function is not implemeted for the GridType: " << typeid(GridType).name() << ERROR_CLOSE;
    return std::numeric_limits<SC>::signaling_NaN();
}

/*!
 * @brief Interpolates to a specified face of the target grid
 * @tparam GridType type of grid
 * @tparam SC type of scalar
 * @tparam N number of components
 * @param target target grid
 * @param ind_target indices of the target cell
 * @param field field with values
 * @param n component ID
 * For some grids, optimization can be gained, if a certain face ID is known, e.g. WEST or EAST
 * on a Cartesian grid. In those cases, the weight determination is trivial.
 * \note This is merely a placeholder, for which an overload should be defined by the grid,
 * otherwise a NaN is returned
 */
template <typename GridType, typename SC, std::size_t N>
[[nodiscard]] dare::utils::Vector<N, SC> InterpolateToFace(const typename GridType::Representation& target,
                                                           const typename GridType::Index& ind_target,
                                                           const typename GridType::NeighborID face,
                                                           const Data::GridVector<GridType, SC, N>& field) {
    ERROR << "This function is not implemeted for the GridType: " << typeid(GridType).name() << ERROR_CLOSE;
    return dare::utils::Vector<N, SC>(std::numeric_limits<SC>::signaling_NaN());
}

/*!
 * @brief Interpolates to a specified point in the field
 * @tparam GridType type of grid
 * @tparam SC type of scalar
 * @tparam N number of components
 * @param point target grid
 * @param field field with values
 * @param n component ID
 * This is the workhorse for interpolation. For each grid, we need to be able to request the values at
 * arbitrary positions inside the grid.
 * \note This is merely a placeholder, for which an overload should be defined by the grid,
 * otherwise a NaN is returned
 */
template <typename GridType, typename SC, std::size_t N>
[[nodiscard]] SC InterpolateToPoint(const typename GridType::VecSC& point,
                      const Data::GridVector<GridType, SC, N>& field,
                      std::size_t n) {
    ERROR << "This function is not implemeted for the GridType: " << typeid(GridType).name() << ERROR_CLOSE;
    return std::numeric_limits<SC>::signaling_NaN();
}

/*!
 * @brief Interpolates to a specified point in the field
 * @tparam GridType type of grid
 * @tparam SC type of scalar
 * @tparam N number of components
 * @param point target grid
 * @param field field with values
 * This is the workhorse for interpolation. For each grid, we need to be able to request the values at
 * arbitrary positions inside the grid.
 * \note This is merely a placeholder, for which an overload should be defined by the grid,
 * otherwise a NaN is returned
 */
template <typename GridType, typename SC, std::size_t N>
[[nodiscard]] dare::utils::Vector<N, SC> InterpolateToPoint(const typename GridType::VecSC& point,
                                                            const Data::GridVector<GridType, SC, N>& field) {
    ERROR << "This function is not implemeted for the GridType: " << typeid(GridType).name() << ERROR_CLOSE;
    return dare::utils::Vector<N, SC>(std::numeric_limits<SC>::signaling_NaN());
}

/*!
 * @brief Computes interpolated value according to indices and weights
 * @tparam GridType type of grid
 * @tparam T type of value
 * @tparam N number of components
 * @param field field with values
 * @param indices a list of indices, correlating with the weights
 * @param weights weight for each index
 * @param n_component component ID
 * The values are computed according to
 * \f[
 *      v \gets \sum_n \omega_n \cdot F(I_n)
 * \f]
 * with Index I and field F.
 */
template<typename GridType, typename T, std::size_t N, std::size_t NUM_VALUES>
[[nodiscard]] T Interpolate(const dare::Data::GridVector<GridType, T, N>& field,
               const dare::utils::Vector<NUM_VALUES, typename GridType::Index>& indices,
               const dare::utils::Vector<NUM_VALUES, T>& weights,
               std::size_t n_component) {
    T v{0};
    for (std::size_t n{0}; n < NUM_VALUES; n++) {
        v = std::fma(weights[n], field.At(indices[n], n_component), v);
    }
    return v;
}

/*!
 * @brief Computes interpolated value according to indices and weights
 * @tparam GridType type of grid
 * @tparam T type of value
 * @tparam N number of components
 * @param field field with values
 * @param indices a list of indices, correlating with the weights
 * @param weights weight for each index
 * The values are computed according to
 * \f[
 *      v \gets \sum_n \omega_n \cdot F(I_n)
 * \f]
 * with Index I and field F.
 */
template <typename GridType, typename T, std::size_t N, std::size_t NUM_VALUES>
[[nodiscard]] dare::utils::Vector<N, T> Interpolate(const dare::Data::GridVector<GridType, T, N>& field,
                                      const dare::utils::Vector<NUM_VALUES, typename GridType::Index>& indices,
                                      const dare::utils::Vector<NUM_VALUES, T>& weights) {
    dare::utils::Vector<N, T> v;
    v.SetAllValues(0);
    for (std::size_t n{0}; n < NUM_VALUES; n++) {
        auto v_field = field.GetValues(indices[n]);
        for (std::size_t n_c{0}; n_c < N; n_c++)
            v[n_c] = std::fma(weights[n], v_field[n_c], v[n_c]);
    }
    return v;
}



}  // end namespace dare::math

#endif  // MATH_INTERPOLATION_H_
