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

#ifndef MATH_INTERPOLATION_CARTESIAN_H_
#define MATH_INTERPOLATION_CARTESIAN_H_

#include "Interpolation.h"
#include "Divisors.h"
#include "../Grid/Cartesian.h"

namespace dare::math {

namespace details::Cartesian {

/*!
 * @brief provides some consistency between the following functions
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 * @tparam SC scalar value type
 * @tparam Dim dimension of grid
 * @tparam DimProj dimension of the projection for the interpolation
 */
template <std::size_t Dim, std::size_t DimProj>
struct THelper {
    static const NUM_VALUES{math::Pow<2, DimProj>()};
    using GridType = dare::Grid::Cartesian<Dim>;
    using Index = typename GridType::Index;
    using IndexList = dare::utils::Vector<NUM_VALUES, Index>;
    using Options = typename GridType::Options;
    template <typename SC, std::size_t N>
    using GridVector = dare::data::GridVector<GridType, SC, N>;
};

/*!
 * @brief compute the indices required to interpolate at special points in the grid
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 * @tparam SC scalar typ
 * @tparam Dim dimension of grid
 * @tparam DimProj dimension of the projecdtion
 * @param ind triplet of indices at the center of the relevant cell
 * @param off_rel relative offset to center in integral values
 * @param dim_aff array with the relevant dimensions
 * @return List of indices for interpolation
 *
 * This is a special function to determine the indices required for the interpolation on a
 * Cartesian grid. It is assumed, that the relative offset is given as integral steps, e.g.
 *
 *  x
 *  --->
 *      -------------
 *      |     *     |
 *      1     O     2
 *      |     *     |
 *      -------------
 *
 *      ^     ^     ^
 *      |     |     |
 *      |     |     neighboring cell center
 *      |     interpolated value at one "half cell step" in x direction
 *      center of origin cell
 * 
 *  In this case, there are 2 values required, one at the center (1) and at the neighbor (2)
 *  The reduced dimension of the projected interpolation domain is 1.
 *
 *  In the more general case, a reduced dimension of the projected interpolation domain has to be
 *  provided. This leads to 2 ^ DimProj values, e.g.
 *  DimProj    NUM_VALUES
 *    0           1         point interpolation, just the value at the center of the cell
 *    1           2         example above with two points on a line
 *    2           4         interpolation plane and therefore 4 points are required
 *    3           8         interpolation in a 3D box with 8 required points
 */
template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t DimProj>
[[nodiscard]] typename THelper<Dim, LO, GO SC>::IndexList
GetInterpolationIndices(const typename THelper<Dim, LO, GO SC>::Index& ind,
                        const typename THelper<Dim, LO, GO SC>::Options& off_rel,
                        const typename dare::utils::Vector<DimProj, std::size_t>& dim_aff) {
    using IndexList = THelper<Dim, DimProj>::IndexList;
    IndexList indices;
    indices.SetAllValues(ind);
    if constexpr (DimProj != 0) {
        for (std::size_t n_aff{0}; n_aff < DimProj; n_aff++) {
            std::size_t dim_aff_v{dim_aff[n_aff]};
            std::size_t n_step{math::pow(2, n_aff)};
            std::size_t n_range{n_step * 2};
            for (std::size_t n_v{0}; n_v < NUM_V; n_v += n_range)
                for (std::size_t pos{n_v}; pos < (pos + n_step); pos++)
                    indices[pos][dim_aff_v] += off_rel[dim_aff_v];
        }
    }
    return indices;
}

/*!
 * @brief interpolates a specified component at a half-cell step (e.g. center, face or node)
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 * @tparam SC scalar type
 * @tparam Dim dimension of grid
 * @tparam N number of components
 * @tparam DimProj dimension of interpolation domain
 * @param ind triplet of indices at the center
 * @param field reference to the field with the values
 * @param off_rel relative offset
 * @param dim_aff affected dimensions (where the relative offset in not zero)
 * @param n component id
 * @return interpolated scalar value
 */
template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N, std::size_t DimProj>
[[nodiscard]] SC InterpolateCartesianLinear(
    const typename THelper<Dim, LO, GO SC>::Index& ind,
    const typename THelper<Dim, LO, GO SC>::GridVector<N>& field,
    const typename THelper<Dim, LO, GO SC>::Options& off_rel,
    const typename dare::utils::Vector<DimProj, std::size_t>& dim_aff,
    std::size_t n) {
    using IndexList = typename THelper<Dim, DimProj>::IndexList;
    const std::size_t NUM_VALUES{THelper<Dim, DimProj>::NUM_VALUES};

    IndexList indices{GetInterpolationIndices(ind, off_rel, dim_aff)};
    SC value{field.At(indices[0], n)};
    for (std::size_t n_ind{1}; n_ind < NUM_VALUES; n_ind++) {
        value += field.At(indices[n_ind], n);
    }
    if constexpr (NUM_VALUES > 1) {
        value = math::Divide<NUM_VALUES>(value);
    }
    return value;
}

/*!
 * @brief interpolates all component at a half-cell step (e.g. center, face or node)
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 * @tparam SC scalar type
 * @tparam Dim dimension of grid
 * @tparam N number of components
 * @tparam DimProj dimension of interpolation domain
 * @param ind triplet of indices at the center
 * @param field reference to the field with the values
 * @param off_rel relative offset
 * @param dim_aff affected dimensions (where the relative offset in not zero)
 * @return vector of length N with interpolated scalar values
 */
template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N, std::size_t DimProj>
[[nodiscard]] dare::utils::Vector<N, SC>
InterpolateCartesianLinear(
    const typename THelper<Dim, LO, GO SC>::Index& ind,
    const typename THelper<Dim, LO, GO SC>::GridVector<N>& field,
    const typename THelper<Dim, LO, GO SC>::Options& off_rel,
    const typename dare::utils::Vector<DimProj, std::size_t>& dim_aff) {
    using IndexList = typename THelper<Dim, DimProj>::IndexList;
    const std::size_t NUM_VALUES{THelper<Dim, DimProj>::NUM_VALUES};

    IndexList indices{GetInterpolationIndices(ind, off_rel, dim_aff)};
    dare::utils::Vector<N, SC> values;
    for (std::size_t n_ind{0}; n_ind < NUM_VALUES; n_ind++) {
        values += field.GetValues(indices[n_ind]);
    }
    value *= Divisor<SC, NUM_VALUES>();
    return values;
}
}  // end namespace details::Cartesian

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC InterpolateToFace(const typename detail::Cartesian::THelper<Dim, LO, GO SC>::GridRepresentation& target,
                     const typename detail::Cartesian::THelper<Dim, LO, GO SC>::Index& ind_target,
                     const typename detail::Cartesian::THelper<Dim, LO, GO SC>::GridType::NeighborID face,
                     const typename detail::Cartesian::THelper<Dim, LO, GO SC>::GridVector<N>& field,
                     std::size_t n) {
    using CNB = dare::Grid::CartesianNeighbor;
    using THelper = detail::Cartesian::THelper<Dim, LO, GO SC>;
    using Options = typename THelper::Options;

    Options off_rel{target.GetOptions()};
    off_rel[ToFace(face) / 2 + 1] += ToNormal(face);
    Options off_field{field.GetGridRepresentation().GetOptions()};
    off_rel -= off_field;
    std::size_t n_dim_aff{0};
    for (auto e : off_rel)
        n_dim_aff += (e != 0);
    switch (n_dim_aff) {
    case 0:
        dare::utils::Vector<0, std::size_t> dim_aff;
        return detail::Cartesian::InterpolateCartesianLinear(ind_target, field, off_rel, dim_aff, n);
    case 1: {
        dare::utils::Vector<1, std::size_t> dim_aff;
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[0] = i;
            }
        }
        return detail::Cartesian::InterpolateCartesianLinear(ind_target, field, off_rel, dim_aff, n);
    }
    case 2: {
        dare::utils::Vector<2, std::size_t> dim_aff;
        std::size_t count{0};
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[count++] = i;
            }
        }
        return detail::Cartesian::InterpolateCartesianLinear(ind_target, field, off_rel, dim_aff, n);
    }
    case 3: {
        dare::utils::Vector<3, std::size_t> dim_aff;
        std::size_t count{0};
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[count++] = i;
            }
        }
        return detail::Cartesian::InterpolateCartesianLinear(ind_target, field, off_rel, dim_aff, n);
    }
    otherwise:
        std::cerr << "In " << __func__ << ": Interpolation not implemented for cases higher than 3D" << std::endl;
        return 0;
    }
}

}  // end namespace dare::math

#endif  // MATH_INTERPOLATION_CARTESIAN_H_
