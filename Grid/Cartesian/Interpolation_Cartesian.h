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

#ifndef GRID_CARTESIAN_INTERPOLATION_CARTESIAN_H_
#define GRID_CARTESIAN_INTERPOLATION_CARTESIAN_H_

#include <limits>

#include "Math/Interpolation.h"
#include "Math/Divisors.h"
#include "Grid/Cartesian.h"
#include "Utilities/Errors.h"

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
    static const std::size_t NUM_VALUES{math::Pow<2, DimProj>()};
    using GridType = dare::Grid::Cartesian<Dim>;
    using Index = typename GridType::Index;
    using IndexList = dare::utils::Vector<NUM_VALUES, Index>;
    using Options = typename GridType::Options;
    template <typename SC, std::size_t N>
    using GridVector = dare::Data::GridVector<GridType, SC, N>;
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
 * @return List of indices for interpolation, where the first index is the one furthest removed
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
 * 
 * \note A special requirement for later use is, that the first value equals the origin cell
 * and the last value corresponds to the cell which is furthes away!
 */
template <std::size_t Dim, std::size_t DimProj>
[[nodiscard]] typename THelper<Dim, DimProj>::IndexList
GetInterpolationIndices(const dare::utils::Vector<Dim, dare::defaults::LocalOrdinalType>& ind,
                        const typename THelper<Dim, DimProj>::Options& off_rel,
                        const typename dare::utils::Vector<DimProj, std::size_t>& dim_aff) {
    using IndexList = typename THelper<Dim, DimProj>::IndexList;
    static const std::size_t NUM_V{THelper<Dim, DimProj>::NUM_VALUES};
    const std::size_t TWO{2};   // To avoid narrowing conversion to int further down
    IndexList indices;
    indices.SetAllValues(ind);
    if constexpr (DimProj != 0) {
        for (std::size_t n_aff{0}; n_aff < DimProj; n_aff++) {
            std::size_t dim_aff_v{dim_aff[n_aff]};
            std::size_t n_step{math::Pow(TWO, n_aff)};
            std::size_t n_range{n_step * TWO};
            for (std::size_t n_v{0}; n_v < NUM_V; n_v += n_range)
                for (std::size_t pos{n_v}; pos < (n_v + n_step); pos++)
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
template <std::size_t Dim, typename SC, std::size_t N, std::size_t DimProj>
[[nodiscard]] SC InterpolateCartesianLinear(
    const dare::utils::Vector<Dim, dare::defaults::LocalOrdinalType>& ind,
    const dare::Data::GridVector<Grid::Cartesian<Dim>, SC, N>& field,
    const typename THelper<Dim, DimProj>::Options& off_rel,
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
template <std::size_t Dim, typename SC, std::size_t N, std::size_t DimProj>
[[nodiscard]] dare::utils::Vector<N, SC>
InterpolateCartesianLinear(
    const dare::utils::Vector<Dim, dare::defaults::LocalOrdinalType>& ind,
    const dare::Data::GridVector<Grid::Cartesian<Dim>, SC, N>& field,
    const typename THelper<Dim, DimProj>::Options& off_rel,
    const typename dare::utils::Vector<DimProj, std::size_t>& dim_aff) {
    using IndexList = typename THelper<Dim, DimProj>::IndexList;
    const std::size_t NUM_VALUES{THelper<Dim, DimProj>::NUM_VALUES};

    IndexList indices{GetInterpolationIndices(ind, off_rel, dim_aff)};
    dare::utils::Vector<N, SC> values;
    for (std::size_t n_ind{0}; n_ind < NUM_VALUES; n_ind++) {
        values += field.GetValues(indices[n_ind]);
    }
    values *= Divisor<SC, NUM_VALUES>();
    return values;
}

/*!
 * @brief computes the relative offset for interpolation depending on the grids and face requested
 * @tparam Dim dimension of grid
 * @tparam LO type of local ordinal
 * @param opt_source information concerning the staggering of the source grid
 * @param opt_target information concerning the staggering of the target grid
 * @param face identifier of the target face
 * @param ind indices of the center cell, which may be adapted
 * @param off_rel relative offset will be set
 */
template <std::size_t Dim, typename LO>
void GetRelativeOffset(typename utils::Vector<Dim, LO> opt_source,
                       typename Grid::Cartesian<Dim>::Options opt_target,
                       typename Grid::Cartesian<Dim>::NeighborID face,
                       typename Grid::Cartesian<Dim>::Index* ind,
                       typename Grid::Cartesian<Dim>::Options* off_rel) {
    opt_source[ToFace(face) / 2] += ToNormal(face);  // accounting for the position of the face

    (*off_rel) = opt_source - opt_target;  // relative offset of the fields

    // Here, we correct for the case, that we should go 2 "half-steps" or more
    for (std::size_t dim{0}; dim < Dim; dim++) {
        LO v_corr = (*off_rel)[dim] / 2;  // 2 "half-steps" means a whole cell
        (*ind)[dim] += v_corr;
        (*off_rel)[dim] = (*off_rel)[dim] % 2;
    }
}

template <std::size_t Dim, typename SC>
dare::utils::Vector<THelper<Dim, Dim>::NUM_VALUES, SC> GetLinearInterpolationWeights(
    const typename dare::Grid::Cartesian<Dim>::Representation& grid,
    const dare::utils::Vector<Dim, SC>& poi,
    const dare::utils::Vector<Dim, SC>& point_center,
    const dare::utils::Vector<Dim, SC>& point_far,
    const dare::utils::Vector<THelper<Dim, Dim>::NUM_VALUES, typename dare::Grid::Cartesian<Dim>::Index>& indices) {
    const std::size_t NUM_VALUES{details::Cartesian::THelper<Dim, Dim>::NUM_VALUES};
    using GridType = dare::Grid::Cartesian<Dim>;
    using WList = dare::utils::Vector<NUM_VALUES, SC>;
    using Index = typename GridType::Index;

#ifndef DARE_NDEBUG
    // check for requirements on the indices
    Index ind_center = grid.GetCell(point_center);
    Index ind_far = grid.GetCell(point_far);
    if (indices[0] != ind_far) {
        ERROR << "Indices are not sorted as expected! "
              << "The first index needs to correspond to the one furthest away from the cell"
              << ERROR_CLOSE;
    }
    if (indices[NUM_VALUES - 1] != ind_center) {
        ERROR << "Indices are not sorted as expected! "
              << "The last index needs to correspond to the one of the center cell"
              << ERROR_CLOSE;
    }
#endif

    // define some relevant values
    WList weights;
    weights.SetAllValues(1);
    dare::utils::Vector<Dim, SC> delta_close{std::abs(point_center - poi)};
    dare::utils::Vector<Dim, SC> delta_far{std::abs(point_far - poi)};
    const Index& ind{indices[NUM_VALUES - 1]};

    // assign weights to each corner, depending on their position
    for (std::size_t n{0}; n < NUM_VALUES; n++) {
        for (std::size_t d{0}; d < Dim; d++) {
            SC w = (indices[n][d] == ind[d]) * delta_far[d];
            w += (indices[n][d] != ind[d]) * delta_close[d];
            weights[n] *= w;
        }
    }

    // normalization of the weights
    SC sum_weights{0};
    for (auto e : weights)
        sum_weights += e;
    weights /= sum_weights;

    return weights;
}

}  // end namespace details::Cartesian

/*!
 * @brief interpolates one component to the face of the target from the source field
 * @tparam SC type of scalar to interpolate
 * @tparam Dim dimension of the grid
 * @tparam N number of components
 * @param target grid of the target
 * @param ind_target triplet of the target cell
 * @param face face of target cell, at which the value is required
 * @param field reference to the field, from which the value is interpolated
 * @return interpolated value the face
 * This is an optimized interpolation function, in the case that a value is required at
 * a face (or center) of a cell of a Cartesian grid. The grids may also be the same.
 * An example for its application:
 * @code{.cpp}
 * // assuming that GridRepresentation g_rep of our target and a field are defined
 * Index ind{3, 4, 5};  // Cartesian Index in 3D
 * dare::Grid::CartesianNeighbor face = dare::Grid::CartesianNeighbor::WEST;
 *
 * // Get the value of component 0 at the west face
 * SC value = InterpolateToFace(g_rep, ind, face, field, 0);
 * @endcode
 */
template <std::size_t Dim, typename SC, std::size_t N>
[[nodiscard]] SC InterpolateToFace(const typename Grid::Cartesian<Dim>::Representation& target,
                                   const typename Grid::Cartesian<Dim>::Index& ind_target,
                                   const typename Grid::Cartesian<Dim>::NeighborID face,
                                   const typename Data::GridVector<Grid::Cartesian<Dim>, SC, N>& field,
                                   std::size_t n) {
    using Options = typename dare::Grid::Cartesian<Dim>::Options;

    Options stagg_source{field.GetGridRepresentation().GetOptions()};     // staggered position of source field
    Options stagg_target{target.GetOptions()};                            // staggered position of target

    Options off_rel;
    typename Grid::Cartesian<Dim>::Index ind_corr{ind_target};

    details::Cartesian::GetRelativeOffset(stagg_source, stagg_target, face, &ind_corr, &off_rel);

    std::size_t n_dim_aff{0};
    for (auto e : off_rel)
        n_dim_aff += (e != 0);
    switch (n_dim_aff) {
    case 0: {
        dare::utils::Vector<0, std::size_t> dim_aff;
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff, n);
    } break;
    case 1: {
        dare::utils::Vector<1, std::size_t> dim_aff;
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[0] = i;
            }
        }
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff, n);
    } break;
    case 2: {
        dare::utils::Vector<2, std::size_t> dim_aff;
        std::size_t count{0};
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[count++] = i;
            }
        }
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff, n);
    } break;
    case 3: {
        dare::utils::Vector<3, std::size_t> dim_aff;
        std::size_t count{0};
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[count++] = i;
            }
        }
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff, n);
    } break;
    }
    ERROR << "Interpolation not implemented for cases higher than 3D" << ERROR_CLOSE;
    return std::numeric_limits<SC>::signaling_NaN();
}

/*!
 * @brief interolates all components to the face of the target from the source field
 * @tparam SC type of scalar to interpolate
 * @tparam Dim dimension of the grid
 * @tparam N number of components
 * @param target grid of the target
 * @param ind_target triplet of the target cell
 * @param face face of target cell, at which the value is required
 * @param field reference to the field, from which the value is interpolated
 * @return vector with all components interpolated to the face
 * This is an optimized interpolation function, in the case that a value is required at
 * a face (or center) of a cell of a Cartesian grid. The grids may also be the same.
 * An example for its application:
 * @code{.cpp}
 * // assuming that GridRepresentation g_rep of our target and a field are defined
 * Index ind{3, 4, 5};  // Cartesian Index in 3D
 * dare::Grid::CartesianNeighbor face = dare::Grid::CartesianNeighbor::WEST;
 *
 * // Get the value of all component at the west face
 * dare::utils::Vector<Dim, SC> value = InterpolateToFace(g_rep, ind, face, field);
 * @endcode
 */
template <std::size_t Dim, typename SC, std::size_t N>
dare::utils::Vector<N, SC> InterpolateToFace(const typename Grid::Cartesian<Dim>::Representation& target,
                                   const typename Grid::Cartesian<Dim>::Index& ind_target,
                                   const typename Grid::Cartesian<Dim>::NeighborID face,
                                   const typename Data::GridVector<Grid::Cartesian<Dim>, SC, N>& field) {
    using Options = typename dare::Grid::Cartesian<Dim>::Options;

    Options stagg_source{field.GetGridRepresentation().GetOptions()};     // staggered position of source field
    Options stagg_target{target.GetOptions()};                            // staggered position of target

    Options off_rel;
    typename Grid::Cartesian<Dim>::Index ind_corr{ind_target};

    details::Cartesian::GetRelativeOffset(stagg_source, stagg_target, face, &ind_corr, &off_rel);

    std::size_t n_dim_aff{0};
    for (auto e : off_rel)
        n_dim_aff += (e != 0);
    switch (n_dim_aff) {
    case 0: {
        dare::utils::Vector<0, std::size_t> dim_aff;
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff);
    } break;
    case 1: {
        dare::utils::Vector<1, std::size_t> dim_aff;
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[0] = i;
            }
        }
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff);
    } break;
    case 2: {
        dare::utils::Vector<2, std::size_t> dim_aff;
        std::size_t count{0};
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[count++] = i;
            }
        }
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff);
    } break;
    case 3: {
        dare::utils::Vector<3, std::size_t> dim_aff;
        std::size_t count{0};
        for (std::size_t i{0}; i < Dim; i++) {
            if (off_rel[i] != 0) {
                dim_aff[count++] = i;
            }
        }
        return details::Cartesian::InterpolateCartesianLinear(ind_corr, field, off_rel, dim_aff);
    } break;
    }
    ERROR << "Interpolation not implemented for cases higher than 3D" << ERROR_CLOSE;
    return dare::utils::Vector<N, SC>(std::numeric_limits<SC>::signaling_NaN());
}

template <std::size_t Dim, typename SC, std::size_t N>
SC InterpolateToPoint(const typename dare::utils::Vector<Dim, SC>& point,
                      const Data::GridVector<dare::Grid::Cartesian<Dim>, SC, N>& field,
                      std::size_t n) {
    const std::size_t NUM_VALUES{details::Cartesian::THelper<Dim, Dim>::NUM_VALUES};
    using GridType = dare::Grid::Cartesian<Dim>;
    using Index = typename GridType::Index;
    using LO = typename GridType::LocalOrdinalType;
    using Point = dare::utils::Vector<Dim, SC>;
    using VecList = dare::utils::Vector<NUM_VALUES, Index>;
    using WList = dare::utils::Vector<NUM_VALUES, SC>;

    Index ind = field.GetGridRepresentation().GetCell(point);
    Point point_center = field.GetGridRepresentation().GetCoordinatesCenter(ind);
    dare::utils::Vector<Dim, LO> off_rel;
    Point point_rel{point - point_center};
    for (std::size_t d{0}; d < Dim; d++) {
        off_rel[d] = static_cast<LO>(point_rel[d] >= 0) - (point_rel[d] < 0);
    }
    dare::utils::Vector<Dim, std::size_t> dim_aff;
    for (std::size_t d{0}; d < Dim; d++)
        dim_aff[d] = d;

    VecList vlist = details::Cartesian::GetInterpolationIndices(
        ind, off_rel, dim_aff);

    Point point_far = field.GetGridRepresentation().GetCoordinatesCenter(vlist[0]);
    WList weights = details::Cartesian::GetLinearInterpolationWeights(field.GetGridRepresentation(),
                                                                      point,
                                                                      point_center,
                                                                      point_far,
                                                                      vlist);

    return Interpolate(field, vlist, weights, n);
}

template <std::size_t Dim, typename SC, std::size_t N>
dare::utils::Vector<N, SC> InterpolateToPoint(const typename dare::utils::Vector<Dim, SC>& point,
                                              const Data::GridVector<dare::Grid::Cartesian<Dim>, SC, N>& field) {
    const std::size_t NUM_VALUES{details::Cartesian::THelper<Dim, Dim>::NUM_VALUES};
    using GridType = dare::Grid::Cartesian<Dim>;
    using Index = typename GridType::Index;
    using LO = typename GridType::LocalOrdinalType;
    using Point = dare::utils::Vector<Dim, SC>;
    using VecList = dare::utils::Vector<NUM_VALUES, Index>;
    using WList = dare::utils::Vector<NUM_VALUES, SC>;

    Index ind = field.GetGridRepresentation().GetCell(point);
    Point point_center = field.GetGridRepresentation().GetCoordinatesCenter(ind);
    dare::utils::Vector<Dim, LO> off_rel;
    Point point_rel{point - point_center};
    for (std::size_t d{0}; d < Dim; d++) {
        off_rel[d] = static_cast<LO>(point_rel[d] >= 0) - (point_rel[d] < 0);
    }
    dare::utils::Vector<Dim, std::size_t> dim_aff;
    for (std::size_t d{0}; d < Dim; d++)
        dim_aff[d] = d;

    VecList vlist = details::Cartesian::GetInterpolationIndices(
        ind, off_rel, dim_aff);

    Point point_far = field.GetGridRepresentation().GetCoordinatesCenter(vlist[0]);
    WList weights = details::Cartesian::GetLinearInterpolationWeights(field.GetGridRepresentation(),
                                                                      point,
                                                                      point_center,
                                                                      point_far,
                                                                      vlist);

    return Interpolate(field, vlist, weights);
}
}  // end namespace dare::math

#endif  // GRID_CARTESIAN_INTERPOLATION_CARTESIAN_H_
