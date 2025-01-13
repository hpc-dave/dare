/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder

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

#ifndef ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_
#define ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_

#include <concepts>
#include "PM_Information.h"
#include "Algorithm/AlgorithmTraits.h"
#include "Data/Field.h"

namespace dare::algorithm {

/*!
 * \brief tagging class for Dijkhuizens approach for the stress Tensor
 */
struct PMDijkhuizenStressTensor {
};

/*!
 * \brief tagging class for default approach for the stress Tensor
 */
struct PMDefaultStressTensor {
};

template <typename T>
concept PMDijkhuizenStressTreatment =
    std::is_base_of_v<PMDijkhuizenStressTensor, T>;

template <typename T>
concept PMDefaultStressTreatment =
    std::is_base_of_v<PMDefaultStressTensor, T>;

template <typename T>
struct is_pm_dijkuizen_stress_tensor : std::bool_constant<PMDijkhuizenStressTreatment<T>> {};

template <typename T>
struct is_pm_default_stress_tensor : std::bool_constant<PMDefaultStressTreatment<T>> {};

template<typename Grid>
struct PMPropertyInfoDefault {
    using FieldType = dare::Data::Field<Grid, typename Grid::ScalarType, 1>;
    using density = PMDensityInfo<FieldType>;
    using viscosity = PMViscosityInfo<FieldType>;
    using porosity = PMPorosityInfo<dare::utils::None>;
    using implicit_force = PMImplicitForceInfo<dare::utils::None>;
    using explicit_force = PMExplicitForceInfo<dare::utils::None>;
    using compressible = PMCompressibleInfo<false>;
};

template <typename Grid>
struct PMSolverInfoDefault {
    using MomentumAlgorithm = dare::algorithm::FixedPoint;
    using ContinuityAlgorithm = dare::algorithm::Newton;
};

namespace detail {

template<typename PropertyInfoUser, typename PropertyInfoDefault>
struct PMAssembledPropertyInfoWithDefaults {
    using density = decltype(PMGetDensity<PropertyInfoUser, PropertyInfoDefault>());
    using viscosity = decltype(PMGetViscosity<PropertyInfoUser, PropertyInfoDefault>());
    using porosity = decltype(PMGetPorosity<PropertyInfoUser, PropertyInfoDefault>());
    using implicit_force = decltype(PMGetImplicitForce<PropertyInfoUser, PropertyInfoDefault>());
    using explicit_force = decltype(PMGetExplicitForce<PropertyInfoUser, PropertyInfoDefault>());
    using compressible = decltype(PMGetCompressibility<PropertyInfoUser, PropertyInfoDefault>());
};

template <typename T>
struct DetermineDensityVariableType {
};

template <std::floating_point T>
struct DetermineDensityVariableType<T> {
    using type = T;
};

template <std::integral T>
struct DetermineDensityVariableType<T> {
    using type = double;
};

template <FieldType T>
    requires(T::NUM_COMPONENTS == 1)  // NOLINT
struct DetermineDensityVariableType<T> {
    using type = const T*;
};

template <typename T>
requires( requires {typename T::type;} )    // NOLINT
using determine_density_variable_type_t = typename DetermineDensityVariableType<typename T::type>::type;

template <typename T>
struct DetermineViscosityVariableType {
};

template <std::floating_point T>
struct DetermineViscosityVariableType<T> {
    using type = T;
};

template <std::integral T>
struct DetermineViscosityVariableType<T> {
    using type = double;
};

template <FieldType T>
requires (T::NUM_COMPONENTS == 1)       // NOLINT
struct DetermineViscosityVariableType<T> {
    using type = const T*;
};

template <typename T>
requires(requires { typename T::type; })    // NOLINT
using determine_viscosity_variable_type_t = typename DetermineViscosityVariableType<typename T::type>::type;

}  // namespace detail

/*!
 * @brief holds all relevant entities to compute the flow field via a two-step projection method
 * @tparam Grid the grid type
 * @tparam AlgorithmInfo compile time information about the algorithm
 * @tparam PropertyInfo compile time information about the properties
 *
 * EXPLAIN THE COMPILE TIME INFORMATION APPROACH
 * Compile time information concerning the properties:
 * - density:
 * - viscosity:
 * - porosity:
 * - implicit_force:
 * - explicit_force:
 *
 *
 */
template <typename Grid, typename PropertyInfo, typename AlgorithmInfo>
class ProjectionMethod {
public:
    // general types based on the grid
    using GridType = Grid;
    using SC = typename Grid::ScalarType;
    using LO = typename Grid::LocalOrdinalType;
    using GO = typename Grid::GlobalOrdinalType;
    using Index = typename Grid::Index;
    using IndexGlobal = typename Grid::IndexGlobal;
    using Field = Data::Field<GridType, SC, 1>;

    // properties determined from the PropertyInfo type
    using PropertyTypeInfo = detail::PMAssembledPropertyInfoWithDefaults<PropertyInfo, PMPropertyInfoDefault<Grid>>;
    using DensityInfo = typename PropertyTypeInfo::density;
    using ViscosityInfo = typename PropertyTypeInfo::viscosity;
    using PorosityInfo = typename PropertyTypeInfo::porosity;
    using ImplicitForceInfo = typename PropertyTypeInfo::implicit_force;
    using ExplicitForceInfo = typename PropertyTypeInfo::explicit_force;
    using CompressibilityInfo = typename PropertyTypeInfo::compressible;
    static const bool compressible = CompressibilityInfo::flag;
    static const std::size_t dimension = Grid::Dimension;
    using DensityVariableType = detail::determine_density_variable_type_t<DensityInfo>;
    using ViscosityVariableType = detail::determine_viscosity_variable_type_t<ViscosityInfo>;

    constexpr bool IsCompressible() const { return compressible; }
    constexpr std::size_t GetDimension() const { return dimension; }



private:
    const Field* rho;
    const Field* mu;
};

}  // namespace dare::algorithm

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_H_
