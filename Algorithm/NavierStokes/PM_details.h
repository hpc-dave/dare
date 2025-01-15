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

#ifndef ALGORITHM_NAVIERSTOKES_PM_DETAILS_H_
#define ALGORITHM_NAVIERSTOKES_PM_DETAILS_H_

#include <concepts>
#include <set>

#include "Algorithm/AlgorithmTraits.h"
#include "Data/Field.h"
#include "PM_Information.h"


namespace dare::algorithm::detail {

template <typename TDict, ContainsViscosity TDictDefault>
struct pm_get_viscosity {
    using type = typename TDictDefault::viscosity;
};

template <ContainsViscosity TDict, ContainsViscosity TDictDefault>
struct pm_get_viscosity<TDict, TDictDefault> {
    using type = typename TDict::viscosity;
};

template <typename TDict, ContainsViscosity TDictDefault>
using pm_get_viscosity_t = pm_get_viscosity<TDict, TDictDefault>::type;

template <typename TDict, ContainsDensity TDictDefault>
struct pm_get_density{
    using type = typename TDictDefault::density;
};

template <ContainsDensity TDict, ContainsDensity TDictDefault>
struct pm_get_density<TDict, TDictDefault> {
    using type = typename TDict::density;
};

template <typename TDict, ContainsDensity TDictDefault>
using pm_get_density_t = pm_get_density<TDict, TDictDefault>::type;

template <typename TDict, ContainsPorosity TDictDefault>
struct pm_get_porosity {
    using type = typename TDictDefault::porosity;
};

template <ContainsPorosity TDict, ContainsPorosity TDictDefault>
struct pm_get_porosity<TDict, TDictDefault> {
    using type = typename TDict::porosity;
};

template <typename TDict, ContainsPorosity TDictDefault>
using pm_get_porosity_t = pm_get_porosity<TDict, TDictDefault>::type;

template <typename TDict, ContainsImplicitForce TDictDefault>
struct pm_get_implicit_force {
    using type = typename TDictDefault::implicit_force;
};

template <ContainsImplicitForce TDict, ContainsImplicitForce TDictDefault>
struct pm_get_implicit_force<TDict, TDictDefault> {
    using type = typename TDict::implicit_force;
};

template <typename TDict, ContainsImplicitForce TDictDefault>
using pm_get_implicit_force_t = pm_get_implicit_force<TDict, TDictDefault>::type;

template <typename TDict, ContainsExplicitForce TDictDefault>
struct pm_get_explicit_force {
    using type = typename TDictDefault::explicit_force;
};

template <ContainsExplicitForce TDict, ContainsExplicitForce TDictDefault>
struct pm_get_explicit_force<TDict, TDictDefault> {
    using type = typename TDict::explicit_force;
};

template <typename TDict, ContainsExplicitForce TDictDefault>
using pm_get_explicit_force_t = pm_get_explicit_force<TDict, TDictDefault>::type;

template <typename TDict, ContainsCompressible TDictDefault>
struct pm_get_compressible {
    static const bool value = TDictDefault::compressible::flag;
    using type = bool;
};

template <ContainsCompressible TDict, ContainsCompressible TDictDefault>
struct pm_get_compressible<TDict, TDictDefault> {
    using type = typename TDict::compressible;
    static const bool value = type::flag;
};

template <ContainsCompressible TDict, ContainsCompressible TDictDefault>
    requires std::is_same_v<decltype(TDict::compressible), const bool>
struct pm_get_compressible<TDict, TDictDefault> {
    using type = bool;
    static const bool value = TDict::compressible;
};

template <typename TDict, ContainsCompressible TDictDefault>
constexpr bool pm_get_compressible_v = pm_get_compressible<TDict, TDictDefault>::value;

template <typename TDict, ContainsCompressible TDictDefault>
using pm_get_compressible_t = pm_get_compressible<TDict, TDictDefault>::type;

template <typename TDict, ContainsTVD TDictDefault>
struct pm_get_tvd {
    using type = TDictDefault::tvd;
};

template <ContainsTVD TDict, ContainsTVD TDictDefault>
struct pm_get_tvd<TDict, TDictDefault> {
    using type = TDict::tvd;
};

template <typename TDict, ContainsTVD TDictDefault>
using pm_get_tvd_t = pm_get_tvd<TDict, TDictDefault>::type;

template <typename TDict, ContainsViscousStress TDictDefault>
struct pm_get_viscous_stress {
    using type = TDictDefault::viscous_stress;
};

template <ContainsViscousStress TDict, ContainsViscousStress TDictDefault>
struct pm_get_viscous_stress<TDict, TDictDefault> {
    using type = TDict::viscous_stress;
};

template <typename TDict, ContainsViscousStress TDictDefault>
using pm_get_viscous_stress_t = pm_get_viscous_stress<TDict, TDictDefault>::type;

template <typename TDict, ContainsMomentumIterations TDictDefault>
struct pm_get_momentum_iterations {
    using type = TDictDefault::momentum_iterations;
};

template <ContainsMomentumIterations TDict, ContainsMomentumIterations TDictDefault>
struct pm_get_momentum_iterations<TDict, TDictDefault> {
    using type = TDict::momentum_iterations;
};

template <typename TDict, ContainsMomentumIterations TDictDefault>
using pm_get_momentum_iterations_t = pm_get_momentum_iterations<TDict, TDictDefault>::type;

template <typename TDict, ContainsContinuityIterations TDictDefault>
struct pm_get_continuity_iterations {
    using type = TDictDefault::continuity_iterations;
};

template <ContainsContinuityIterations TDict, ContainsContinuityIterations TDictDefault>
struct pm_get_continuity_iterations<TDict, TDictDefault> {
    using type = TDict::continuity_iterations;
};

template <typename TDict, ContainsMomentumIterations TDictDefault>
using pm_get_continuity_iterations_t = pm_get_continuity_iterations<TDict, TDictDefault>::type;

template <typename TDict, ContainsTimeSchemeConvective TDictDefault>
struct pm_get_time_scheme_convective {
    using type = TDictDefault::time_scheme_convective;
};

template <ContainsTimeSchemeConvective TDict, ContainsTimeSchemeConvective TDictDefault>
struct pm_get_time_scheme_convective<TDict, TDictDefault> {
    using type = TDict::time_scheme_convective;
};

template <typename TDict, ContainsMomentumIterations TDictDefault>
using pm_get_time_scheme_convective_t = pm_get_continuity_iterations<TDict, TDictDefault>::type;

template <typename PropertyInfoUser, typename PropertyInfoDefault>
struct PMAssembledPropertyInfoWithDefaults {
    using _density_t = pm_get_density_t<PropertyInfoUser, PropertyInfoDefault>;
    using _viscosity_t = pm_get_viscosity_t<PropertyInfoUser, PropertyInfoDefault>;
    using _porosity_t = pm_get_porosity_t<PropertyInfoUser, PropertyInfoDefault>;
    using _implicit_force_t = pm_get_implicit_force_t<PropertyInfoUser, PropertyInfoDefault>;
    using _explicit_force_t = pm_get_explicit_force_t<PropertyInfoUser, PropertyInfoDefault>;
    static const bool _compressibility_v = pm_get_compressible_v<PropertyInfoUser, PropertyInfoDefault>;
    using density = utils::default_convert_to_tagged_info_t<_density_t, PMProperties, PMProperties::Density>;
    using viscosity = utils::default_convert_to_tagged_info_t<_viscosity_t, PMProperties, PMProperties::Viscosity>;
    using porosity = utils::default_convert_to_tagged_info_t<_porosity_t, PMProperties, PMProperties::Porosity>;
    using implicit_force = utils::default_convert_to_tagged_info_t<_implicit_force_t, PMProperties, PMProperties::ImplicitForce>;   // NOLINT
    using explicit_force = utils::default_convert_to_tagged_info_t<_explicit_force_t, PMProperties, PMProperties::ExplicitForce>;   // NOLINT
    using compressible = utils::FlaggedInfo<_compressibility_v, PMProperties, PMProperties::Compressible>;    // NOLINT
};

template <typename NumericalInfoUser, typename NumericalInfoDefault>
struct PMAssembledNumericalInfoWithDefaults {
    using tvd = pm_get_tvd_t<NumericalInfoUser, NumericalInfoDefault>;
    using viscous_stress = pm_get_viscous_stress_t<NumericalInfoUser, NumericalInfoDefault>;
    using momentum_iterations = pm_get_momentum_iterations_t<NumericalInfoUser, NumericalInfoDefault>;
    using continuity_iterations = pm_get_continuity_iterations_t<NumericalInfoUser, NumericalInfoDefault>;
    using time_scheme_convective = pm_get_time_scheme_convective_t<NumericalInfoUser, NumericalInfoDefault>;
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
    requires(requires { typename T::type; })  // NOLINT
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
    requires(T::NUM_COMPONENTS == 1)  // NOLINT
struct DetermineViscosityVariableType<T> {
    using type = const T*;
};

template <typename T>
    requires(requires { typename T::type; })  // NOLINT
using determine_viscosity_variable_type_t = typename DetermineViscosityVariableType<typename T::type>::type;

template <typename T>
struct DeterminePorosityVariableType {
};

template <dare::utils::NoneType T>
struct DeterminePorosityVariableType<T> {
    using type = T;
};

template <std::floating_point T>
struct DeterminePorosityVariableType<T> {
    using type = T;
};

template <std::integral T>
struct DeterminePorosityVariableType<T> {
    using type = double;
};

template <FieldType T>
    requires(T::NUM_COMPONENTS == 1)  // NOLINT
struct DeterminePorosityVariableType<T> {
    using type = const T*;
};

template <typename T>
    requires(requires { typename T::type; })  // NOLINT
using determine_porosity_variable_type_t = typename DeterminePorosityVariableType<typename T::type>::type;

template <typename T>
struct DetermineImplicitForceVariableType {
};

template <dare::utils::NoneType T>
struct DetermineImplicitForceVariableType<T> {
    using type = T;
    using member_type = T;
};

template <FieldType T>
    requires(T::NUM_COMPONENTS == 1)  // NOLINT
struct DetermineImplicitForceVariableType<T> {
    using type = const T*;
    using member_type = std::set<type>;
};

template <typename T>
    requires(requires { typename T::type; })  // NOLINT
using determine_implicit_force_variable_type_t = typename DetermineImplicitForceVariableType<typename T::type>::type;

template <typename T>
using determine_implicit_force_member_variable_type_t
    = typename DetermineImplicitForceVariableType<typename T::type>::member_type;

template <typename T>
struct DetermineExplicitForceVariableType {
};

template <dare::utils::NoneType T>
struct DetermineExplicitForceVariableType<T> {
    using type = T;
    using member_type = T;
};

template <FieldType T>
    requires(T::NUM_COMPONENTS == 1)  // NOLINT
struct DetermineExplicitForceVariableType<T> {
    using type = const T*;
    using member_type = std::set<type>;
};

template <typename T>
    requires(requires { typename T::type; })  // NOLINT
using determine_explicit_force_variable_type_t = typename DetermineExplicitForceVariableType<typename T::type>::type;

template <typename T>
    requires(requires { typename T::type; })  // NOLINT
using determine_explicit_force_member_variable_type_t
    = typename DetermineExplicitForceVariableType<typename T::type>::member_type;

}  // namespace dare::algorithm::detail

#endif  // ALGORITHM_NAVIERSTOKES_PM_DETAILS_H_
