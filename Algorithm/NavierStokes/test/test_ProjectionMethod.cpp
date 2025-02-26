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

#include <gtest/gtest.h>
#include <type_traits>
#include <set>
#include "Algorithm/NavierStokes/ProjectionMethod.h"
#include "Grid/Cartesian.h"

namespace dare::test {
template <typename T>
using DensityInfo = dare::PMDensityInfo<T>;

template<typename T>
using ViscosityInfo = dare::PMViscosityInfo<T>;

template <typename T>
using PorosityInfo = dare::PMImplicitForceInfo<T>;

template <typename T>
using ImplicitForceInfo = dare::PMImplicitForceInfo<T>;

template <typename T>
using ExplicitForceInfo = dare::PMExplicitForceInfo<T>;

template <bool Flag>
using CompressibleInfo = dare::PMCompressibleInfo<Flag>;

// a dummy for the boundary strategy
struct BStrat{
    template <typename T>
    void operator()(T t) {}
};

template <typename GridType, typename PDict, typename SDict>
using PM = dare::ProjectionMethod<GridType, BStrat, PDict, SDict>;

// that one needs to be here, as it cannot be declared locally
struct PDict_compressible_only_raw {
    static const bool compressible = true;
};

}  // namespace dare::test

TEST(ProjectionMethodTest, PropertyInfo) {
    // Static testing of the property propagation
    using GridType = dare::Cartesian<1>;
    using FieldType = dare::Field<GridType, typename GridType::ScalarType, 1>;
    using PDefault = dare::PMPropertyInfoDefault<GridType>;
    using SDefault = dare::PMNumericalInfoDefault;
    using PMProperties = dare::PMProperties;

    struct PDict_empty {
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::DensityInfo,
                  PDefault::density>);

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::ViscosityInfo,
                  PDefault::viscosity
                  >);

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::PorosityInfo,
                  PDefault::porosity>);

    // testing explicit force with default value None
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::ExplicitForceInfo,
                  PDefault::explicit_force>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::ExplicitForceInfo::type,
                  dare::None>);

    // testing implicit force with default value None
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::ImplicitForceInfo,
                  PDefault::implicit_force>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::ImplicitForceInfo::type,
                  dare::None>);

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_empty, SDefault>::CompressibilityInfo,
                  dare::FlaggedInfo<PDefault::compressible::flag, PMProperties, PMProperties::Compressible>>);

    using dare::test::DensityInfo;
    struct PDict_density_only_double_raw {
        using density = double;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_double_raw, SDefault>::DensityInfo,
                  dare::TaggedTypeInfo<double, PMProperties, PMProperties::Density>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_double_raw, SDefault>::DensityVariableType,
                  double>);

    struct PDict_density_only_double {
        using density = DensityInfo<double>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_double, SDefault>::DensityInfo,
                  DensityInfo<double>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_double, SDefault>::DensityVariableType,
                  double>);

    struct PDict_density_only_int_raw {
        using density = int;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_int_raw, SDefault>::DensityInfo,
                  dare::TaggedTypeInfo<int, PMProperties, PMProperties::Density>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_int_raw, SDefault>::DensityVariableType,
                  double>);

    struct PDict_density_only_int {
        using density = DensityInfo<int>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_int, SDefault>::DensityInfo,
                  DensityInfo<int>>);
    // integer is converted to int
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_int, SDefault>::DensityVariableType,
                  double>);

    struct PDict_density_only_field_raw {
        using density = FieldType;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_field_raw, SDefault>::DensityInfo,
                  dare::TaggedTypeInfo<FieldType, PMProperties, PMProperties::Density>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_field_raw, SDefault>::DensityVariableType,
                  const FieldType*>);

    struct PDict_density_only_field {
        using density = DensityInfo<FieldType>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_field, SDefault>::DensityInfo,
                  DensityInfo<FieldType>>);
    // the field is converted to a pointer
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_density_only_field, SDefault>::DensityVariableType,
                  const FieldType*>);

    using dare::test::ViscosityInfo;
    struct PDict_viscosity_only_double {
        using viscosity = ViscosityInfo<double>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_double, SDefault>::ViscosityInfo,
                  ViscosityInfo<double>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_double, SDefault>::ViscosityVariableType,    // NOLINT
                  double>);

    struct PDict_viscosity_only_double_raw {
        using viscosity = double;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_double_raw, SDefault>::ViscosityInfo,
                  dare::TaggedTypeInfo<double, PMProperties, PMProperties::Viscosity>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_double_raw, SDefault>::ViscosityVariableType,
                  double>);

    struct PDict_viscosity_only_int {
        using viscosity = ViscosityInfo<int>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_int, SDefault>::ViscosityInfo,
                  ViscosityInfo<int>>);
    // integer is converted to int
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_int, SDefault>::ViscosityVariableType,   // NOLINT
                  double>);

    struct PDict_viscosity_only_int_raw {
        using viscosity = int;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_int_raw, SDefault>::ViscosityInfo,
                  dare::TaggedTypeInfo<int, PMProperties, PMProperties::Viscosity>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_int_raw, SDefault>::ViscosityVariableType,
                  double>);

    struct PDict_viscosity_only_field {
        using viscosity = ViscosityInfo<FieldType>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_field, SDefault>::ViscosityInfo,
                  ViscosityInfo<FieldType>>);
    // the field is converted to a pointer
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_field, SDefault>::ViscosityVariableType, // NOLINT
                  const FieldType*>);

    struct PDict_viscosity_only_field_raw {
        using viscosity = FieldType;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_field_raw, SDefault>::ViscosityInfo,
                  dare::TaggedTypeInfo<FieldType, PMProperties, PMProperties::Viscosity>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_viscosity_only_field_raw, SDefault>::ViscosityVariableType,
                  const FieldType*>);

    using dare::test::PorosityInfo;
    struct PDict_porosity_only_none {
        using porosity = PorosityInfo<dare::None>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_none, SDefault>::PorosityInfo,
                  PorosityInfo<dare::None>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_none, SDefault>::PorosityVariableType,  // NOLINT
                  dare::None>);

    struct PDict_porosity_only_double {
        using porosity = PorosityInfo<double>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_double, SDefault>::PorosityInfo,
                  PorosityInfo<double>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_double, SDefault>::PorosityVariableType,  // NOLINT
                  double>);

    struct PDict_porosity_only_double_raw {
        using porosity = double;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_double_raw, SDefault>::PorosityInfo,
                  dare::TaggedTypeInfo<double, PMProperties, PMProperties::Porosity>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_double_raw, SDefault>::PorosityVariableType,
                  double>);

    struct PDict_porosity_only_int {
        using porosity = PorosityInfo<int>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_int, SDefault>::PorosityInfo,
                  PorosityInfo<int>>);
    // integer is converted to int
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_int, SDefault>::PorosityVariableType,  // NOLINT
                  double>);

    struct PDict_porosity_only_int_raw {
        using porosity = int;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_int_raw, SDefault>::PorosityInfo,
                  dare::TaggedTypeInfo<int, PMProperties, PMProperties::Porosity>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_int_raw, SDefault>::PorosityVariableType,
                  double>);

    struct PDict_porosity_only_field {
        using porosity = PorosityInfo<FieldType>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_field, SDefault>::PorosityInfo,
                  PorosityInfo<FieldType>>);
    // the field is converted to a pointer
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_field, SDefault>::PorosityVariableType,  // NOLINT
                  const FieldType*>);

    struct PDict_porosity_only_field_raw {
        using porosity = FieldType;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_field_raw, SDefault>::PorosityInfo,
                  dare::TaggedTypeInfo<FieldType, PMProperties, PMProperties::Porosity>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_porosity_only_field_raw, SDefault>::PorosityVariableType,
                  const FieldType*>);

    using dare::test::ImplicitForceInfo;
    struct PDict_imforce_only_none {
        using implicit_force = ImplicitForceInfo<dare::None>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_none, SDefault>::ImplicitForceInfo,
                  ImplicitForceInfo<dare::None>>);
    // the field is converted to a pointer
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_none, SDefault>::ImplicitForceVariableType,  // NOLINT
                  dare::None>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_none, SDefault>::ImplicitForceMemberType,  // NOLINT
                  dare::None>);

    struct PDict_imforce_only_field {
        using implicit_force = ImplicitForceInfo<FieldType>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_field, SDefault>::ImplicitForceInfo,
                  ImplicitForceInfo<FieldType>>);
    // the field is converted to a pointer
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_field, SDefault>::ImplicitForceVariableType,  // NOLINT
                  const FieldType*>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_field, SDefault>::ImplicitForceMemberType,  // NOLINT
                  std::set<const FieldType*>>);

    struct PDict_imforce_only_field_raw {
        using implicit_force = FieldType;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_field_raw, SDefault>::ImplicitForceInfo,
                  dare::TaggedTypeInfo<FieldType, PMProperties, PMProperties::ImplicitForce>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_field_raw, SDefault>::ImplicitForceVariableType,   // NOLINT
                  const FieldType*>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_imforce_only_field_raw, SDefault>::ImplicitForceMemberType,     // NOLINT
                  std::set<const FieldType*>>);

    using dare::test::ExplicitForceInfo;
    struct PDict_exforce_only_none {
        using explicit_force = dare::None;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_none, SDefault>::ExplicitForceInfo,
                  dare::TaggedTypeInfo<dare::None, PMProperties, PMProperties::ExplicitForce>>);
    // the field is converted to a pointer
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_none, SDefault>::ExplicitForceVariableType,  // NOLINT
                  dare::None>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_none, SDefault>::ExplicitForceMemberType,  // NOLINT
                  dare::None>);

    struct PDict_exforce_only_field {
        using explicit_force = ExplicitForceInfo<FieldType>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_field, SDefault>::ExplicitForceInfo,
                  ExplicitForceInfo<FieldType>>);
    // the field is converted to a pointer
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_field, SDefault>::ExplicitForceVariableType,  // NOLINT
                  const FieldType*>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_field, SDefault>::ExplicitForceMemberType,  // NOLINT
                  std::set<const FieldType*>>);

    struct PDict_exforce_only_field_raw {
        using explicit_force = FieldType;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_field_raw, SDefault>::ExplicitForceInfo,
                  dare::TaggedTypeInfo<FieldType, PMProperties, PMProperties::ExplicitForce>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_field_raw, SDefault>::ExplicitForceVariableType,  // NOLINT
                  const FieldType*>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_exforce_only_field_raw, SDefault>::ExplicitForceMemberType,  // NOLINT
                  std::set<const FieldType*>>);

    using dare::test::CompressibleInfo;
    struct PDict_compressible_only {
        using compressible = CompressibleInfo<true>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDict_compressible_only, SDefault>::CompressibilityInfo,
                  dare::FlaggedInfo<true, PMProperties, PMProperties::Compressible>>);
    static_assert(dare::test::PM<GridType, PDict_compressible_only, SDefault>::compressible);

    static_assert(std::is_same_v<
                  dare::test::PM<GridType,dare::test::PDict_compressible_only_raw, SDefault>::CompressibilityInfo,    // NOLINT
                  dare::FlaggedInfo<true, PMProperties, PMProperties::Compressible>>);
    static_assert(dare::test::PM<GridType, dare::test::PDict_compressible_only_raw, SDefault>::compressible);         // NOLINT
}

TEST(ProjectionMethodTest, NumericalInfo) {
    // Static testing of the property propagation
    using GridType = dare::Cartesian<1>;
    using PDefault = dare::PMPropertyInfoDefault<GridType>;
    using SDefault = dare::PMNumericalInfoDefault;

    struct SDict_empty {
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::TVDInfo,
                  SDefault::tvd>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::TVDScheme,
                  SDefault::tvd::type>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::ViscousStressInfo,
                  SDefault::viscous_stress>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::ViscousStressTreatment,
                  SDefault::viscous_stress::type>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::MomentumIterationInfo,
                  SDefault::momentum_iterations>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::MomentumIterationType,
                  SDefault::momentum_iterations::type>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::ContinuityIterationInfo,
                  SDefault::continuity_iterations>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::ContinuityIterationType,
                  SDefault::continuity_iterations::type>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::ConvectiveTimeSchemeInfo,
                  SDefault::time_scheme_convective>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_empty>::ConvectiveTimeSchemeType,
                  SDefault::time_scheme_convective::type>);

    struct SDict_tvd_info_only {
        using tvd = dare::PMTVDInfo<dare::VANALBADA>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_tvd_info_only>::TVDInfo,
                  dare::PMTVDInfo<dare::VANALBADA>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_tvd_info_only>::TVDScheme,
                  dare::VANALBADA>);

    struct SDict_tvd_info_only_raw {
        using tvd = dare::VANALBADA;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_tvd_info_only_raw>::TVDInfo,
                  dare::TaggedTypeInfo<dare::VANALBADA, dare::PMNumerics, dare::PMNumerics::TVD>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_tvd_info_only_raw>::TVDScheme,
                  dare::VANALBADA>);

    struct SDict_momentum_iteration_info_only {
        using momentum_iterations = dare::PMMomentumIterationInfo<dare::Newton>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_momentum_iteration_info_only>::MomentumIterationInfo,
                  dare::PMMomentumIterationInfo<dare::Newton>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_momentum_iteration_info_only>::MomentumIterationType,
                  dare::Newton>);

    struct SDict_momentum_iteration_info_only_raw {
        using momentum_iterations = dare::Newton;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_momentum_iteration_info_only_raw>::MomentumIterationInfo,
                  dare::TaggedTypeInfo<dare::Newton, dare::PMNumerics, dare::PMNumerics::MomentumIterations>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_momentum_iteration_info_only_raw>::MomentumIterationType,
                  dare::Newton>);

    struct SDict_continuity_iteration_info_only {
        using continuity_iterations = dare::PMContinuityIterationInfo<dare::Newton>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_continuity_iteration_info_only>::ContinuityIterationInfo,
                  dare::PMContinuityIterationInfo<dare::Newton>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_continuity_iteration_info_only>::ContinuityIterationType,
                  dare::Newton>);

    struct SDict_continuity_iteration_info_only_raw {
        using continuity_iterations = dare::Newton;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_continuity_iteration_info_only_raw>::ContinuityIterationInfo,
                  dare::TaggedTypeInfo<dare::Newton, dare::PMNumerics, dare::PMNumerics::ContinuityIterations>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_continuity_iteration_info_only_raw>::ContinuityIterationType,
                  dare::Newton>);

    struct SDict_ts_convective_info_only {
        using time_scheme_convective = dare::PMTimeSchemeConvectiveInfo<dare::EULER_FORWARD>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_ts_convective_info_only>::ConvectiveTimeSchemeInfo,
                  dare::PMTimeSchemeConvectiveInfo<dare::EULER_FORWARD>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_ts_convective_info_only>::ConvectiveTimeSchemeType,
                  dare::EULER_FORWARD>);

    struct SDict_ts_convective_info_only_raw {
        using time_scheme_convective = dare::EULER_FORWARD;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_ts_convective_info_only_raw>::ConvectiveTimeSchemeInfo,
                  dare::TaggedTypeInfo<dare::EULER_FORWARD, dare::PMNumerics, dare::PMNumerics::TimeSchemeConvective>>);    // NOLINT
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_ts_convective_info_only_raw>::ConvectiveTimeSchemeType,
                  dare::EULER_FORWARD>);

    struct SDict_vstress_info_only {
        using viscous_stress = dare::PMViscousStressInfo<dare::PMDijkhuizenStressTensor>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_vstress_info_only>::ViscousStressInfo,
                  dare::PMViscousStressInfo<dare::PMDijkhuizenStressTensor>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_vstress_info_only>::ViscousStressTreatment,
                  dare::PMDijkhuizenStressTensor>);

    struct SDict_vstress_info_only_raw {
        using viscous_stress = dare::PMDijkhuizenStressTensor;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_vstress_info_only_raw>::ViscousStressInfo,
                  dare::TaggedTypeInfo<dare::PMDijkhuizenStressTensor, dare::PMNumerics, dare::PMNumerics::ViscousStress>>);  // NOLINT
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_vstress_info_only_raw>::ViscousStressTreatment,
                  dare::PMDijkhuizenStressTensor>);

    struct SDict_mom_normalizer_info_only {
        using momentum_normalizer = dare::PMMomentumNormalizerInfo<int>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_mom_normalizer_info_only>::MomentumNormalizerInfo,
                  dare::PMMomentumNormalizerInfo<int>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_mom_normalizer_info_only>::MomentumNormalizerType,
                  int>);

    struct SDict_mom_normalizer_info_only_raw {
        using momentum_normalizer = int;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_mom_normalizer_info_only_raw>::MomentumNormalizerInfo,
                  dare::TaggedTypeInfo<int, dare::PMNumerics, dare::PMNumerics::MomentumNormalizer>>);  // NOLINT
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_mom_normalizer_info_only_raw>::MomentumNormalizerType,
                  int>);

    struct SDict_cont_normalizer_info_only {
        using continuity_normalizer = dare::PMContinuityNormalizerInfo<int>;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_cont_normalizer_info_only>::ContinuityNormalizerInfo,
                  dare::PMContinuityNormalizerInfo<int>>);
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_cont_normalizer_info_only>::ContinuityNormalizerType,
                  int>);

    struct SDict_cont_normalizer_info_only_raw {
        using continuity_normalizer = int;
    };

    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_cont_normalizer_info_only_raw>::ContinuityNormalizerInfo,
                  dare::TaggedTypeInfo<int, dare::PMNumerics, dare::PMNumerics::ContinuityNormalizer>>);  // NOLINT
    static_assert(std::is_same_v<
                  dare::test::PM<GridType, PDefault, SDict_cont_normalizer_info_only_raw>::ContinuityNormalizerType,
                  int>);
}
