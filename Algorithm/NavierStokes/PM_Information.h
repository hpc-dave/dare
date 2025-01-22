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

#ifndef ALGORITHM_NAVIERSTOKES_PM_INFORMATION_H_
#define ALGORITHM_NAVIERSTOKES_PM_INFORMATION_H_

#include "Equations/TimeDiscretizationSchemes.h"
#include "Utilities/PropertyInformation.h"
#include "Algorithm/StaticInformation.h"

namespace dare {

/*!
 * \brief information related to the properties of the system simulated via the projection method
 */
enum class PMProperties {
    Density,
    Viscosity,
    Porosity,
    ImplicitForce,
    ExplicitForce,
    Compressible
};

/*!
 * \brief information related to the numerical approach and discretization of the projection method
 */
enum class PMNumerics {
    TVD,
    ViscousStress,
    TimeSchemeConvective,
    MomentumIterations,
    ContinuityIterations
};

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

template <bool Flag, PMProperties Tag>
using PMFlaggedInfo = FlaggedInfo<Flag, PMProperties, Tag>;

template <typename Type, PMProperties Tag>
using PMTypeInfo = TaggedTypeInfo<Type, PMProperties, Tag>;

template <typename Type, PMProperties Tag, std::size_t N>
using PMTypeCountedInfo = TaggedCountedTypeInfo<Type, PMProperties, Tag, N>;

template <typename Type, PMNumerics Tag>
using PMNumericsInfo = TaggedTypeInfo<Type, PMNumerics, Tag>;

template <typename Type>
struct PMDensityInfo : PMTypeInfo<Type, PMProperties::Density> {};

template <typename Type>
struct PMViscosityInfo : PMTypeInfo<Type, PMProperties::Viscosity> {};

template <typename Type>
struct PMPorosityInfo : PMTypeInfo<Type, PMProperties::Porosity> {};

template <typename Type>
struct PMImplicitForceInfo : PMTypeInfo<Type, PMProperties::ImplicitForce> {};

template <typename Type>
struct PMExplicitForceInfo : PMTypeInfo<Type, PMProperties::ExplicitForce> {};

template <bool Flag>
struct PMCompressibleInfo : PMFlaggedInfo<Flag, PMProperties::Compressible> {};

template <typename Type>
struct PMTVDInfo : PMNumericsInfo<Type, PMNumerics::TVD> {};

template <typename Type>
struct PMViscousStressInfo : PMNumericsInfo<Type, PMNumerics::ViscousStress> {};

template <typename Type>
struct PMMomentumIterationInfo : PMNumericsInfo<Type, PMNumerics::MomentumIterations> {};

template <typename Type>
struct PMContinuityIterationInfo : PMNumericsInfo<Type, PMNumerics::ContinuityIterations> {};

template <typename Type>
struct PMTimeSchemeConvectiveInfo : PMNumericsInfo<Type, PMNumerics::TimeSchemeConvective> {};

/*!
 * \brief concept for determining Dijkhuizens stress treatment
 * @tparam T type to test
 */
template <typename T>
concept PMDijkhuizenStressTreatment =
    std::is_base_of_v<PMDijkhuizenStressTensor, T>;

/*!
 * \brief concept for determining the default stress treatment
 * @tparam T type to test
 */
template <typename T>
concept PMDefaultStressTreatment =
    std::is_base_of_v<PMDefaultStressTensor, T>;

template <typename T>
struct is_pm_dijkuizen_stress_tensor : std::bool_constant<PMDijkhuizenStressTreatment<T>> {};

template <typename T>
struct is_pm_default_stress_tensor : std::bool_constant<PMDefaultStressTreatment<T>> {};

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PM_INFORMATION_H_
