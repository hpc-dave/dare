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

#include "Utilities/PropertyInformation.h"

namespace dare::algorithm {

enum class PMProperties {
    Density,
    Viscosity,
    Porosity,
    ImplicitForce,
    ExplicitForce,
    Compressible
};

enum class PMSolver {
    Momentum
};


template <bool Flag, PMProperties Tag>
using PMFlaggedInfo = utils::FlaggedInfo<Flag, PMProperties, Tag>;

template <typename Type, PMProperties Tag>
using PMTypeInfo = utils::TaggedTypeInfo<Type, PMProperties, Tag>;

template <typename Type, PMProperties Tag, std::size_t N>
using PMTypeCountedInfo = utils::TaggedCountedTypeInfo<Type, PMProperties, Tag, N>;

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

namespace detail {

template <typename T>
concept ContainsViscosity = requires { typename T::viscosity; };    // NOLINT

template <typename T>
concept ContainsDensity = requires { typename T::density; };  // NOLINT

template <typename T>
concept ContainsPorosity = requires { typename T::porosity; };  // NOLINT

template <typename T>
concept ContainsImplicitForce = requires { typename T::implicit_force; };  // NOLINT

template <typename T>
concept ContainsExplicitForce = requires { typename T::explicit_force; };  // NOLINT

template <typename T>
concept ContainsCompressible = requires { typename T::compressible; };  // NOLINT

template <typename TDict, ContainsViscosity TDictDefault>
constexpr auto PMGetViscosity() {
    return typename TDictDefault::viscosity();
}

template <ContainsViscosity TDict, ContainsViscosity TDictDefault>
    requires utils::TaggedTypeInfoType<typename TDict::viscosity>
constexpr auto PMGetViscosity() {
    return typename TDict::viscosity();
}

template <typename TDict, ContainsDensity TDictDefault>
constexpr auto PMGetDensity() {
    return typename TDictDefault::density();
}

template <ContainsViscosity TDict, ContainsDensity TDictDefault>
    requires utils::TaggedTypeInfoType<typename TDict::density>
constexpr auto PMGetDensity() {
    return typename TDict::density();
}

template <typename TDict, ContainsPorosity TDictDefault>
constexpr auto PMGetPorosity() {
    return typename TDictDefault::porosity();
}

template <ContainsPorosity TDict, ContainsPorosity TDictDefault>
    requires utils::TaggedTypeInfoType<typename TDict::porosity>
constexpr auto PMGetPorosity() {
    return typename TDict::porosity();
}

template <typename TDict, ContainsImplicitForce TDictDefault>
constexpr auto PMGetImplicitForce() {
    return typename TDictDefault::implicit_force();
}

template <ContainsImplicitForce TDict, ContainsImplicitForce TDictDefault>
    requires utils::TaggedCountedTypeInfoType<typename TDict::implicit_force>
constexpr auto PMGetImplicitForce() {
    return typename TDict::implicit_force();
}

template <typename TDict, ContainsExplicitForce TDictDefault>
constexpr auto PMGetExplicitForce() {
    return typename TDictDefault::explicit_force();
}

template <ContainsExplicitForce TDict, ContainsExplicitForce TDictDefault>
    requires utils::TaggedCountedTypeInfoType<typename TDict::explicit_force>
constexpr auto PMGetExplicitForce() {
    return typename TDict::explicit_force();
}

template <typename TDict, ContainsCompressible TDictDefault>
constexpr auto PMGetCompressibility() {
    return typename TDictDefault::compressible();
}

template <ContainsCompressible TDict, ContainsCompressible TDictDefault>
    requires std::is_same_v<decltype(TDict::compressible), bool>
constexpr auto PMGetCompressibility() {
    return typename TDict::compressible();
}

}  // namespace detail

}  // namespace dare::algorithm

#endif  // ALGORITHM_NAVIERSTOKES_PM_INFORMATION_H_
