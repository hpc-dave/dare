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

#ifndef ALGORITHM_STATICINFORMATION_H_
#define ALGORITHM_STATICINFORMATION_H_

namespace dare {

template <typename T>
concept ContainsViscosity = requires { typename T::viscosity; };  // NOLINT

template <typename T>
concept ContainsDensity = requires { typename T::density; };  // NOLINT

template <typename T>
concept ContainsPorosity = requires { typename T::porosity; };  // NOLINT

template <typename T>
concept ContainsImplicitForce = requires { typename T::implicit_force; };  // NOLINT

template <typename T>
concept ContainsExplicitForce = requires { typename T::explicit_force; };  // NOLINT

template <typename T>
concept ContainsCompressible = requires { typename T::compressible; } || std::is_same_v<std::remove_cv_t<decltype(T::compressible)>, bool>;  // NOLINT

template <typename T>
concept ContainsTVD = requires { typename T::tvd; };  // NOLINT

template <typename T>
concept ContainsViscousStress = requires { typename T::viscous_stress; };  // NOLINT

template <typename T>
concept ContainsMomentumIterations = requires { typename T::momentum_iterations; };  // NOLINT

template <typename T>
concept ContainsContinuityIterations = requires { typename T::continuity_iterations; };  // NOLINT

template <typename T>
concept ContainsTimeSchemeConvective = requires { typename T::time_scheme_convective; };  // NOLINT
}  // namespace dare

#endif  // ALGORITHM_STATICINFORMATION_H_
