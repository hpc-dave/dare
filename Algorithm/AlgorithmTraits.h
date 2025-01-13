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

#ifndef ALGORITHM_ALGORITHMTRAITS_H_
#define ALGORITHM_ALGORITHMTRAITS_H_

namespace dare::algorithm {

/*!
 * \brief an tagging struct for Newton iterations
*/
struct Newton {
};

template <typename T>
concept NewtonIterations =
    std::is_base_of<Newton, std::remove_cv_t<T>>::value;

template <typename T>
struct uses_newton_iterations : std::bool_constant<NewtonIterations<T>>{};

template <typename T>
constexpr bool uses_newton_iterations_v = uses_newton_iterations<T>::value;

/*!
 * \brief an tagging class for FixedPoint iterations
 */
struct FixedPoint {
};

template <typename T>
concept FixedPointIterations =
    std::is_base_of<FixedPoint, std::remove_cv_t<T>>::value;

template <typename T>
struct uses_fixed_point_iterations : std::bool_constant<FixedPointIterations<T>> {};

template <typename T>
constexpr bool uses_fixed_point_iterations_v = uses_fixed_point_iterations<T>::value;

}  // namespace dare::algorithm

#endif  // ALGORITHM_ALGORITHMTRAITS_H_
