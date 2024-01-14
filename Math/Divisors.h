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

#ifndef MATH_DIVISORS_H_
#define MATH_DIVISORS_H_

namespace dare::math {

/*!
 * \brief computes divisor at compile time
 * @tparam DIV integral value by which is divided
 * @tparam T type of variable
 * @tparam TEnable make sure the divisor is a multiple of 2 and not integral
 */
template <typename T,
          std::size_t DIV,
          typename TEnable = std::enable_if_t<!std::is_integral_v<T> && (DIV % 2 == 0 || DIV == 1)>>
constexpr T Divisor() {
    if constexpr (DIV == 1) {
        return static_cast<T>(1.);
    } else {
        return Divisor<T, DIV / 2>() * static_cast<T>(0.5);
    }
}

/*!
 * \brief Divides the provided floating point value by an integral value which is pow(2)
 * @tparam DIV integral value by which is divided
 * @tparam T type of variable
 * @tparam TEnable make sure the divisor is a multiple of 2
 */
template <std::size_t DIV,
          typename T,
          typename TEnable = std::enable_if_t<(DIV % 2 == 0) || DIV == 1>>
[[nodiscard]] T Divide(T value) {
    if constexpr (std::is_integral_v<T>) {
        return value / DIV;
    } else {
        return value * Divisor<T, DIV>();
    }
}

}  // end namespace dare::math


#endif  // MATH_DIVISORS_H_
