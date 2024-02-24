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
 * \brief tests if an integral number is a integral root of another
 * @tparam NUM number which is tested
 * @tparam ROOT potential root of the number
 * Note, that root here refers to the possibility to retain the original number
 * with an exponent
 */
template<int NUM, int ROOT>
constexpr bool IsRootOf() {
    if constexpr (NUM == 0)
        return false;
    else if constexpr (ROOT == 1)
        return true;
    else if constexpr (NUM % ROOT > 0)
        return false;
    else if constexpr (NUM == ROOT)
        return true;
    else
        return IsRootOf<NUM / ROOT, ROOT>();
}

/*!
 * \brief computes divisor at compile time
 * @tparam DIV integral value by which is divided
 * @tparam T type of variable
 * @tparam TEnable check to make sure we don't try to get the divisor for an integral value
 * Note, that division of integral values is heavily optimized by bitshifting. A divisor will
 * not lead to the desired result and may lead to a loss in accuracy!
 */
template<typename T, int Denominator, typename TEnable = std::enable_if_t<!std::is_integral_v<T>>>
constexpr T Divisor() {
    static_assert(Denominator != 0);
    return static_cast<T>(1) / static_cast<T>(Denominator);
}

/*!
 * \brief Divides the provided value by an integral value
 * @tparam DIV integral value by which is divided
 * @tparam T type of variable
 * @tparam TEnable make sure the divisor has the root 2
 */
template <int DIV,
          typename T>
[[nodiscard]] T Divide(T value) {
    if constexpr (std::is_integral_v<T>) {
        return value / DIV;
    } else {
        return value * Divisor<T, DIV>();
    }
}

/*!
 * @brief provides division value at compile time
 * @tparam T type of value we want to get
 * @tparam Nominator integer nominator
 * @tparam Denominator integer denominator
 */
template <typename T, int Nominator, int Denominator>
[[nodiscard]] constexpr T Divide() {
    if constexpr (std::is_integral_v<T>) {
        static_assert(IsRootOf<Nominator, Denominator>,
                      "I won't allow the loss of accuracy by dividing integer by not a root!");
    }
    return static_cast<T>(Nominator) / static_cast<T>(Denominator);
}

}  // end namespace dare::math


#endif  // MATH_DIVISORS_H_
