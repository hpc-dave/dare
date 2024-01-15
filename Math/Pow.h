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

#ifndef MATH_POW_H_
#define MATH_POW_H_

namespace dare::math {

// template <std::size_t BASE, std::size_t EXP, typename TEnable = std::enable_if_t<std::is_integral_v<BASE>>>
// [[nodiscard]] constexpr std::size_t Pow() {
//     if constexpr (EXP == 0) {
//         return 1;
//     } else {
//         return BASE * Pow<BASE, EXP - 1>();
//     }
// }
template <int BASE, std::size_t EXP>
[[nodiscard]] constexpr std::size_t Pow() {
    if constexpr (EXP == 0) {
        return 1;
    } else {
        return BASE * Pow<BASE, EXP - 1>();
    }
}

template <int EXP, typename BASE>
[[nodiscard]] BASE Pow(BASE base) {
    if constexpr (EXP == 0) {
        return 1;
    } else if constexpr(EXP < 0) {
        return Pow<-EXP>(1. / base);
    } else {
        return base * Pow<EXP - 1>(base);
    }
}

template <typename TBase, typename TExp, typename TEnable = std::enable_if_t<std::is_integral_v<TExp>>>
[[nodiscard]] TBase Pow(TBase base, TExp exp) {
    TBase res{1};
    if constexpr (std::is_signed_v<TExp>) {
        if (exp < 0) {
            base = static_cast<TBase>(1) / base;
            exp = -exp;
        }
    }
    for (TExp n{0}; n < exp; ++n) {
        res *= base;
    }
    return res;
}
}  // end namespace dare::math

#endif  // MATH_POW_H_
