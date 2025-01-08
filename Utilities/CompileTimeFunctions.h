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

#ifndef UTILITIES_COMPILETIMEFUNCTIONS_H_
#define UTILITIES_COMPILETIMEFUNCTIONS_H_
#include <tuple>
#include <array>
#include <utility>

namespace dare::utils {

/*!
 * @brief converts a homogeneous tuple into an array at compile time
 * @tparam tuple_t type of the tuple
 * @param tuple actual tuple
 * @return std::array
 */
template <typename tuple_t>
constexpr auto convert_tuple_to_array(tuple_t&& tuple) {
    constexpr auto get_array = [](auto&&... x) { return std::array{std::forward<decltype(x)>(x)...}; };
    return std::apply(get_array, std::forward<tuple_t>(tuple));
}

}  // namespace dare::utils

#endif  // UTILITIES_COMPILETIMEFUNCTIONS_H_
