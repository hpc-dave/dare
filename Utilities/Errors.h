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

#ifndef UTILITIES_ERRORS_H_
#define UTILITIES_ERRORS_H_

#include <iostream>
#include <type_traits>

#define ERROR \
    std::cerr << "ERROR! In " << __FILE__ << ": " << __LINE__ << " -> "

#define ERROR_CLOSE ""<< std::endl

namespace dare {

/*!
 * @brief a workaround to allow for compile time false values
 * @tparam ...Ts Any type of template argument
 * 
 * \note The variable cannot be named 'AlwaysFalse', as it conflicts with
 * a function definition of google test
 */
template <typename... Ts>
static const bool always_false = std::false_type::value;

}  // end namespace dare

#endif  // UTILITIES_ERRORS_H_
