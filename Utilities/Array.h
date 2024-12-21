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

#ifndef UTILITIES_ARRAY_H_
#define UTILITIES_ARRAY_H_

#include "Utilities/Vector.h"

namespace dare::utils {

template<std::size_t M, std::size_t N, typename T>
class Array {
public:
    using InnerVec = typename dare::utils::Vector<N, T>;            //!< holds data for each m
    using OuterVec = typename dare::utils::Vector<M, InnerVec>;     //!< holds m vectors
    using SelfType = Array<M, N, T>;                                //!< its own type

    /*!
     * @brief access to inner vectors
     * @param m vector at position m
     */
    InnerVec& operator[](std::size_t m);

    /*!
     * @brief const access to inner vectors
     * @param m vector at position m
     */
    const InnerVec& operator[](std::size_t m) const;

    /*!
     * @brief access to inner vectors
     * @param m vector at position m
     */
    InnerVec& Get(std::size_t m);

    /*!
     * @brief const access to inner vectors
     * @param m vector at position m
     */
    const InnerVec& Get(std::size_t m) const;

    /*!
     * @brief access to data
     * @param m position in outer vector
     * @param n position in inner vector
     */
    T& Get(std::size_t m, std::size_t n);

    /*!
     * @brief const access to data
     * @param m position in outer vector
     * @param n position in inner vector
     */
    const T& Get(std::size_t m, std::size_t n) const;

private:
    OuterVec data;
};

}  // namespace dare::utils

#include "Array.inl"
#endif  // UTILITIES_ARRAY_H_
