/*
 * MIT License
 *
 * Copyright (c) 2022 David Rieder

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

#ifndef UTILITIES_UTILS_HASHES_H_
#define UTILITIES_UTILS_HASHES_H_
#include <cstddef>
namespace dare::utils {

/*\struct _FNVparam
 * \brief Helpers for hashing with the Fowler-Noll-Vu (FNV) hash function
 *
 * Background of the FNV-hash function:
 *  hash <-- offset (depending on size of hash)
 *  for each byte
 *     hash = hash * FNV-prime
 *     hash = hash XOR byte_of_data
 *  end for
 *
 */
template <int bytes>
struct _FNVparam {
};

/*
 * specialization for 32 bit FNV parameters
 */
template <>
struct _FNVparam<4> {
    static const std::size_t prime{0x01000193};
    static const std::size_t offset{0x811c9dc5};
};

/*
 * specialization for 64 bit FNV parameters
 */
template <>
struct _FNVparam<8> {
    static const std::size_t prime{0x00000100000001B3};
    static const std::size_t offset{0xcbf29ce484222325};
};
}  // namespace dare::utils

#endif  // UTILITIES_UTILS_HASHES_H_
