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

#ifndef MPI_OBLACKHOLESTREAM_H_
#define MPI_OBLACKHOLESTREAM_H_

#include <iostream>
#include <string>

namespace dare::mpi {

/*!
 * \brief A convenient way to swallow output
 * In certain situations (e.g. when working in parallel environments)
 * the output needs to be managed. This class can be used instead of
 * std::cout if output should be suppressed.
 */
class BlackHoleOStream
    : virtual public std::basic_ostream<char, std::char_traits<char>> {
public:
    /*!
     * \brief constructor
     * the underlying ostream will never print anything
     */
    BlackHoleOStream() : std::basic_ostream<char, std::char_traits<char>>(NULL) {}
};

}  // namespace dare::mpi

#endif  // MPI_OBLACKHOLESTREAM_H_
