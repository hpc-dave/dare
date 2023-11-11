/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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

#ifndef DATA_STENCIL_H_
#define DATA_STENCIL_H_

#include <string>
#include <Kokkos_Core.hpp>

namespace dare::Data {

template <typename T>
class Stencil {
public:
    using Ordinal = int;

    /*!
     * @brief constructor
     * @param identifier name provided for debugging
     * @param size size of the stencil
     */
    explicit Stencil(std::string identifier = "unspecified", Ordinal size = 0);

    /*!
     * @brief copy constructor
     * @param other instance to copy from
     * \note This is required, since the underlying data structure otherwise will
     * do a shallow copy of the data
     */
    Stencil(const Stencil<T>& other);

    /*!
     * @brief copy assignment
     * @param other instance to copy from
     */
    Stencil<T>& operator=(const Stencil<T>& other);

    /*!
     * @brief resized the data structure
     * @param stencil_size size of the stencil
     */
    void Resize(Ordinal stencil_size);

    /*!
     * @brief returns size
     */
    Ordinal GetSize() const;

    /*!
     * @brief sets value at certain position
     * @param pos position
     * @param value value to set
     * No range check will be conducted!
     */
    void InsertValue(Ordinal pos, const T& value);

    /*!
     * @brief getter
     * @param pos position in array
     */
    T& GetValue(Ordinal pos);

    /*!
     * @brief constant getter
     * @param pos position in array
     */
    const T& GetValue(Ordinal pos) const;

    /*!
     * @brief returns raw data set
     */
    const Kokkos::View<T*>& GetData() const;

private:
    Kokkos::View<T*> values;  //!< data array
};

}  // namespace dare::Data

#include "Stencil.inl"

#endif  // DATA_STENCIL_H_
