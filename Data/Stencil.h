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

#ifndef DATA_STENCIL_H_
#define DATA_STENCIL_H_

#include <string>

namespace dare::Data {

/*!
 * @brief dummy class for SFINAE
 * @tparam Grid type of grid
 */
template <typename Grid>
class CenterMatrixStencil {
    /*!
     * @brief dummy operator for compilation
     * @param v 
     */
    CenterMatrixStencil operator*(typename Grid::ScalarType v) const {
        return *this;
    }
};

/*!
 * @brief dummy class for SFINAE
 * @tparam Grid type of grid
 */
template <typename Grid>
class CenterValueStencil {
};

/*!
 * @brief dummy class for SFINAE
 * @tparam Grid type of grid
 */
template <typename Grid>
class FaceMatrixStencil {
};

/*!
 * @brief dummy class for SFINAE
 * @tparam Grid type of grid
 */
template <typename Grid>
class FaceValueStencil {
};

// operators for improved use

/*!
 * @brief multiplication operator with doubles
 * @tparam Grid 
 * @param v 
 * @param s 
 * @return 
 */
template <typename Grid>
CenterMatrixStencil<Grid> operator*(typename Grid::ScalarType v, const CenterMatrixStencil<Grid>& s) {
    return s * v;
}

/*!
 * @brief multiplication operator with doubles
 * @tparam Grid
 * @param v
 * @param s
 * @return
 */
template <typename Grid>
CenterValueStencil<Grid> operator*(typename Grid::ScalarType v, const CenterValueStencil<Grid>& s) {
    return s * v;
}

}  // namespace dare::Data

#endif  // DATA_STENCIL_H_
