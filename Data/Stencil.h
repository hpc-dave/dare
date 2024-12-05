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
template <typename Grid, typename SC, std::size_t N>
class CenterMatrixStencil {
    /*!
     * @brief dummy operator for compilation
     * @param v 
     */
    CenterMatrixStencil operator*(SC v) const {
        return *this;
    }
};

/*!
 * @brief dummy class for SFINAE
 * @tparam Grid type of grid
 */
template <typename Grid, typename SC, std::size_t N>
class CenterValueStencil {
};

/*!
 * @brief dummy class for SFINAE
 * @tparam Grid type of grid
 */
template <typename Grid, typename SC, std::size_t N>
class FaceMatrixStencil {
};

/*!
 * @brief dummy class for SFINAE
 * @tparam Grid type of grid
 */
template <typename Grid, typename SC, std::size_t N>
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
template <typename Grid, typename SC, std::size_t N>
CenterMatrixStencil<Grid, SC, N> operator*(SC v, const CenterMatrixStencil<Grid, SC, N>& s) {
    return s * v;
}

/*!
 * @brief multiplication operator with doubles
 * @tparam Grid
 * @param v
 * @param s
 * @return
 */
template <typename Grid, typename SC, std::size_t N>
CenterValueStencil<Grid, SC, N> operator*(SC v, const CenterValueStencil<Grid, SC, N>& s) {
    return s * v;
}

/*!
 * @brief multiplication operator with doubles
 * @tparam Grid
 * @param v
 * @param s
 * @return
 */
template <typename Grid, typename SC, std::size_t N>
FaceMatrixStencil<Grid, SC, N> operator*(SC v, const FaceMatrixStencil<Grid, SC, N>& s) {
    return s * v;
}

/*!
 * @brief multiplication operator with doubles
 * @tparam Grid
 * @param v
 * @param s
 * @return
 */
template <typename Grid, typename SC, std::size_t N>
FaceValueStencil<Grid, SC, N> operator*(SC v, const FaceValueStencil<Grid, SC, N>& s) {
    return s * v;
}

}  // namespace dare::Data


namespace dare {

template<typename T>
struct is_face_matrix_stencil_helper : std::false_type {
};
template <typename GridType, typename SC, std::size_t N>
struct is_face_matrix_stencil_helper<Data::FaceMatrixStencil<GridType, SC, N>> : std::true_type {
};
template <typename T>
struct is_face_matrix_stencil : is_face_matrix_stencil_helper<std::remove_cv_t<T>> {
};
template <typename T>
static const bool is_face_matrix_stencil_v = is_face_matrix_stencil<T>::value;

template <typename T>
struct is_center_matrix_stencil_helper : std::false_type {
};
template <typename GridType, typename SC, std::size_t N>
struct is_center_matrix_stencil_helper<Data::CenterMatrixStencil<GridType, SC, N>> : std::true_type {
};
template <typename T>
struct is_center_matrix_stencil : is_center_matrix_stencil_helper<std::remove_cv_t<T>> {
};
template <typename T>
static const bool is_center_matrix_stencil_v = is_center_matrix_stencil<T>::value;

template <typename T>
struct is_face_value_stencil_helper : std::false_type {
};
template <typename GridType, typename SC, std::size_t N>
struct is_face_value_stencil_helper<Data::FaceValueStencil<GridType, SC, N>> : std::true_type {
};
template <typename T>
struct is_face_value_stencil : is_face_value_stencil_helper<std::remove_cv_t<T>> {
};
template <typename T>
static const bool is_face_value_stencil_v = is_face_value_stencil<T>::value;

template <typename T>
struct is_center_value_stencil_helper : std::false_type {
};
template <typename GridType, typename SC, std::size_t N>
struct is_center_value_stencil_helper<Data::CenterValueStencil<GridType, SC, N>> : std::true_type {
};
template <typename T>
struct is_center_value_stencil : is_center_value_stencil_helper<std::remove_cv_t<T>> {
};
template <typename T>
static const bool is_center_value_stencil_v = is_center_value_stencil<T>::value;

template <typename T>
static const bool is_matrix_stencil_v = is_face_matrix_stencil_v<T> || is_center_matrix_stencil_v<T>;
template <typename T>
static const bool is_value_stencil_v = is_face_value_stencil_v<T> || is_center_value_stencil_v<T>;
}  // end namespace dare
#endif  // DATA_STENCIL_H_
