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

#ifndef UTILITIES_VECTOR_TRAITS_H_
#define UTILITIES_VECTOR_TRAITS_H_

#include <cstddef>
#include <type_traits>

namespace dare::utils {

/*! \struct _vector_functor
 * @tparam N number of elements in the data set
 * @tparam T type of variables
 * @tparam Enable SFINAE type placeholder
 * The main target of this class is to provide specialized functionality to the
 * Vector class, depending on the provided type. Good examples are the
 * getters and setters, which might be different.
 */
template <std::size_t N, typename T, typename Enable = void>
struct _vector_functor {
    explicit _vector_functor(T* ptr) {}
};

/*! \struct _vector_functor
 * \brief specialization for integer types
 * @tparam N number of elements in the data set
 * @tparam T type of integer
 */
template <std::size_t N, typename T>
struct _vector_functor<N, T, std::enable_if_t<std::is_integral_v<T>>> : public _vector_functor<N-1, T> {
    explicit _vector_functor(T* ptr) : _vector_functor<N - 1, T>(ptr) {}
};

/*! \struct _vector_functor
 * \brief specialization for integer types and N == 1
 * @tparam T type of integer
 */
template <typename T>
struct _vector_functor<1, T, std::enable_if_t<std::is_integral_v<T>>> {
    explicit _vector_functor(T* ptr) : p_data(ptr) {}

    /*!
     * \brief getter for first variable
     */
    T& i() { return p_data[0]; }

    /*!
     * \brief const getter for first variable
     */
    T i() const { return p_data[0]; }

    T* p_data;  //!< pointer to data
};

/*! \struct _vector_functor
 * \brief specialization for integer types and N == 2
 * @tparam T type of integer
 */
template <typename T>
struct _vector_functor<2, T, std::enable_if_t<std::is_integral_v<T>>> : public _vector_functor<1, T> {
    explicit _vector_functor(T* ptr) : _vector_functor<1, T>(ptr) {}

    /*!
     * \brief getter for second variable
     */
    T& j() { return this->p_data[1]; }

    /*!
     * \brief const getter for second variable
     */
    T j() const { return this->p_data[1]; }
};

/*! \struct _vector_functor
 * \brief specialization for integer types and N == 2
 * @tparam T type of integer
 */
template <typename T>
struct _vector_functor<3, T, std::enable_if_t<std::is_integral_v<T>>> : public _vector_functor<2, T> {
    explicit _vector_functor(T* ptr) : _vector_functor<2, T>(ptr) {}

    /*!
     * \brief getter for second variable
     */
    T& k() { return this->p_data[2]; }

    /*!
     * \brief const getter for second variable
     */
    T k() const { return this->p_data[2]; }
};

/*! \struct _vector_functor
 * \brief specialization for floating point types
 * @tparam N number of elements in the data set
 * @tparam T type of floating point numbers
 */
template <std::size_t N, typename T>
struct _vector_functor<N, T, std::enable_if_t<std::is_floating_point_v<T>>> : public _vector_functor<N - 1, T> {
    /*!
     * \brief constructor
     * @param ptr pointer to beginning of data
     */
    explicit _vector_functor(T* ptr) : _vector_functor<N - 1, T>(ptr) {}
};

/*! \struct _vector_functor
 * \brief specialization for floating point types and N == 1
 * @tparam T type of floating point numbers
 */
template <typename T>
struct _vector_functor<1, T, std::enable_if_t<std::is_floating_point_v<T>>> {
    explicit _vector_functor(T* ptr) : p_data(ptr) {}

    /*!
     * \brief getter for first variable
     */
    T& x() { return p_data[0]; }

    /*!
     * \brief const getter for first variable
     */
    T x() const { return p_data[0]; }

    T* p_data;  //! pointer to data
};

/*! \struct _vector_functor
 * \brief specialization for floating point types and N == 2
 * @tparam T type of floating point numbers
 */
template <typename T>
struct _vector_functor<2, T, std::enable_if_t<std::is_floating_point_v<T>>> : public _vector_functor<1, T> {
    explicit _vector_functor(T* ptr) : _vector_functor<1, T>(ptr) {}

    /*!
     * \brief getter for second variable
     */
    T& y() { return this->p_data[1]; }

    /*!
     * \brief const getter for second variable
     */
    T y() const { return this->p_data[1]; }
};

/*! \struct _vector_functor
 * \brief specialization for floating point types and N == 3
 * @tparam T type of floating point numbers
 */
template <typename T>
struct _vector_functor<3, T, std::enable_if_t<std::is_floating_point_v<T>>> : public _vector_functor<2, T> {
    explicit _vector_functor(T* ptr) : _vector_functor<2, T>(ptr) {}

    /*!
     * \brief getter for third variable
     */
    T& z() { return this->p_data[2]; }

    /*!
     * \brief const getter for third variable
     */
    T z() const { return this->p_data[2]; }
};

}  // namespace dare::utils

#endif  // UTILITIES_VECTOR_TRAITS_H_
