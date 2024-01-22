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

#ifndef UTILITIES_VECTOR_TRAITS_H_
#define UTILITIES_VECTOR_TRAITS_H_

#include <cstddef>
#include <type_traits>

namespace dare::utils {

/*! \class VectorBase
 * @tparam N number of elements in the data set
 * @tparam T type of variables
 * The main target of this class is to provide basic functionality and memory allocation
 * for Vector class. It serves as parent class for the decorator below.
 * \note all common functionality is currently implemented by the Vector class
 */
template <std::size_t N, typename T>
class VectorBase {
public:
    /*!
     * @brief default constructor
     */
    VectorBase() {}
protected:
    T _data[N];
};

/*! \class VectorDecorator
 * @tparam N specialization for n-th element
 * @tparam Dim number of elements in the data set
 * @tparam T type of variables
 * @tparam Enable SFINAE type placeholder
 * This class uses the decorator design pattern to provide specialized functionality to the
 * Vector class, depending on the provided type. Good examples are the
 * getters and setters, which might be different.
 * The default version doesn't add anything.
 */
template <std::size_t N, std::size_t Dim, typename T, typename Enable = void>
class VectorDecorator : public VectorBase<Dim, T> {
public:
    /*!
     * @brief default constructor
     */
    VectorDecorator() : VectorBase<Dim, T>(){}
};

/*! \class VectorDecorator
 * \brief specialization for integer types
 * @tparam N specialization for n-th element
 * @tparam Dim number of elements in the data set
 * @tparam T type of integer
 */
template <std::size_t N, std::size_t Dim, typename T>
class VectorDecorator<N, Dim, T, std::enable_if_t<std::is_integral_v<T>>> : public VectorDecorator<N - 1, Dim, T> {
public:
    VectorDecorator() : VectorDecorator<N - 1, Dim, T>() {}
};

/*! \class VectorDecorator
 * \brief specialization for integer types and N == 0
 * @tparam Dim number of elements in the data set
 * @tparam T type of integer
 */
template <std::size_t Dim, typename T>
class VectorDecorator<0, Dim, T, std::enable_if_t<std::is_integral_v<T>>> : public VectorBase<Dim, T> {
public:
    VectorDecorator() : VectorBase<Dim, T>() {}
};

/*! \class VectorDecorator
 * \brief specialization for integer types and N == 1
 * @tparam Dim number of elements in the data set
 * @tparam T type of integer
 */
template <std::size_t Dim, typename T>
class VectorDecorator<1, Dim, T, std::enable_if_t<std::is_integral_v<T>>> : public VectorDecorator<0, Dim, T> {
public:
    VectorDecorator() : VectorDecorator<0, Dim, T>() {}

    /*!
     * \brief getter for first variable
     */
    T& i() { return this->_data[0]; }

    /*!
     * \brief const getter for first variable
     */
    T i() const { return this->_data[0]; }
};

/*! \class VectorDecorator
 * \brief specialization for integer types and N == 2
 * @tparam Dim number of elements in the data set
 * @tparam T type of integer
 */
template <std::size_t Dim, typename T>
class VectorDecorator<2, Dim, T, std::enable_if_t<std::is_integral_v<T>>> : public VectorDecorator<1, Dim, T> {
public:
    VectorDecorator() : VectorDecorator<1, Dim, T>() {}

    /*!
     * \brief getter for second variable
     */
    T& j() { return this->_data[1]; }

    /*!
     * \brief const getter for second variable
     */
    T j() const { return this->_data[1]; }
};

/*! \class VectorDecorator
 * \brief specialization for integer types and N == 2
 * @tparam Dim number of elements in the data set
 * @tparam T type of integer
 */
template <std::size_t Dim, typename T>
class VectorDecorator<3, Dim, T, std::enable_if_t<std::is_integral_v<T>>> : public VectorDecorator<2, Dim, T> {
public:
    VectorDecorator() : VectorDecorator<2, Dim, T>() {}

    /*!
     * \brief getter for second variable
     */
    T& k() { return this->_data[2]; }

    /*!
     * \brief const getter for second variable
     */
    T k() const { return this->_data[2]; }
};

/*! \class VectorDecorator
 * \brief specialization for floating point types
 * @tparam N specialization for n-th element
 * @tparam Dim number of elements in the data set
 * @tparam T type of floating point numbers
 */
template <std::size_t N, std::size_t Dim, typename T>
class VectorDecorator<N, Dim, T, std::enable_if_t<std::is_floating_point_v<T>>>
        : public VectorDecorator<N - 1, Dim, T> {
public:
    /*!
     * \brief constructor
     * @param  pointer to beginning of data
     */
    VectorDecorator() : VectorDecorator<N - 1, Dim, T>() {}
};

/*! \class VectorDecorator
 * \brief specialization for floating point types and N == 1
 * @tparam Dim number of elements in the data set
 * @tparam T type of floating point numbers
 */
template <std::size_t Dim, typename T>
class VectorDecorator<1, Dim, T, std::enable_if_t<std::is_floating_point_v<T>>>
    : public VectorBase<Dim, T>{
public:
    VectorDecorator() {}

    /*!
     * \brief getter for first variable
     */
    T& x() { return this->_data[0]; }

    /*!
     * \brief const getter for first variable
     */
    T x() const { return this->_data[0]; }
};

/*! \class VectorDecorator
 * \brief specialization for floating point types and N == 2
 * @tparam Dim number of elements in the data set
 * @tparam T type of floating point numbers
 */
template <std::size_t Dim, typename T>
class VectorDecorator<2, Dim, T, std::enable_if_t<std::is_floating_point_v<T>>> : public VectorDecorator<1, Dim, T> {
public:
    VectorDecorator() : VectorDecorator<1, Dim, T>() {}

    /*!
     * \brief getter for second variable
     */
    T& y() { return this->_data[1]; }

    /*!
     * \brief const getter for second variable
     */
    T y() const { return this->_data[1]; }
};

/*! \class VectorDecorator
 * \brief specialization for floating point types and N == 3
 * @tparam Dim number of elements in the data set
 * @tparam T type of floating point numbers
 */
template <std::size_t Dim, typename T>
class VectorDecorator<3, Dim, T, std::enable_if_t<std::is_floating_point_v<T>>> : public VectorDecorator<2, Dim, T> {
public:
    VectorDecorator() : VectorDecorator<2, Dim, T>() {}

    /*!
     * \brief getter for third variable
     */
    T& z() { return this->_data[2]; }

    /*!
     * \brief const getter for third variable
     */
    T z() const { return this->_data[2]; }
};

}  // namespace dare::utils

#endif  // UTILITIES_VECTOR_TRAITS_H_
