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

#ifndef UTILITIES_VECTOR_H_
#define UTILITIES_VECTOR_H_

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>
#include <array>
#include <string>

#include "Hashes.h"
#include "Vector_traits.h"
#include "Errors.h"

namespace dare {

template <typename T, typename... Ts>
using AllConvertible = std::enable_if_t<std::conjunction_v<std::is_convertible<T, Ts>...>>;

/*! \class Vector
 * \tparam N number of values in the data set
 * \tparam T type of data
 * \brief general class for vector like data
 * In principle, this class stores a 1D data set of arbitrary type,
 * For integer and floating point values it provides a set of optimization
 * and convenience functions, i.e. determining the length of a vector or
 * computing a dot/cross product. Used data types have to support
 * ==, !=, +, -, *, /, +=, -=, *=, /=-operators, as well as being copyable
 */
template <std::size_t N, typename T = double>
class Vector : public VectorDecorator<N, N, T>{
public:
    using InternalType = T;
    using ContainerType = typename VectorBase<N, T>::ContainerType;
    using difference_type = typename ContainerType::difference_type;
    using iterator = typename ContainerType::iterator;
    using const_iterator = typename ContainerType::const_iterator;
    using reverse_iterator = typename ContainerType::reverse_iterator;
    using const_reverse_iterator = typename ContainerType::const_reverse_iterator;

    /*!
     * \brief constructor
     * @param args variable number of initialization arguments
     * Takes a variable amount of input values for construction. Those values
     * have to be the same type as the specified template type T and a maximum
     * of N values can be provided, missing values will be given the value 0
     * \warning implicit conversion is NOT allowed, depending on the compiler you
     * might receive strange error-message
     */
    template <typename... Ts,
              typename = AllConvertible<T, Ts...>,
              typename = std::enable_if_t<sizeof...(Ts) <= N>>
    explicit Vector(const Ts&... args);

    /*!
     * copy constructor
     * @param other object to copy from
     */
    template <typename A,
              typename = std::enable_if_t<std::is_convertible_v<A, T>> >
    Vector(const Vector<N, A>& other);

    /*!
     * returns direct access to raw data
     * @return Value type pointer
     */
    T* data();

    /*!
     * provides direct access to the data
     */
    const T* data() const;

    /*!
     * assignment constructor
     *
     * @param other object to copy from
     */
    template <typename A,
              typename = std::enable_if_t<std::is_convertible_v<A, T>>>
    Vector<N, T>& operator=(const Vector<N, A>& other);

    /*!
     * access operator
     *
     * @param n position to access
     * \note without NDEBUG, a bounds check will be conducted
     */
    T& operator[](std::size_t n);

    /*!
     *  access operator
     *
     * @param n position to access
     * \note without NDEBUG, a bounds check will be conducted
     */
    const T& operator[](std::size_t n) const;

    /*!
     * addition of other vector
     *
     * @param other addition partner
     */
    Vector<N, T> operator+(const Vector<N, T>& other) const;

    /*!
     * addition of a single value to all internal values
     *
     * @param other addition partner
     */
    Vector<N, T> operator+(const T& val) const;

    /*!
     * addition of another vector to the current instance
     *
     * @param other addition partner
     */
    void operator+=(const Vector<N, T>& other);

    /*!
     * addition of a single value to all internal values
     *
     * @param other addition partner
     */
    void operator+=(const T& val);

    /*!
     * subtraction of other vector
     * @param other addition partner
     */
    Vector<N, T> operator-(const Vector<N, T>& other) const;

    /*!
     * subtraction of single value from all components
     * @param val value to subtract
     */
    Vector<N, T> operator-(const T& val) const;

    /*!
     * provides negative of the values
     */
    Vector<N, T> operator-() const;

    /*!
     * -= operator
     * @param other other vector
     */
    void operator-=(const Vector<N, T>& other);

    /*!
     * subtraction of single value from all components
     * @param val value to subtract
     */
    void operator-=(const T& val);

    /*!
     * multiplication
     * @param other vector to mulitply with
     */
    Vector<N, T> operator*(const Vector<N, T>& other) const;

    /*!
     * elementwise multiplication
     * @param val value to multiply with
     */
    Vector<N, T> operator*(const T& val) const;

    /*!
     * multiplication
     * @param other vector for multiplication
     */
    void operator*=(const Vector<N, T>& other);

    /*!
     * elementwise multiplication
     * @param val value to multiply with
     */
    void operator*=(const T& val);

    /*!
     * division
     * @param other vector to divide with
     */
    Vector<N, T> operator/(const Vector<N, T>& other) const;

    /*!
     * elementwise division
     * @param val value to divide with
     */
    Vector<N, T> operator/(const T& val) const;

    /*!
     * division
     * @param other vector to divide with
     */
    void operator/=(const Vector<N, T>& other);

    /*!
     * division
     * @param val value to divide with
     */
    void operator/=(const T& val);

    /*!
     * comparison operator
     * @param other vector to compare with
     */
    bool operator==(const Vector<N, T>& other) const;

    /*!
     * non-equal operator
     * @param other vector to compare with
     */
    bool operator!=(const Vector<N, T>& other) const;

    /*!
     * returns number of elements
     */
    constexpr std::size_t size() const;

    /*!
     * returns length of the vector
     * \note only sensible for floating point data, but also enabled for other types
     */
    T length() const;

    /*!
     * return forward iterator to first object
     */
    iterator begin();

    /*!
     * returns constant forward iterator
     */
    const_iterator begin() const;

    /*!
     * returns constant forward iterator
     */
    const_iterator cbegin() const;

    /*!
     * return reverse iterator to first object
     */
    reverse_iterator rbegin();

    /*!
     * returns constant reverse iterator
     */
    const_reverse_iterator rbegin() const;

    /*!
     * returns constant reverse iterator
     */
    const_reverse_iterator crbegin() const;

    /*!
     * returns iterator to end of data
     */
    iterator end();

    /*!
     * constant iterator to end of data
     */
    const_iterator end() const;

    /*!
     * constant iterator to end of data
     */
    const_iterator cend() const;

    /*!
     * returns iterator to end of data
     */
    reverse_iterator rend();

    /*!
     * constant iterator to end of data
     */
    const_reverse_iterator rend() const;

    /*!
     * constant iterator to end of data
     */
    const_reverse_iterator crend() const;

    /*!
     * \tparam I value to access in the data set
     * \tparam Ts parameter .pack with remainder
     * @param arg argument to set at I
     * @param args parameter pack with remaining values
     */
    template <std::size_t I = 0,
              typename A,
              typename... Ts,
              typename = std::is_convertible<A, T>,
              typename = AllConvertible<T, Ts...>,
              typename = std::enable_if_t<sizeof...(Ts) <= N>>
    void SetValues(const A& arg, const Ts&... args);

    /*!
     * sets all values to default values
     * \tparam I value to start with setting the default values
     */
    template <std::size_t I>
    void SetValues();

    /*!
     * @brief sets all values to specified value
     * @param val value to apply
     */
    template<typename A, typename = std::is_convertible<A, T>>
    void SetAllValues(const A& val);

    /*!
     * computes dot product with another vector
     * @param other vector to compute the dot product with
     */
    T dot(const Vector<N, T>& other) const;

    /*!
     * computes cross product
     * @param other vector to compute the cross product with
     * \note only enabled for 3D vectors, for higher dimensions a more general algorithm is required
     */
    template <std::size_t Ns = N>
    typename std::enable_if<(Ns == 3), Vector<N, T>>::type cross(const Vector<N, T>& other) const;

    /*!
     * outputs the vector
     */
    template <typename OS>
    friend OS& operator<<(OS& os, const Vector<N, T>& v) {
        os  << v[0];
        for (std::size_t i{1}; i < N; ++i)
            os << ' ' << v[i];

        return os;
    }

    /*!
     * @brief provides pretty printing as a string value
     * 
     * \note mostly required for clang
     */
    std::string ToString() const {
        std::ostringstream os;
        os << *this;
        return os.str();
    }

    /*
     * calculates hash for hash-maps
     */
    std::size_t GetHash() const;

    /*!
     * @brief returns sum of all elements
     */
    T AllSum() const;

    /*!
     * @brief returns sum of all absolute values
     */
    T AllAbsSum() const;

private:
    /*!
     * iterates over the internal values and executes arbitrary manipulation
     * @param lambda operation to execute per data entry
     * @param op join operation for reduction
     */
    template <std::size_t I = 0, typename Expr, typename Op>
    auto IterateValues(Expr lambda, Op op);

    /*!
     * iterates over the internal values and executes arbitrary manipulation
     * @param lambda operation to execute per data entry
     * @param op join operation for reduction
     */
    template <std::size_t I = 0, typename Expr, typename Op>
    auto IterateValues(Expr lambda, Op op) const;
};

/*!
 * @brief multiplication operator for convenience
 * @tparam T basic type
 * @tparam N number of elements
 * @param v1 value to multiply with
 * @param v2 vector
 * @return vector
 */
template <std::size_t N, typename T>
dare::Vector<N, T> operator*(const T& v1, const dare::Vector<N, T>& v2) {
    return v2 * v1;
}

/*!
 * @brief division operator for convenience
 * @tparam T basic type
 * @tparam N number of elements
 * @param v1 value to multiply with
 * @param v2 vector
 * @return vector
 */
template <std::size_t N, typename T>
dare::Vector<N, T> operator/(const T& v1, const dare::Vector<N, T>& v2) {
    dare::Vector<N, T> v;
    v.SetAllValues(v1);
    for (std::size_t n{0}; n < N; n++) {
        v[n] /= v2[n];
    }
    return v;
}

}  // namespace dare

namespace std {
/*
 * \brief specialization of the hash-function for Vector
 */
template <std::size_t N, typename T>
class hash<dare::Vector<N, T>> {
public:
    /*
     * \brief returns the hash for a Vector
     */
    [[nodiscard]] std::size_t operator()(const dare::Vector<N, T>& v) const {
        return v.GetHash();
    }
};

/*!
 * \brief specialization for dare::Vector: computes the absolute value
 * @tparam N number of elements in vector
 * @tparam T type of variable
 * @param v input vector
 * @return vector with all values absolute
 */
template<std::size_t N, typename T>
[[nodiscard]] dare::Vector<N, T> abs(dare::Vector<N, T> v) {
    for (auto& e : v)
        e = std::abs(e);
    return v;
}

}  // namespace std

#include "Vector.inl"
#endif  // UTILITIES_VECTOR_H_
