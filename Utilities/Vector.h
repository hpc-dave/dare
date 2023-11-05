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

#ifndef UTILITIES_VECTOR_H_
#define UTILITIES_VECTOR_H_

#include <iostream>
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <type_traits>

#include "Hashes.h"
#include "Vector_traits.h"

namespace dare::utils {

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
    typedef T InternalType;
    typedef std::random_access_iterator_tag iterator_tag;
    typedef std::ptrdiff_t iterator_diff_type;

    /*! \struct Iterator
     * \brief random access forward iterator
     */
    struct Iterator {
        using iterator_category = iterator_tag;
        using difference_type = iterator_diff_type;
        using value_type = T;
        using pointer = T*;
        using reference = T&;

        /*!
         * \brief constructor
         * @param _ptr Pointer to data
         */
        explicit Iterator(pointer _ptr) : ptr(_ptr) {}
        reference operator*() const {
            return *ptr;
        }

        /*!
         * \brief arrow access operator
         */
        pointer operator->() {
            return ptr;
        }

        /*!
         * \brief prefix increment operator
         */
        Iterator& operator++() {
            ptr++;
            return *this;
        }

        /*!
         * \brief postfix increment operator
         */
        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        /*!
         * \brief prefix decrement operator
         */
        Iterator& operator--() {
            ptr--;
            return *this;
        }

        /*!
         * \brief postfix increment operator
         */
        Iterator operator--(int) {
            Iterator tmp = *this;
            --(*this);
            return tmp;
        }

        /*!
         * \brief equal comparison operator
         */
        bool operator==(const Iterator& a) const {
            return ptr == a.ptr;
        }

        /*!
         * \brief non-equal comparison operator
         */
        bool operator!=(const Iterator& a) const {
            return ptr != a.ptr;
        }

        bool operator<=(const Iterator& a) const {
            return ptr <= a.ptr;
        }

        bool operator>=(const Iterator& a) const {
            return ptr >= a.ptr;
        }

        bool operator<(const Iterator& a) const {
            return ptr < a.ptr;
        }

        bool operator>(const Iterator& a) const {
            return ptr > a.ptr;
        }

        /*!
         * \brief random increment operator
         */
        Iterator operator+(int n) const {
            return Iterator(ptr + n);
        }

        /*!
         * \brief random decrement operator
         */
        Iterator operator-(int n) const {
            return Iterator(ptr - n);
        }

        /*!
         * \brief distance operator
         */
        int operator-(const Iterator& other) {
            return ptr - other.ptr;
        }

        /*!
         * \brief relative access operator
         */
        reference operator[](int n) {
            return *(ptr + n);
        }

    private:
        pointer ptr;  // pointer to data
    };

    /*! \struct Iterator
     * \brief random access constant forward iterator
     */
    struct ConstIterator {
        using iterator_category = iterator_tag;
        using difference_type = iterator_diff_type;
        using value_type = T;
        using pointer = const T*;
        using reference = const T&;
        explicit ConstIterator(pointer _ptr) : ptr(_ptr) {}
        reference operator*() const {
            return *ptr;
        }
        pointer operator->() {
            return ptr;
        }
        // Prefix increment
        ConstIterator& operator++() {
            ptr++;
            return *this;
        }

        // Postfix increment
        ConstIterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        // Prefix decrement
        ConstIterator& operator--() {
            ptr--;
            return *this;
        }

        // Postfix decrement
        ConstIterator operator--(int) {
            ConstIterator tmp = *this;
            --(*this);
            return tmp;
        }

        bool operator==(const ConstIterator& a) {
            return ptr == a.ptr;
        }

        bool operator!=(const ConstIterator& a) {
            return ptr != a.ptr;
        }

        bool operator<=(const ConstIterator& a) {
            return ptr <= a.ptr;
        }

        bool operator>=(const ConstIterator& a) const {
            return ptr >= a.ptr;
        }

        bool operator<(const ConstIterator& a) const {
            return ptr < a.ptr;
        }

        bool operator>(const ConstIterator& a) const {
            return ptr > a.ptr;
        }

        ConstIterator operator+(int n) const {
            return ConstIterator(ptr + n);
        }

        ConstIterator operator-(int n) const {
            return ConstIterator(ptr - n);
        }

        int operator-(const ConstIterator& other) {
            return ptr - other.ptr;
        }

        reference operator[](int n) {
            return *(ptr + n);
        }

    private:
        pointer ptr;  // pointer to data
    };

    /*!
     * \brief reverse iterator
     */
    struct ReverseIterator {
        using iterator_category = iterator_tag;
        using difference_type = iterator_diff_type;
        using value_type = T;
        using pointer = T*;
        using reference = T&;
        explicit ReverseIterator(pointer _ptr) : ptr(_ptr) {}
        reference operator*() const {
            return *ptr;
        }
        pointer operator->() {
            return ptr;
        }
        // Prefix increment
        ReverseIterator& operator++() {
            ptr--;
            return *this;
        }

        // Postfix increment
        ReverseIterator operator++(int) {
            ReverseIterator tmp = *this;
            --(*this);
            return tmp;
        }

        // Prefix decrement
        ReverseIterator& operator--() {
            ptr++;
            return *this;
        }

        // Postfix decrement
        ReverseIterator operator--(int) {
            ReverseIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const ReverseIterator& a) {
            return ptr == a.ptr;
        }

        bool operator!=(const ReverseIterator& a) {
            return ptr != a.ptr;
        }

        bool operator<=(const ReverseIterator& a) {
            return ptr >= a.ptr;
        }

        bool operator>=(const ReverseIterator& a) const {
            return ptr <= a.ptr;
        }

        bool operator<(const ReverseIterator& a) const {
            return ptr > a.ptr;
        }

        bool operator>(const ReverseIterator& a) const {
            return ptr < a.ptr;
        }

        ReverseIterator operator+(int n) const {
            return ReverseIterator(ptr - n);
        }

        ReverseIterator operator-(int n) const {
            return ReverseIterator(ptr + n);
        }

        int operator-(const ReverseIterator& other) {
            return other.ptr - ptr;
        }

        reference operator[](int n) {
            return *(ptr - n);
        }

    private:
        pointer ptr;  // pointer to data
    };

    /*!
     * \brief constant reverse iterator
     */
    struct ConstReverseIterator {
        using iterator_category = iterator_tag;
        using difference_type = iterator_diff_type;
        using value_type = T;
        using pointer = const T*;
        using reference = const T&;
        explicit ConstReverseIterator(pointer _ptr) : ptr(_ptr) {}
        reference operator*() const {
            return *ptr;
        }
        pointer operator->() {
            return ptr;
        }
        // Prefix increment
        ConstReverseIterator& operator++() {
            ptr--;
            return *this;
        }

        // Postfix increment
        ConstReverseIterator operator++(int) {
            ConstReverseIterator tmp = *this;
            --(*this);
            return tmp;
        }

        // Prefix decrement
        ConstReverseIterator& operator--() {
            ptr++;
            return *this;
        }

        // Postfix decrement
        ConstReverseIterator operator--(int) {
            ConstReverseIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const ConstReverseIterator& a) {
            return ptr == a.ptr;
        }

        bool operator!=(const ConstReverseIterator& a) {
            return ptr != a.ptr;
        }

        bool operator<=(const ConstReverseIterator& a) {
            return ptr >= a.ptr;
        }

        bool operator>=(const ConstReverseIterator& a) const {
            return ptr <= a.ptr;
        }

        bool operator<(const ConstReverseIterator& a) const {
            return ptr > a.ptr;
        }

        bool operator>(const ConstReverseIterator& a) const {
            return ptr < a.ptr;
        }

        ConstReverseIterator operator+(int n) const {
            return ConstReverseIterator(ptr - n);
        }

        ConstReverseIterator operator-(int n) const {
            return ConstReverseIterator(ptr + n);
        }

        int operator-(const ConstReverseIterator& other) {
            return other.ptr - ptr;
        }

        reference operator[](int n) {
            return *(ptr - n);
        }

    private:
        pointer ptr;  // pointer to data
    };

    /*!
     * \brief constructor
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
     * \brief copy constructor
     * @param other object to copy from
     */
    template <typename A,
              typename = std::enable_if_t<std::is_convertible_v<A, T>> >
    Vector(const Vector<N, A>& other);

    /*!
     * \brief provides direct access to the data
     */
    T* data();

    /*!
     * \brief provides direct access to the data
     */
    const T* data() const;

    /*!
     * \brief assignment constructor
     * @param other object to copy from
     */
    template <typename A,
              typename = std::enable_if_t<std::is_convertible_v<A, T>>>
    Vector<N, T>& operator=(const Vector<N, A>& other);

    /*!
     * \brief access operator
     * @param n position to access
     * \note without NDEBUG, a bounds check will be conducted
     */
    T& operator[](std::size_t n);

    /*!
     * \brief access operator
     * @param n position to access
     * \note without NDEBUG, a bounds check will be conducted
     */
    const T& operator[](std::size_t n) const;

    /*!
     * \brief addition of other vector
     * @param other addition partner
     */
    Vector<N, T> operator+(const Vector<N, T>& other) const;

    /*!
     * \brief addition of a single value to all internal values
     * @param other addition partner
     */
    Vector<N, T> operator+(const T& val) const;

    /*!
     * \brief addition of another vector to the current instance
     * @param other addition partner
     */
    void operator+=(const Vector<N, T>& other);

    /*!
     * \brief addition of a single value to all internal values
     * @param other addition partner
     */
    void operator+=(const T& val);

    /*!
     * \brief subtraction of other vector
     * @param other addition partner
     */
    Vector<N, T> operator-(const Vector<N, T>& other) const;

    /*!
     * \brief subtraction of single value from all components
     * @param val value to subtract
     */
    Vector<N, T> operator-(const T& val) const;

    /*!
     * \brief -= operator
     * @param other other vector
     */
    void operator-=(const Vector<N, T>& other);

    /*!
     * \brief subtraction of single value from all components
     * @param val value to subtract
     */
    void operator-=(const T& val);

    /*!
     * \brief multiplication
     * @param other vector to mulitply with
     */
    Vector<N, T> operator*(const Vector<N, T>& other) const;

    /*!
     * \brief elementwise multiplication
     * @param val value to multiply with
     */
    Vector<N, T> operator*(const T& val) const;

    /*!
     * \brief multiplication
     * @param other vector for multiplication
     */
    void operator*=(const Vector<N, T>& other);

    /*!
     * \brief elementwise multiplication
     * @param val value to multiply with
     */
    void operator*=(const T& val);

    /*!
     * \brief division
     * @param other vector to divide with
     */
    Vector<N, T> operator/(const Vector<N, T>& other) const;

    /*!
     * \brief elementwise division
     * @param val value to divide with
     */
    Vector<N, T> operator/(const T& val) const;

    /*!
     * \brief division
     * @param other vector to divide with
     */
    void operator/=(const Vector<N, T>& other);

    /*!
     * \brief division
     * @param val value to divide with
     */
    void operator/=(const T& val);

    /*!
     * \brief comparison operator
     * @param other vector to compare with
     */
    bool operator==(const Vector<N, T>& other) const;

    /*!
     * \brief non-equal operator
     * @param other vector to compare with
     */
    bool operator!=(const Vector<N, T>& other) const;

    /*!
     * \brief returns number of elements
     */
    constexpr std::size_t size() const;

    /*!
     * \brief returns length of the vector
     * \note only sensible for floating point data, but also enabled for other types
     */
    T length() const;

    /*!
     * \brief return forward iterator to first object
     */
    Iterator begin();

    /*!
     * \brief returns constant forward iterator
     */
    ConstIterator begin() const;

    /*!
     * \brief returns constant forward iterator
     */
    ConstIterator cbegin() const;

    /*!
     * \brief return reverse iterator to first object
     */
    ReverseIterator rbegin();

    /*!
     * \brief returns constant reverse iterator
     */
    ConstReverseIterator rbegin() const;

    /*!
     * \brief returns constant reverse iterator
     */
    ConstReverseIterator crbegin() const;

    /*!
     * \brief returns iterator to end of data
     */
    Iterator end();

    /*!
     * \brief constant iterator to end of data
     */
    ConstIterator end() const;

    /*!
     * \brief constant iterator to end of data
     */
    ConstIterator cend() const;

    /*!
     * \brief returns iterator to end of data
     */
    ReverseIterator rend();

    /*!
     * \brief constant iterator to end of data
     */
    ConstReverseIterator rend() const;

    /*!
     * \brief constant iterator to end of data
     */
    ConstReverseIterator crend() const;

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
     * \brief sets all values to default values
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
     * \brief computes dot product with another vector
     * @param other vector to compute the dot product with
     */
    T dot(const Vector<N, T>& other) const;

    /*!
     * \brief computes cross product
     * @param other vector to compute the cross product with
     * \note only enabled for 3D vectors, for higher dimensions a more general algorithm is required
     */
    template <std::size_t Ns = N>
    typename std::enable_if<(Ns == 3), Vector<N, T>>::type cross(const Vector<N, T>& other) const;

    /*!
     * \brief outputs the vector
     */
    template <typename OS>
    friend OS& operator<<(OS& os, const Vector<N, T>& v) {
        os  << v[0];
        for (std::size_t i{1}; i < N; ++i)
            os << ' ' << v[i];

        return os;
    }

    /*
     * \brief calculates hash for hash-maps
     */
    std::size_t GetHash() const;

private:
    /*!
     * \brief iterates over the internal values and executes arbitrary manipulation
     * @param lambda operation to execute per data entry
     * @param op join operation for reduction
     */
    template <std::size_t I = 0, typename Expr, typename Op>
    auto IterateValues(Expr lambda, Op op);

    /*!
     * \brief iterates over the internal values and executes arbitrary manipulation
     * @param lambda operation to execute per data entry
     * @param op join operation for reduction
     */
    template <std::size_t I = 0, typename Expr, typename Op>
    auto IterateValues(Expr lambda, Op op) const;

    T _data[N];  //!< internal data set
};

typedef Vector<3, double> Vector3;  //!< convenient typedef

}  // namespace dare::utils

namespace std {
/*
 * \brief specialization of the hash-function for Vector
 */
template <std::size_t N, typename T>
class hash<dare::utils::Vector<N, T>> {
public:
    /*
     * \brief returns the hash for a Vector
     */
    std::size_t operator()(const dare::utils::Vector<N, T>& v) const {
        return v.GetHash();
    }
};

}  // namespace std

#include "Vector.inl"
#endif  // UTILITIES_VECTOR_H_
