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

#include <cmath>
namespace dare::utils {

template <std::size_t N, typename T>
template <typename... Ts,
          typename A,
          typename B>
Vector<N, T>::Vector(const Ts&... args) : VectorDecorator<N, N, T>() {
    SetValues<0>(args...);
}

template <std::size_t N, typename T>
template <typename A, typename B>
Vector<N, T>::Vector(const Vector<N, A>& other) : VectorDecorator<N, N, T>() {
    for (std::size_t n{0}; n < N; n++)
        this->_data[n] = static_cast<T>(other[n]);
}

template <std::size_t N, typename T>
template <typename A, typename B>
Vector<N, T>& Vector<N, T>::operator=(const Vector<N, A>& other) {
    if constexpr (std::is_same_v<A, T>)
        if (&other == this)
            return *this;

    for (std::size_t n{0}; n < N; n++)
        this->_data[n] = static_cast<T>(other[n]);

    return *this;
}

template <std::size_t N, typename T>
T* Vector<N, T>::data() {
    return this->_data;
}

template <std::size_t N, typename T>
const T* Vector<N, T>::data() const {
    return this->_data;
}

template <std::size_t N, typename T>
T& Vector<N, T>::operator[](std::size_t n) {
#ifndef NDEBUG
    if (n >= N) {
        std::cerr << __func__ << ": access @" << n << " out of bounds, size is " << N << std::endl;
    }
#endif
    return this->_data[n];
}

template <std::size_t N, typename T>
const T& Vector<N, T>::operator[](std::size_t n) const {
#ifndef NDEBUG
    if (n >= N) {
        std::cerr << __func__ << ": access @" << n << " out of bounds, size is " << N << std::endl;
    }
#endif
    return this->_data[n];
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator+(const Vector<N, T>& other) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] += other[i];
    }
    return vec;
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator+(const T& val) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] += val;
    }
    return vec;
}

template <std::size_t N, typename T>
void Vector<N, T>::operator+=(const Vector<N, T>& other) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] += other[i];
    }
}

template <std::size_t N, typename T>
void Vector<N, T>::operator+=(const T& val) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] += val;
    }
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator-(const Vector<N, T>& other) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] -= other[i];
    }
    return vec;
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator-(const T& val) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] -= val;
    }
    return vec;
}

template <std::size_t N, typename T>
void Vector<N, T>::operator-=(const Vector<N, T>& other) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] -= other[i];
    }
}

template <std::size_t N, typename T>
void Vector<N, T>::operator-=(const T& val) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] -= val;
    }
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator*(const Vector<N, T>& other) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] *= other[i];
    }
    return vec;
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator*(const T& val) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] *= val;
    }
    return vec;
}

template <std::size_t N, typename T>
void Vector<N, T>::operator*=(const Vector<N, T>& other) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] *= other[i];
    }
}

template <std::size_t N, typename T>
void Vector<N, T>::operator*=(const T& val) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] *= val;
    }
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator/(const Vector<N, T>& other) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] /= other[i];
    }
    return vec;
}

template <std::size_t N, typename T>
Vector<N, T> Vector<N, T>::operator/(const T& val) const {
    Vector<N, T> vec(*this);
    for (std::size_t i{0}; i < N; ++i) {
        vec[i] /= val;
    }
    return vec;
}

template <std::size_t N, typename T>
void Vector<N, T>::operator/=(const Vector<N, T>& other) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] /= other[i];
    }
}

template <std::size_t N, typename T>
void Vector<N, T>::operator/=(const T& val) {
    for (std::size_t i{0}; i < N; ++i) {
        this->_data[i] /= val;
    }
}

template <std::size_t N, typename T>
bool Vector<N, T>::operator==(const Vector<N, T>& other) const {
    bool is_same{true};
    for (std::size_t n{0}; n < N; n++)
        is_same &= (this->_data[n] == other[n]);
    return is_same;
}

template <std::size_t N, typename T>
bool Vector<N, T>::operator!=(const Vector<N, T>& other) const {
    return !(other == *this);
}

template <std::size_t N, typename T>
constexpr std::size_t Vector<N, T>::size() const {
    return N;
}

template <std::size_t N, typename T>
T Vector<N, T>::length() const {
    auto lambda = [this](const std::size_t i) {
        return this->_data[i] * this->_data[i];
    };

    auto op = [](const T& val1, const T& val2) {
        return val1 + val2;
    };

    return sqrt(IterateValues(lambda, op));
}

template <std::size_t N, typename T>
typename Vector<N, T>::Iterator Vector<N, T>::begin() {
    return Iterator(&this->_data[0]);
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstIterator Vector<N, T>::begin() const {
    return cbegin();
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstIterator Vector<N, T>::cbegin() const {
    return ConstIterator(&this->_data[0]);
}

template <std::size_t N, typename T>
typename Vector<N, T>::ReverseIterator Vector<N, T>::rbegin() {
    return ReverseIterator(&this->_data[N - 1]);
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstReverseIterator Vector<N, T>::rbegin() const {
    return cbegin();
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstReverseIterator Vector<N, T>::crbegin() const {
    return ConstReverseIterator(&this->_data[N - 1]);
}

template <std::size_t N, typename T>
typename Vector<N, T>::Iterator Vector<N, T>::end() {
    return Iterator(&this->_data[N]);
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstIterator Vector<N, T>::end() const {
    return cend();
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstIterator Vector<N, T>::cend() const {
    return ConstIterator(&this->_data[N]);
}

template <std::size_t N, typename T>
typename Vector<N, T>::ReverseIterator Vector<N, T>::rend() {
    return ReverseIterator(&this->_data[0])--;
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstReverseIterator Vector<N, T>::rend() const {
    return cend();
}

template <std::size_t N, typename T>
typename Vector<N, T>::ConstReverseIterator Vector<N, T>::crend() const {
    return ConstReverseIterator((&this->_data[0])--);
}

template <std::size_t N, typename T>
T Vector<N, T>::dot(const Vector<N, T>& other) const {
    T res{static_cast<T>(0)};
    for (std::size_t n{0}; n < N; ++n)
        res += this->_data[n] * other[n];
    return res;
}

template <std::size_t N, typename T>
template <std::size_t Ns>
typename std::enable_if<(Ns == 3), Vector<N, T>>::type Vector<N, T>::cross(const Vector<N, T>& other) const {
    Vector<N, T> vec(this->_data[1] * other[2] - this->_data[2] * other[1],
                     this->_data[2] * other[0] - this->_data[0] * other[2],
                     this->_data[0] * other[1] - this->_data[1] * other[0]);
    return vec;
}

template <std::size_t N, typename T>
std::size_t Vector<N, T>::GetHash() const {
    // computes hash using a variant of the Fowler-Noll-Vu hash function
    // details can be found on:
    // https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
    // or in the doxygen documentation of _FNVparam

    size_t result{_FNVparam<sizeof(std::size_t)>::offset};      // FNV offset
    const size_t prime{_FNVparam<sizeof(std::size_t)>::prime};  // FNV prime
    const char* ptr{reinterpret_cast<const char*>(this->_data)};      // ptr to the beginning of the data
    for (size_t count{0}; count < (sizeof(T) * N); ++count)
        result = (result * prime) ^ (*(ptr + count));  // (hash * FNV-prime) XOR (byte_of_data)

    return result;
}

template <std::size_t N, typename T>
template <std::size_t I,
          typename Tin,
          typename... Ts,
          typename A,
          typename B,
          typename C>
void Vector<N, T>::SetValues(const Tin& arg, const Ts&... args) {
    this->_data[I] = static_cast<T>(arg);
    if constexpr (sizeof...(args) > 0 && (I + 1 < N))
        SetValues<I + 1>(args...);
    else if constexpr (I + 1 < N)
        SetValues<I + 1>(static_cast<T>(0));
}

template <std::size_t N, typename T>
template <std::size_t I>
void Vector<N, T>::SetValues() {
    SetValues<0>(static_cast<T>(0));
}

template <std::size_t N, typename T>
template <typename A, typename B>
void Vector<N, T>::SetAllValues(const A& val) {
    for (std::size_t n{0}; n < N; n++)
        this->_data[n] = static_cast<T>(val);
}

template <std::size_t N, typename T>
template <std::size_t I, typename Expr, typename Op>
auto Vector<N, T>::IterateValues(Expr lambda, Op op) {
    if constexpr ((I + 1) < N) {
        return op(lambda(I), IterateValues<I + 1>(lambda, op));
    } else {
        return lambda(I);
    }
}

template <std::size_t N, typename T>
template <std::size_t I, typename Expr, typename Op>
auto Vector<N, T>::IterateValues(Expr lambda, Op op) const {
    if constexpr ((I + 1) < N) {
        return op(lambda(I), IterateValues<I + 1>(lambda, op));
    } else {
        return lambda(I);
    }
}

}  // namespace dare::utils
