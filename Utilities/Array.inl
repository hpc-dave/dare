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

namespace dare::utils {

template <std::size_t M, std::size_t N, typename T>
typename Array<M, N, T>::InnerVec& Array<M, N, T>::operator[](std::size_t m) {
    return Get(m);
}

template <std::size_t M, std::size_t N, typename T>
const typename Array<M, N, T>::InnerVec& Array<M, N, T>::operator[](std::size_t m) const {
    return Get(m);
}

template <std::size_t M, std::size_t N, typename T>
typename Array<M, N, T>::InnerVec& Array<M, N, T>::Get(std::size_t m) {
#ifndef NDEBUG
    if (m >= M) {
        ERROR << "access out of bounds (" << m << " >= " << M << ERROR_CLOSE;
    }
#endif
    return data[m];
}

template <std::size_t M, std::size_t N, typename T>
const typename Array<M, N, T>::InnerVec& Array<M, N, T>::Get(std::size_t m) const {
#ifndef NDEBUG
    if (m >= M) {
        ERROR << "access out of bounds (" << m << " >= " << M << ERROR_CLOSE;
    }
#endif
    return data[m];
}

template <std::size_t M, std::size_t N, typename T>
T& Array<M, N, T>::Get(std::size_t m, std::size_t n) {
#ifndef NDEBUG
    if (m >= M) {
        ERROR << "access out of bounds (" << m << " >= " << M << ERROR_CLOSE;
    }
#endif
    return data[m][n];
}

template <std::size_t M, std::size_t N, typename T>
const T& Array<M, N, T>::Get(std::size_t m, std::size_t n) const {
#ifndef NDEBUG
    if (m >= M) {
        ERROR << "access out of bounds (" << m << " >= " << M << ERROR_CLOSE;
    }
#endif
    return data[m][n];
}

}  // namespace dare::utils
