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

namespace dare::Data {

template <typename T>
Stencil<T>::Stencil(std::string identifier, Ordinal stencil_size) : values(identifier, stencil_size) {
}

template <typename T>
Stencil<T>::Stencil(const Stencil<T>& other) {
    *this = other;
}

template <typename T>
Stencil<T>& Stencil<T>::operator=(const Stencil<T>& other) {
    if (&other == this)
        return *this;
    Resize(other.GetSize());
    Kokkos::deep_copy(values, other.values);
    return *this;
}

template <typename T>
void Stencil<T>::Resize(Ordinal stencil_size) {
    Kokkos::resize(values, stencil_size);
}

template <typename T>
typename Stencil<T>::Ordinal Stencil<T>::GetSize() const {
    return values.size();
}

template <typename T>
void Stencil<T>::InsertValue(Ordinal id, const T& value) {
    values[id] = value;
}

template <typename T>
T& Stencil<T>::GetValue(Ordinal pos) {
    return values[pos];
}

template <typename T>
const T& Stencil<T>::GetValue(Ordinal pos) const {
    return values[pos];
}

template <typename T>
const typename Stencil<T>::ArrayType& Stencil<T>::GetData() const {
    return values;
}

}  // namespace dare::Data
