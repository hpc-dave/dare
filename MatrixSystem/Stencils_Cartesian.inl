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

namespace dare::Data {

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::CenterMatrixStencil() {
    static_assert(static_cast<char>(Positions::CENTER) == 0, "The center position needs to be at 0!");
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::~CenterMatrixStencil() {
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::CenterMatrixStencil(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other)
    : coefficients(other.coefficients) {
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) {
    if (this == &other)
        return *this;

    coefficients = other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator*=(SC v) {
    coefficients *= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator*(SC v) const {
    CenterMatrixStencil<GridType> s(*this);
    s *= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator/=(SC v) {
    coefficients /= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator/(SC v) const {
    CenterMatrixStencil<GridType> s(*this);
    s /= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator+=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) {
    coefficients += other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator+(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) const {
    CenterMatrixStencil<GridType> s(*this);
    s += other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator-=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) {
    coefficients -= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator-(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) const {
    CenterMatrixStencil<GridType> s(*this);
    s -= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator*=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) {
    coefficients *= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator*(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) const {
    CenterMatrixStencil<GridType> s(*this);
    s *= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator/=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) {
    coefficients /= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::operator/(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>& other) const {
    CenterMatrixStencil<GridType> s(*this);
    s /= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::SetValue(Positions pos, SC v) {
    RangeCheck(__func__, pos);
    coefficients[static_cast<char>(pos)] = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::SetAll(SC v) {
    for (auto& e : coefficients)
        e = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
SC& CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::GetValue(Positions pos) {
    RangeCheck(__func__, pos);
    return coefficients[static_cast<char>(pos)];
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
SC CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::GetValue(Positions pos) const {
    RangeCheck(__func__, pos);
    return coefficients[static_cast<char>(pos)];
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
typename CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::DataArray&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::GetData() {
    return coefficients;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
const typename CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::DataArray&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::GetData() const {
    return coefficients;
}

template <std::size_t Dim, typename LO, typename GO, typename SC>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>::RangeCheck(std::string func, Positions pos) const {
#ifndef DARE_NDEBUG
    if (static_cast<char>(pos) >= NUM_ENTRIES) {
        std::cerr << "In " << func << ": Provided position (" << std::to_string(static_cast<char>(pos)) << ") "
                  << "is out of range of " << std::to_string(NUM_ENTRIES) << "-point center based stencil\n";
    }
#endif
}

}  // end namespace dare::Data
