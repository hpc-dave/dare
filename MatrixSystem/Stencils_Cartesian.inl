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

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::CenterMatrixStencil() {
    static_assert(static_cast<char>(Positions::CENTER) == 0, "The center position needs to be at 0!");
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::~CenterMatrixStencil() {
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::CenterMatrixStencil(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other)
    : coefficients(other.coefficients) {
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    if (this == &other)
        return *this;

    coefficients = other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*=(SC v) {
    for (auto& c : coefficients)
        c *= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*(SC v) const {
    CenterMatrixStencil<GridType, N> s(*this);
    s *= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/=(SC v) {
    for (auto& c : coefficients)
        c /= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/(SC v) const {
    CenterMatrixStencil<GridType, N> s(*this);
    s /= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator+=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients += other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator+(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    CenterMatrixStencil<GridType, N> s(*this);
    s += other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator-=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients -= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator-(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    CenterMatrixStencil<GridType, N> s(*this);
    s -= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients *= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    CenterMatrixStencil<GridType, N> s(*this);
    s *= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients /= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    CenterMatrixStencil<GridType, N> s(*this);
    s /= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetValue(Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients[n][static_cast<char>(pos)] = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetAll(SC v) {
    for (auto& a : coefficients)
        for (auto& e : a)
            e = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC& CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValue(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients[n][static_cast<char>(pos)];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValue(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients[n][static_cast<char>(pos)];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
typename CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetData() {
    return coefficients;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
const typename CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetData() const {
    return coefficients;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::RangeCheck(
    std::string func, Positions pos, std::size_t n) const {
#ifndef DARE_NDEBUG
    if (static_cast<char>(pos) >= NUM_ENTRIES) {
        std::cerr << "In " << func << ": Provided position (" << std::to_string(static_cast<char>(pos)) << ") "
                  << "is out of range of " << std::to_string(NUM_ENTRIES) << "-point center based stencil\n";
    }
    if (n >= NUM_COMPONENTS) {
        std::cerr << "In " << func << ": Provided component ID (" << std::to_string(n) << ") "
                  << "is out of range of " << std::to_string(NUM_COMPONENTS) << " components\n";
    }
#endif
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::FaceMatrixStencil() {
    static_assert(static_cast<char>(Positions::CENTER) == 0, "The center position needs to be at 0!");
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::~FaceMatrixStencil() {
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::FaceMatrixStencil(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other)
    : coefficients_nb(other.coefficients_nb), coefficients_c(other.coefficients_c) {
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    if (this == &other)
        return *this;

    coefficients_nb = other.coefficients_nb;
    coefficients_c = other.coefficients_c;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*=(SC v) {
    for (auto& c : coefficients_nb)
        c *= v;
    for (auto& c : coefficients_c)
        c *= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*(SC v) const {
    FaceMatrixStencil<GridType, N> s(*this);
    s *= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/=(SC v) {
    for (auto& c : coefficients_c)
        c /= v;
    for (auto& c : coefficients_nb)
        c /= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/(SC v) const {
    FaceMatrixStencil<GridType, N> s(*this);
    s /= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator+=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients_nb += other.coefficients_nb;
    coefficients_c += other.coefficients_c;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator+(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceMatrixStencil<GridType, N> s(*this);
    s += other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator-=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients_c -= other.coefficients_c;
    coefficients_nb -= other.coefficients_nb;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator-(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceMatrixStencil<GridType, N> s(*this);
    s -= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients_c *= other.coefficients_c;
    coefficients_nb *= other.coefficients_nb;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceMatrixStencil<GridType, N> s(*this);
    s *= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients_c /= other.coefficients_c;
    coefficients_nb /= other.coefficients_nb;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceMatrixStencil<GridType, N> s(*this);
    s /= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetValueNeighbor(
    Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients_nb[n][static_cast<char>(pos) - 1] = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetValueCenter(Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients_c[n][static_cast<char>(pos) - 1] = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetValues(
    Positions pos, std::size_t n, SC v_nb, SC v_c) {
    SetValueNeighbor(pos, n, v_nb);
    SetValueCenter(pos, n, v_c);
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetAll(SC v) {
    for (auto& a : coefficients_c)
        for (auto& e : a)
            e = v;
    for (auto& a : coefficients_nb)
        for (auto& e : a)
            e = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC& FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValueCenter(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients_c[n][static_cast<char>(pos) - 1];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC& FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValueNeighbor(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients_nb[n][static_cast<char>(pos) - 1];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValueCenter(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients_c[n][static_cast<char>(pos) - 1];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValueNeighbor(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients_nb[n][static_cast<char>(pos) - 1];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
typename FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetDataNeighbor() {
    return coefficients_nb;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
typename FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetDataCenter() {
    return coefficients_c;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
const typename FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetDataNeighbor() const {
    return coefficients_nb;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
const typename FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetDataCenter() const {
    return coefficients_c;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::RangeCheck(
    std::string func, Positions pos, std::size_t n) const {
#ifndef DARE_NDEBUG
    if (static_cast<char>(pos) > NUM_FACES) {
        std::cerr << "In " << func << ": Provided position (" << std::to_string(static_cast<char>(pos)) << ") "
                  << "is out of range of " << std::to_string(NUM_FACES+1) << "-point center based stencil\n";
    } else if (static_cast<char>(pos) == 0) {
        std::cerr << "In " << func << ": The face stencil cannot take CENTER as an argument!\n";
    }
    if (n >= NUM_COMPONENTS) {
        std::cerr << "In " << func << ": Provided component ID (" << std::to_string(n) << ") "
                  << "is out of range of " << std::to_string(NUM_COMPONENTS) << " components\n";
    }
#endif
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::FaceValueStencil() {
    static_assert(static_cast<char>(Positions::CENTER) == 0, "The center position needs to be at 0!");
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::~FaceValueStencil() {
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::FaceValueStencil(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other)
    : coefficients(other.coefficients) {
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    if (this == &other)
        return *this;

    coefficients = other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*=(SC v) {
    for (auto& c : coefficients)
        c *= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*(SC v) const {
    FaceValueStencil<GridType, N> s(*this);
    s *= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/=(SC v) {
    for (auto& c : coefficients)
        c /= v;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/(SC v) const {
    FaceValueStencil<GridType, N> s(*this);
    s /= v;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator+=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients += other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator+(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceValueStencil<GridType, N> s(*this);
    s += other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator-=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients -= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator-(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceValueStencil<GridType, N> s(*this);
    s -= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients *= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceValueStencil<GridType, N> s(*this);
    s *= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) {
    coefficients /= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator/(
    const FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>& other) const {
    FaceValueStencil<GridType, N> s(*this);
    s /= other;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::operator*(
    const FaceMatrixStencil<GridType, N>& other) const {
    FaceMatrixStencil<GridType, N> s(other);
    s.GetDataCenter() *= coefficients;
    s.GetDataNeighbor() *= coefficients;
    return s;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetValue(Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients[n][static_cast<char>(pos) - 1] = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::SetAll(SC v) {
    for (auto& a : coefficients)
        for (auto& e : a)
            e = v;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC& FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValue(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients[n][static_cast<char>(pos) - 1];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
SC FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetValue(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients[n][static_cast<char>(pos) - 1];
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
typename FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetData() {
    return coefficients;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
const typename FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::DataArray&
FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::GetData() const {
    return coefficients;
}

template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
void FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>::RangeCheck(
    std::string func, Positions pos, std::size_t n) const {
#ifndef DARE_NDEBUG
    if (static_cast<char>(pos) > NUM_FACES) {
        std::cerr << "In " << func << ": Provided position (" << std::to_string(static_cast<char>(pos)) << ") "
                  << "is out of range of " << std::to_string(NUM_FACES + 1) << "-point center based stencil\n";
    } else if (static_cast<char>(pos) == 0) {
        std::cerr << "In " << func << ": The face stencil cannot take CENTER as an argument!\n";
    }
    if (n >= NUM_COMPONENTS) {
        std::cerr << "In " << func << ": Provided component ID (" << std::to_string(n) << ") "
                  << "is out of range of " << std::to_string(NUM_COMPONENTS) << " components\n";
    }
#endif
}
}  // end namespace dare::Data
