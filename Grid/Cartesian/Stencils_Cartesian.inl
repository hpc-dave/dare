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

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::CenterMatrixStencil() {
    static_assert(static_cast<char>(Positions::CENTER) == 0, "The center position needs to be at 0!");
    SetAll(0);
    rhs.SetAllValues(0);
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::~CenterMatrixStencil() {
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::CenterMatrixStencil(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other)
    : coefficients(other.coefficients), rhs(other.rhs) {
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    if (this == &other)
        return *this;

    coefficients = other.coefficients;
    rhs = other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*=(SC v) {
    for (auto& c : coefficients)
        c *= v;
    for (auto& e : rhs)
        e *= v;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*(SC v) const {
    CenterMatrixStencil<GridType, SC, N> s(*this);
    s *= v;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/=(SC v) {
    for (auto& c : coefficients)
        c /= v;
    for (auto& e : rhs)
        e /= v;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/(SC v) const {
    CenterMatrixStencil<GridType, SC, N> s(*this);
    s /= v;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator+=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients += other.coefficients;
    rhs += other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator+(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    CenterMatrixStencil<GridType, SC, N> s(*this);
    s += other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator-=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients -= other.coefficients;
    rhs -= other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator-(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    CenterMatrixStencil<GridType, SC, N> s(*this);
    s -= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients *= other.coefficients;
    rhs += other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    CenterMatrixStencil<GridType, SC, N> s(*this);
    s *= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/=(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients /= other.coefficients;
    rhs += other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/(
    const CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    CenterMatrixStencil<GridType, SC, N> s(*this);
    s /= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
SC& CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::Center(std::size_t n) {
    return GetValue(Positions::CENTER, n);
}

template <std::size_t Dim, typename SC, std::size_t N>
SC CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::Center(std::size_t n) const {
    return GetValue(Positions::CENTER, n);
}

template <std::size_t Dim, typename SC, std::size_t N>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValue(Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients[n][Grid::ToNum(pos)] = v;
}

template <std::size_t Dim, typename SC, std::size_t N>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetRHS(
    std::size_t n, SC v) {
    rhs[n] = v;
}

template <std::size_t Dim, typename SC, std::size_t N>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetAll(SC v) {
    for (auto& a : coefficients)
        for (auto& e : a)
            e = v;
}

template <std::size_t Dim, typename SC, std::size_t N>
SC& CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValue(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients[n][Grid::ToNum(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
SC CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValue(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients[n][Grid::ToNum(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
dare::utils::Vector<N, SC>
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValues(Positions pos) const {
    dare::utils::Vector<N, SC> v;
    for (std::size_t n{0}; n < N; n++) {
        v[n] = GetValue(pos, n);
    }
    return v;
}

template <std::size_t Dim, typename SC, std::size_t N>
typename CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetData() {
    return coefficients;
}

template <std::size_t Dim, typename SC, std::size_t N>
const typename CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetData() const {
    return coefficients;
}

template <std::size_t Dim, typename SC, std::size_t N>
typename CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::RHSType&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS() {
    return rhs;
}

template <std::size_t Dim, typename SC, std::size_t N>
const typename CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::RHSType&
CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS() const {
    return rhs;
}

template <std::size_t Dim, typename SC, std::size_t N>
SC& CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS(std::size_t n) {
    return rhs[n];
}

template <std::size_t Dim, typename SC, std::size_t N>
SC CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS(std::size_t n) const {
    return rhs[n];
}

template <std::size_t Dim, typename SC, std::size_t N>
void CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::RangeCheck(
    std::string func, Positions pos, std::size_t n) const {
#ifndef DARE_NDEBUG
    if (Grid::ToNum(pos) >= static_cast<char>(NUM_ENTRIES)) {
        ERROR << "In " << func << ": Provided position (" << std::to_string(Grid::ToNum(pos)) << ") "
                  << "is out of range of " << std::to_string(NUM_ENTRIES) << "-point center based stencil"
                  << ERROR_CLOSE;
    }
    if (n >= NUM_COMPONENTS) {
        ERROR << "In " << func << ": Provided component ID (" << std::to_string(n) << ") "
                  << "is out of range of " << std::to_string(NUM_COMPONENTS) << " components"
                  << ERROR_CLOSE;
    }
#endif
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::FaceMatrixStencil() {
    static_assert(static_cast<char>(Positions::CENTER) == 0, "The center position needs to be at 0!");
    SetAll(0);
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::~FaceMatrixStencil() {
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::FaceMatrixStencil(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other)
    : coefficients_nb(other.coefficients_nb), coefficients_c(other.coefficients_c), rhs(other.rhs) {
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    if (this == &other)
        return *this;

    coefficients_nb = other.coefficients_nb;
    coefficients_c = other.coefficients_c;
    rhs = other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*=(SC v) {
    for (auto& c : coefficients_nb)
        c *= v;
    for (auto& c : coefficients_c)
        c *= v;
    rhs *= v;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*(SC v) const {
    FaceMatrixStencil<GridType, SC, N> s(*this);
    s *= v;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/=(SC v) {
    for (auto& c : coefficients_c)
        c /= v;
    for (auto& c : coefficients_nb)
        c /= v;
    rhs /= v;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/(SC v) const {
    FaceMatrixStencil<GridType, SC, N> s(*this);
    s /= v;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator+=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients_nb += other.coefficients_nb;
    coefficients_c += other.coefficients_c;
    rhs += other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator+(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceMatrixStencil<GridType, SC, N> s(*this);
    s += other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator-=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients_c -= other.coefficients_c;
    coefficients_nb -= other.coefficients_nb;
    rhs -= other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator-(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceMatrixStencil<GridType, SC, N> s(*this);
    s -= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients_c *= other.coefficients_c;
    coefficients_nb *= other.coefficients_nb;
    rhs += other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceMatrixStencil<GridType, SC, N> s(*this);
    s *= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/=(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients_c /= other.coefficients_c;
    coefficients_nb /= other.coefficients_nb;
    rhs += other.rhs;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/(
    const FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceMatrixStencil<GridType, SC, N> s(*this);
    s /= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValueNeighbor(
    Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients_nb[n][Grid::ToFace(pos)] = v;
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValueCenter(Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients_c[n][Grid::ToFace(pos)] = v;
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValues(
    Positions pos, std::size_t n, SC v_nb, SC v_c) {
    SetValueNeighbor(pos, n, v_nb);
    SetValueCenter(pos, n, v_c);
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetRHS(
    Positions pos, std::size_t n, SC v) {
    rhs.SetValue(pos, n, v);
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValueNeighbor(
    Positions pos, const dare::utils::Vector<N, SC>& v) {
    for (std::size_t n{0}; n < N; n++)
        SetValueNeighbor(pos, n, v[n]);
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValueCenter(
    Positions pos, const dare::utils::Vector<N, SC>& v) {
    for (std::size_t n{0}; n < N; n++)
        SetValueCenter(pos, n, v[n]);
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValues(
    Positions pos, const dare::utils::Vector<N, SC>& v_nb, const dare::utils::Vector<N, SC>& v_c) {
    SetValueNeighbor(pos, v_nb);
    SetValueCenter(pos, v_c);
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetRHS(
    Positions pos, const dare::utils::Vector<N, SC>& v) {
    for (std::size_t n{0}; n < N; n++)
        SetRHS(pos, n, v[n]);
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetAll(SC v) {
    for (auto& a : coefficients_c)
        for (auto& e : a)
            e = v;
    for (auto& a : coefficients_nb)
        for (auto& e : a)
            e = v;
}

template <std::size_t Dim, typename SC, std::size_t N>
SC& FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValueCenter(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients_c[n][Grid::ToFace(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
SC& FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValueNeighbor(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients_nb[n][Grid::ToFace(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
SC FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValueCenter(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients_c[n][Grid::ToFace(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
SC FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValueNeighbor(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients_nb[n][Grid::ToFace(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
SC& FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return rhs.GetValue(pos, n);
}

template <std::size_t Dim, typename SC, std::size_t N>
SC FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return rhs.GetValue(pos, n);
}

template <std::size_t Dim, typename SC, std::size_t N>
typename FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetDataNeighbor() {
    return coefficients_nb;
}

template <std::size_t Dim, typename SC, std::size_t N>
typename FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetDataCenter() {
    return coefficients_c;
}

template <std::size_t Dim, typename SC, std::size_t N>
const typename FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetDataNeighbor() const {
    return coefficients_nb;
}

template <std::size_t Dim, typename SC, std::size_t N>
const typename FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetDataCenter() const {
    return coefficients_c;
}

template <std::size_t Dim, typename SC, std::size_t N>
typename FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::RHSType&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS() {
    return rhs;
}

template <std::size_t Dim, typename SC, std::size_t N>
const typename FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::RHSType&
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetRHS() const {
    return rhs;
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>::RangeCheck(
    std::string func, Positions pos, std::size_t n) const {
#ifndef DARE_NDEBUG
    if (Grid::ToNum(pos) > static_cast<char>(NUM_FACES)) {
        ERROR << "In " << func << ": Provided position (" << std::to_string(Grid::ToNum(pos)) << ") "
                  << "is out of range of " << std::to_string(GridType::STENCIL_SIZE) << "-point center based stencil"
                  << ERROR_CLOSE;
    } else if (Grid::ToNum(pos) == 0) {
        ERROR << "In " << func << ": The face stencil cannot take CENTER as an argument!" << ERROR_CLOSE;
    }
    if (n >= NUM_COMPONENTS) {
        ERROR << "In " << func << ": Provided component ID (" << std::to_string(n) << ") "
                  << "is out of range of " << std::to_string(NUM_COMPONENTS) << " components" << ERROR_CLOSE;
    }
#endif
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::FaceValueStencil() {
    static_assert(static_cast<char>(Positions::CENTER) == 0, "The center position needs to be at 0!");
    SetAll(0);
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::~FaceValueStencil() {
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::FaceValueStencil(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other)
    : coefficients(other.coefficients) {
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    if (this == &other)
        return *this;

    coefficients = other.coefficients;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*=(SC v) {
    for (auto& c : coefficients)
        c *= v;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*(SC v) const {
    FaceValueStencil<GridType, SC, N> s(*this);
    s *= v;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/=(SC v) {
    for (auto& c : coefficients)
        c /= v;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/(SC v) const {
    FaceValueStencil<GridType, SC, N> s(*this);
    s /= v;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator+=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients += other.coefficients;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator+(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceValueStencil<GridType, SC, N> s(*this);
    s += other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator-=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients -= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator-(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceValueStencil<GridType, SC, N> s(*this);
    s -= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients *= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceValueStencil<GridType, SC, N> s(*this);
    s *= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/=(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) {
    coefficients /= other.coefficients;
    return *this;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator/(
    const FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>& other) const {
    FaceValueStencil<GridType, SC, N> s(*this);
    s /= other;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::operator*(
    const FaceMatrixStencil<GridType, SC, N>& other) const {
    FaceMatrixStencil<GridType, SC, N> s(other);
    s.GetDataCenter() *= coefficients;
    s.GetDataNeighbor() *= coefficients;
    return s;
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValue(Positions pos, std::size_t n, SC v) {
    RangeCheck(__func__, pos, n);
    coefficients[n][Grid::ToFace(pos)] = v;
}
template <std::size_t Dim, typename SC, std::size_t N>
void FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetValues(Positions pos,
                                                                    const dare::utils::Vector<N, SC>& v) {
    for (std::size_t n{0}; n < N; n++) {
        RangeCheck(__func__, pos, n);
        coefficients[n][Grid::ToFace(pos)] = v[n];
    }
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::SetAll(SC v) {
    for (auto& a : coefficients)
        for (auto& e : a)
            e = v;
}

template <std::size_t Dim, typename SC, std::size_t N>
SC& FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValue(Positions pos, std::size_t n) {
    RangeCheck(__func__, pos, n);
    return coefficients[n][Grid::ToFace(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
SC FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetValue(Positions pos, std::size_t n) const {
    RangeCheck(__func__, pos, n);
    return coefficients[n][Grid::ToFace(pos)];
}

template <std::size_t Dim, typename SC, std::size_t N>
typename FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetData() {
    return coefficients;
}

template <std::size_t Dim, typename SC, std::size_t N>
const typename FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::DataArray&
FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::GetData() const {
    return coefficients;
}

template <std::size_t Dim, typename SC, std::size_t N>
void FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>::RangeCheck(
    std::string func, Positions pos, std::size_t n) const {
#ifndef DARE_NDEBUG
    if (Grid::ToNum(pos) > static_cast<char>(NUM_FACES)) {
        ERROR << "In " << __FILE__ << ": "<< func << ": Provided position ("
                  << std::to_string(Grid::ToNum(pos)) << ") "
                  << "is out of range of " << std::to_string(NUM_FACES + 1)
                  << "-point center based stencil" << ERROR_CLOSE;
    } else if (Grid::ToNum(pos) == 0) {
        ERROR << "In " << func << ": The face stencil cannot take CENTER as an argument!"
              << ERROR_CLOSE;
    }
    if (n >= NUM_COMPONENTS) {
        ERROR << "In " << func << ": Provided component ID (" << std::to_string(n) << ") "
                  << "is out of range of " << std::to_string(NUM_COMPONENTS) << " components"
                  << ERROR_CLOSE;
    }
#endif
}
}  // end namespace dare::Data
