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

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::GridVector() : GridVector("not_specified", 0, GridRepresentation()) {
}

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::GridVector(std::string identifier, GridRepresentation _grid)
    : GridVector(identifier, _grid.GetNumberLocalCells(), _grid) {
}

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::GridVector(std::string identifier, LO num_cells, GridRepresentation _grid)
    : ident_string(identifier),
      grid(_grid),
      data(identifier, num_cells * N) {
}

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::~GridVector() {}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::Initialize(std::string identifier, GridRepresentation _grid) {
    ident_string = identifier;
    grid = _grid;
    data = DualViewType(identifier, _grid.GetNumberLocalCells() * N);
}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::Resize(LO n) {
    data.resize(n);
}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::ResizeByGrid(LO n) {
    Resize(n * N);
}

template <typename Grid, typename T, std::size_t N>
T& GridVector<Grid, T, N>::At(std::size_t n) {
    return operator[](n);
}

template <typename Grid, typename T, std::size_t N>
T& GridVector<Grid, T, N>::At(LO n, std::size_t c) {
    return operator[](n * N + c);
}

template <typename Grid, typename T, std::size_t N>
T& GridVector<Grid, T, N>::At(const Index& ind, std::size_t c) {
    return At(grid.MapIndexToOrdinalLocal(ind), c);
}

template <typename Grid, typename T, std::size_t N>
T GridVector<Grid, T, N>::At(std::size_t n) const {
    return operator[](n);
}

template <typename Grid, typename T, std::size_t N>
T GridVector<Grid, T, N>::At(LO n, std::size_t c) const {
    return operator[](n * N + c);
}

template <typename Grid, typename T, std::size_t N>
T GridVector<Grid, T, N>::At(const Index& ind, std::size_t c) const {
    return At(grid.MapIndexToOrdinalLocal(ind), c);
}

template <typename Grid, typename T, std::size_t N>
dare::utils::Vector<N, T> GridVector<Grid, T, N>::GetValues(const LO n) const {
    dare::utils::Vector<N, T> values;
    LO n_start{n * N};
    for (std::size_t c{0}; c < N; c++) {
        values[c] = operator[](n_start + c);
    }
    return values;
}

template <typename Grid, typename T, std::size_t N>
dare::utils::Vector<N, T> GridVector<Grid, T, N>::GetValues(const Index& ind) const {
    return GetValues(grid.MapIndexToOrdinalLocal(ind));
}

template <typename Grid, typename T, std::size_t N>
T& GridVector<Grid, T, N>::operator[](LO n) {
    return data.h_view[n];
}

template <typename Grid, typename T, std::size_t N>
T GridVector<Grid, T, N>::operator[](LO n) const {
    return data.h_view[n];
}

template <typename Grid, typename T, std::size_t N>
std::size_t GridVector<Grid, T, N>::GetSize() const {
    return data.h_view.size();
}

template <typename Grid, typename T, std::size_t N>
template <typename TargetSpace>
void GridVector<Grid, T, N>::Synchronize() {
    if constexpr (std::is_same_v<TargetSpace, HostSpace>) {
        data.template modify<ExecutionSpace>();
        data.template sync<HostSpace>();
    } else if (std::is_same_v<TargetSpace, ExecutionSpace>) {
        data.template modify<HostSpace>();
        data.template sync<ExecutionSpace>();
    }
}

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N> GridVector<Grid, T, N>::GetDeepCopy() const {
    GridVector<Grid, T, N> other(ident_string, data.h_view.size(), grid);
    GetDeepCopy(&other);
    return other;
}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::GetDeepCopy(GridVector<Grid, T, N>* other) const {
    other->ident_string = ident_string;
    other->grid = grid;
    other->data.resize(GetSize());
    Kokkos::deep_copy(other->data, data);
}

template <typename Grid, typename T, std::size_t N>
typename GridVector<Grid, T, N>::DeviceViewType& GridVector<Grid, T, N>::GetDeviceView() {
    return data.d_view;
}

template <typename Grid, typename T, std::size_t N>
const typename GridVector<Grid, T, N>::DeviceViewType& GridVector<Grid, T, N>::GetDeviceView() const {
    return data.h_view;
}

template <typename Grid, typename T, std::size_t N>
const typename GridVector<Grid, T, N>::GridRepresentation& GridVector<Grid, T, N>::GetGridRepresentation() const {
    return grid;
}

template <typename Grid, typename T, std::size_t N>
typename GridVector<Grid, T, N>::GridRepresentation& GridVector<Grid, T, N>::GetGridRepresentation() {
    return grid;
}

template <typename Grid, typename T, std::size_t N>
std::string GridVector<Grid, T, N>::GetIdentifier() const {
    return ident_string;
}

}  // namespace dare::Data
