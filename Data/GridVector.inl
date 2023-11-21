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

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::GridVector() : GridVector("not_specified", 0, GridRepresentation()) {
}

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::GridVector(std::string identifier, GridRepresentation _grid)
    : GridVector(identifier, _grid.GetNumberLocalCells(), _grid) {
}

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::GridVector(std::string identifier, LO num_cells, GridRepresentation _grid)
    : data(identifier, num_cells * N), grid(_grid) {
}

template <typename Grid, typename T, std::size_t N>
GridVector<Grid, T, N>::~GridVector() {}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::Resize(LO n) {
    Kokkos::resize(data, n);
    Kokkos::resize(data_h, n);
}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::ResizeByGrid(LO n) {
    Resize(n * N);
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
T GridVector<Grid, T, N>::At(LO n, std::size_t c) const {
    return operator[](n * N + c);
}

template <typename Grid, typename T, std::size_t N>
T GridVector<Grid, T, N>::At(const Index& ind, std::size_t c) const {
    return At(grid.MapIndexToOrdinalLocal(ind), c);
}

template <typename Grid, typename T, std::size_t N>
T& GridVector<Grid, T, N>::operator[](LO n) {
    return data[n];
}

template <typename Grid, typename T, std::size_t N>
T GridVector<Grid, T, N>::operator[](LO n) const {
    return data[n];
}

template <typename Grid, typename T, std::size_t N>
std::size_t GridVector<Grid, T, N>::GetSize() const {
    return data.size();
}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::CopyToHost() const {
    Kokkos::deep_copy(data_h, data);
}

template <typename Grid, typename T, std::size_t N>
void GridVector<Grid, T, N>::CopytoDevice() const {
    Kokkos::deep_copy(data, data_h);
}

}  // namespace dare::Data
