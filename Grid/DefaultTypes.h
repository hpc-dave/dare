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

#ifndef GRID_DEFAULTTYPES_H_
#define GRID_DEFAULTTYPES_H_
#include <Tpetra_Map_decl.hpp>

namespace dare::Grid::details {
using LocalOrdinalType = Tpetra::Map<>::local_ordinal_type;
using GlobalOrdinalType = Tpetra::Map<>::global_ordinal_type;
using ExecutionSpace = Tpetra::Map<>::execution_space;
using MemorySpace = Tpetra::Map<>::memory_space;
using NodeType = Tpetra::Map<>::node_type;
}  // namespace dare::Grid::details

namespace dare::defaults {
using LocalOrdinalType = Grid::details::LocalOrdinalType;
using GlobalOrdinalType = Grid::details::GlobalOrdinalType;
using ScalarType = double;
}  // end namespace dare::defaults

#endif  // GRID_DEFAULTTYPES_H_
