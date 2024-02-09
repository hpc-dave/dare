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

#ifndef IO_VTKWRITER_H_
#define IO_VTKWRITER_H_

#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkNew.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <map>
#include <string>
#include <tuple>

#include "MPI/ExecutionManager.h"
#include "VTKOptions.h"

namespace dare::io {

template<typename VTKGridType>
struct VTKWriterMapper {
};

template <>
struct VTKWriterMapper<vtkStructuredGrid> {
    using type = vtkXMLStructuredGridWriter;
};

template <>
struct VTKWriterMapper<vtkUnstructuredGrid> {
    using type = vtkXMLUnstructuredGridWriter;
};

template <typename Grid>
class VTKWriter {
public:
    using GridType = typename VTKOptions<Grid>::GridType;
    using WriterType = typename VTKWriterMapper<GridType>::type;

    explicit VTKWriter(mpi::ExecutionManager* ex_man);

    template<typename... Data>
    bool Write(const std::string& base_path, const Data&... data);

private:
    template<std::size_t I, typename Lambda, typename... Data>
    void LoopThroughData(Lambda lambda, std::tuple<const Data&...> data);

    mpi::ExecutionManager* exec_man;
};

}  // end namespace dare::io

#include "VTKWriter.inl"

#endif  // IO_VTKWRITER_H_
