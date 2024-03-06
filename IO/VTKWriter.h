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

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>

#include <list>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <iomanip>

#include "FileSystemManager.h"
#include "MPI/ExecutionManager.h"
#include "VTKOptions.h"
#include "Utilities/Vector.h"

namespace dare::io {

namespace details {

/*!
 * @brief concatenates file name for parallel format
 * @param exman execution manager
 * @param output_path path to parallel folder
 * @param grid_name name of the grid
 * @param step step to print
 * @param ext extension of the file
 * Note, that VTK stores the data per process and requires one file to coordinate those.
 * This if the filename for the coordinated one!
 */
[[nodiscard]] std::string VTKGetParallelOutputFileName(dare::mpi::ExecutionManager* exman,
                                                       const std::string& parallel_data_path,
                                                       const std::string& grid_name,
                                                       int step,
                                                       const std::string& ext);

struct VTKXMLPStructuredGridComponentData {
    std::string Name;
    std::string NumberOfComponents;
    std::string OutputType;
    std::string DataAgglomerateType;
};

}  // end namespace details

template <typename Grid>
class VTKWriter {
public:
    using GridType = typename VTKOptions<Grid>::GridType;
    using Writer   = typename VTKWriterMapper<GridType>::type;

    explicit VTKWriter(mpi::ExecutionManager* ex_man,
                       double time,
                       int step);

    template<typename... Data>
    bool Write(const std::string& base_path,
               const Data&... data);

    template <typename... Data>
    bool Write(const FileSystemManager& fman, const Data&... data);

private:
    template<std::size_t I, typename Lambda, typename... Data>
    void LoopThroughData(Lambda lambda, std::tuple<const Data&...> data);

    void AddTimeStamp(vtkNew<GridType>& data_set);  // NOLINT

    mpi::ExecutionManager* exec_man;
    double time;
    double step;
};

}  // end namespace dare::io

#include "VTKWriter.inl"

#endif  // IO_VTKWRITER_H_
