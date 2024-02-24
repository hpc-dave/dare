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
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPStructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <bit>  // for c++20 std::endian, once available
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

template<typename VTKGridType>
struct VTKWriterMapper {
};

template <>
struct VTKWriterMapper<vtkStructuredGrid> {
    using type = vtkXMLStructuredGridWriter;
    using ptype = vtkXMLPStructuredGridWriter;
};

template <>
struct VTKWriterMapper<vtkUnstructuredGrid> {
    using type = vtkXMLUnstructuredGridWriter;
    using ptype = vtkXMLPUnstructuredGridWriter;
};

namespace details {
/*!
 * @brief concatenates and checks output folder path for parallel output
 * @param exman execution manager
 * @param base_path basic output path, where all output is directed to
 * @param pfolder_name name of the parallel folder
 * @return correct string for parallel subfolders
 */
[[nodiscard]] std::string VTKGetParallelOutputPath(const dare::mpi::ExecutionManager& exman,
                                                   const std::string& base_path,
                                                   const std::string& pfolder_name);

/*!
 * @brief concatenates specific file name for output
 * @param exman execution manager
 * @param parallel_data_path path to parallel folder
 * @param grid_name name of the grid
 * @param step step to print
 * @param ext extension of the file
 */
[[nodiscard]] std::string VTKGetOutputFileName(dare::mpi::ExecutionManager* exman,
                                         const std::string& parallel_data_path,
                                         const std::string& grid_name,
                                         int step,
                                         const std::string& ext);
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

/*!
 * @brief writes root file for parallel vtk format, overload for structured grid
 * @param path path of root file
 * @param extent_domain total extent of the domain
 * @param ghost_level number of ghost cells of total domain
 * @param path_distributed_files paths to all the distributed files
 * @param comp_data global data of each component on the grid
 * @param extent_subdomains extent of all subdomains
 * @param time time stamp
 */
bool VTKWritePXMLRootFile(const std::string& path,
                          const VTKExtent& extent_domain,
                          int ghost_level,
                          const std::list<std::string>& path_distributed_files,
                          const std::list<VTKXMLPStructuredGridComponentData>& comp_data,
                          const std::list<VTKExtent>& extent_subdomains,
                          double time = -1.);
}  // end namespace details

template <typename Grid>
class VTKWriter {
public:
    using GridType = typename VTKOptions<Grid>::GridType;
    using WriterType = typename VTKWriterMapper<GridType>::type;

    explicit VTKWriter(mpi::ExecutionManager* ex_man,
                       double time,
                       int step);

    template<typename... Data>
    bool Write(const std::string& base_path,
               const std::string& parallel_folder_name,
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
