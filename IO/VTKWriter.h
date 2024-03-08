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
#include <vtkPointData.h>
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

}  // end namespace details

/*!
 * \brief high-level object for creating output in VTKformat
 * Each time output should be created, an instance of this class should be
 * constructed. Do not try to reuse it.
 */
template <typename Grid>
class VTKWriter {
public:
    using GridType = typename VTKOptions<Grid>::GridType;
    using Writer   = typename VTKWriterMapper<GridType>::type;
    template<typename T>
    using GroupedData = std::map<std::string, T>;
    using GroupedIndices = GroupedData<std::list<std::size_t>>;
    using GroupedRepresentations = GroupedData<typename Grid::Representation*>;
    using LO = typename Grid::LocalOrdinalType;
    using Options = VTKOptions<Grid>;
    using vtkOrdinal = typename Options::vtkOrdinal;
    using Mapper = typename Options::Mapper;

    /*!
     * @brief constructor
     * @param ex_man instance of execution manager
     * @param time current simulation time (ignored if <0)
     * @param step time/simulation step
     */
    explicit VTKWriter(mpi::ExecutionManager* ex_man,
                       double time,
                       int step);

    /*!
     * @brief writes data into a specified folder
     * @tparam ...Data type of parameter pack with data
     * @param base_path path to folder
     * @param ...data parameter pack with output data
     * @return true, if successful
     */
    template<typename... Data>
    bool Write(const std::string& base_path,
               const Data&... data);

    /*!
     * @brief writes data into folder according to file system manager
     * @tparam ...Data type of parameter pack with data
     * @param fman file system manager
     * @param ...data parameter pack with output data
     * @return true, if successful
     * The use of a FileSystemManager allows the dynamic piping of data
     * into dedicated folders. The FileSystemManager will also create the
     * output folder, if required.
     */
    template <typename... Data>
    bool Write(const FileSystemManager& fman, const Data&... data);

private:
    /*!
     * @brief Convenient loop to access an instance in the parameter pack
     * @tparam Lambda anonymous function/functor to work with the data
     * @tparam ...Data type parameter pack
     * @tparam I instance to access
     * @param lambda functor of form (std::size_t pos, auto instance)
     * @param data tuple with references to data
     */
    template<std::size_t I, typename Lambda, typename... Data>
    void LoopThroughData(Lambda lambda, std::tuple<const Data&...> data);

    /*!
     * @brief sort and group the data in the tuple
     * @param data tupe with references to the data
     * @param positions positions of the data sets of a grid
     * @param rep representations of the grid
     * Loop through provided data to group all instances with the same underlying grid
     * As an example, this allows to reuse the scalar grid for all scalar transport quantities,
     * while the staggered fields need dedicated handling
     */
    template<typename... Data>
    void GroupData(std::tuple<const Data&...> data,
                   GroupedIndices* positions,
                   GroupedRepresentations* rep);

    template <typename... Data>
    void PopulateVTKArray(std::tuple<const Data&...> data,
                          std::size_t num_instance,
                          LO num_cells_loc,
                          const Mapper& mapToVtk,
                          vtkDoubleArray* data_array);

    /*!
     * @brief adds a timestamp, if the time is >= 0
     * @param data_set vtkGridType
     */
    void AddTimeStamp(GridType* data_set);

    mpi::ExecutionManager* exec_man;    //!< reference to execution manager
    double time;                        //!< timestamp
    int step;                           //!< time/simulation step
};

}  // end namespace dare::io

#include "VTKWriter.inl"

#endif  // IO_VTKWRITER_H_
