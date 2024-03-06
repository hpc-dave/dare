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

#include "VTKPXMLStructuredGridWriter.h"
#include <vtkInformation.h>
#include <vtkMPIController.h>
#include <vtkProgrammableFilter.h>
#include <vtkXMLPStructuredGridWriter.h>

namespace dare::io {

template<typename Grid>
VTKWriter<Grid>::VTKWriter(mpi::ExecutionManager* ex_man, double _time, int _step)
    : exec_man(ex_man), time(_time), step(_step) {
}

template <typename Grid>
template <typename... PairLike>
bool VTKWriter<Grid>::Write(const std::string& base_path,
                            const PairLike&... data) {
    static_assert(sizeof...(PairLike) > 0, "at least one data set has to be provided for writing!");
    using LO = typename Grid::LocalOrdinalType;
    using Options = VTKOptions<Grid>;
    using vtkOrdinal = typename Options::vtkOrdinal;
    using Mapper = typename Options::Mapper;

    /*
     * To do:
     * - move data assignment into subfunction
     * - test in 3D
     */
    // convert parameter pack to tuple to use STL-functionality
    auto data_tuple = std::forward_as_tuple(data...);

    // Loop through provided data to group all instances with the same underlying grid
    // As an example, this allows to reuse the scalar grid for all scalar transport quantities,
    // while the staggered fields need dedicated handling
    std::map<std::string, std::list<std::size_t>> grouped_data;
    std::map<std::string, const typename Grid::Representation*> representations;
    auto GroupData = [&](std::size_t pos, auto instance) {
        std::string grid_name = instance.second->GetGridRepresentation().GetName();
        grouped_data[grid_name].push_back(pos);
        if (representations.find(grid_name) == representations.end())
            representations[grid_name] = &instance.second->GetGridRepresentation();
    };
    LoopThroughData<0>(GroupData, data_tuple);

    // Loop through each determined grid and write the associated data
    // to the file
    for (auto grid : grouped_data) {
        // basic setup for a grid
        // get a representation for this grid, the total number of cells/entries to print
        // and a mapper for local grid IDs to vtk grid IDs
        std::string grid_name = grid.first;
        const typename Grid::Representation* grep{representations[grid_name]};
        LO num_cells_loc = grep->GetNumberLocalCells();
        Mapper mapToVTK(grep);

        // allocate vtkGrid
        vtkNew<GridType> vtkDataSet;
        bool success = Options::AllocateGrid(*grep, vtkDataSet);
        if (!success) {
            exec_man->Print(dare::mpi::Verbosity::Low)
                << "Error during grid allocation, cannot write " << grep->GetName() << " to file!" << std::endl;
            continue;
        }

        // add time stamp, if time is negative this will be skipped
        AddTimeStamp(vtkDataSet);

        // add data per data set
        for (auto num_instance : grid.second) {
            // access vtk output type
            vtkNew<vtkDoubleArray> data_array;
            int num_components{0};
            std::string data_name;
            std::vector<std::string> cnames;
            // VTKDataAgglomerateType otype{VTKDataAgglomerateType::SCALARS};
            VTKOutputType otype{VTKOutputType::CELL_DATA};
            // allocate field according to the specified type
            auto SetInformation = [&](std::size_t pos, auto instance) {
                if (pos == num_instance) {
                    otype = instance.first;
                    data_name = instance.second->GetName();
                    num_components = static_cast<int>(instance.second->GetNumComponents());
                    cnames.resize(num_components);
                    for (int n{0}; n < num_components; n++) {
                        cnames[n] = instance.second->GetComponentName(n);
                    }
                }
            };
            LoopThroughData<0>(SetInformation, data_tuple);
            data_array->SetName(data_name.c_str());
            data_array->SetNumberOfComponents(num_components);
            data_array->SetNumberOfTuples(num_cells_loc);
            if (num_components > 1) {
                for (int i{0}; i < num_components; i++) {
                    std::string cname = cnames[i];
                    if (cname.empty())
                        cname = std::to_string(i);
                    data_array->SetComponentName(i, cname.c_str());
                }
            }
            auto SetData = [&](std::size_t pos, auto instance) {
                if (pos == num_instance) {
                    std::vector<double> tuple_like(num_components);
                    for (LO cell_id{0}; cell_id < num_cells_loc; cell_id++) {
                        // vtkOrdinal vtk_id = VTKOptions<Grid>::Map(*grep, cell_id);
                        vtkOrdinal vtk_id = mapToVTK(cell_id);
                        // tuple
                        for (int n{0}; n < num_components; n++) {
                            tuple_like[n] = instance.second->At(cell_id, n);
                        }
                        data_array->SetTuple(vtk_id, tuple_like.data());
                    }
                }
            };
            LoopThroughData<0>(SetData, data_tuple);
            switch (otype) {
            case VTKOutputType::CELL_DATA:
                vtkDataSet->GetCellData()->AddArray(data_array);
                break;
            case VTKOutputType::POINT_DATA:
                vtkDataSet->GetPointData()->AddArray(data_array);
                break;
            }
        }  // end loop through instances

        // Create the parallel writer and call some functions
        vtkNew<Writer> writer;
        vtkNew<vtkMPIController> contr;
        std::string root_path = details::VTKGetParallelOutputFileName(exec_man, base_path,
                                                                      grep->GetName(),
                                                                      step,
                                                                      writer->GetDefaultFileExtension());
        writer->SetController(contr);
        auto persistent_obj = Options::AddDataToWriter(*grep, exec_man, vtkDataSet, writer);
        writer->SetFileName(root_path.c_str());
        writer->SetUseSubdirectory(true);
        writer->SetDataModeToBinary();
        writer->Update();
        writer->Write();
    }
    return true;
}

template <typename Grid>
template <typename... PairLike>
bool VTKWriter<Grid>::Write(const FileSystemManager& fman, const PairLike&... data) {
    if (!fman.CreateCommonDirectory(fman.GetOutputPath()) && exec_man->AmIRoot()) {
        ERROR << "Could not create the output folder!" << ERROR_CLOSE;
        return false;
    }

    return Write(fman.GetOutputPath().native(), data...);
}

template <typename Grid>
template <std::size_t I, typename Lambda, typename... Data>
void VTKWriter<Grid>::LoopThroughData(Lambda lambda, std::tuple<const Data&...> data) {
    lambda(I, std::get<I>(data));
    if constexpr (I < (sizeof...(Data) - 1)) {
        LoopThroughData<I + 1>(lambda, data);
    }
}

template <typename Grid>
void VTKWriter<Grid>::AddTimeStamp(GridType* data_set) {    // NOLINT
    if (time < 0)
        return;

    vtkNew<vtkDoubleArray> time_array;
    time_array->SetName("TimeValue");
    time_array->SetNumberOfTuples(1);
    time_array->SetTuple1(0, time);
    data_set->GetFieldData()->AddArray(time_array);
}

}  // end namespace dare::io
