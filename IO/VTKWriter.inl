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
struct Args {
    vtkProgrammableFilter* pf;
    int local_extent[6];
};

// function to operate on the point attribute data
void inline execute(void* arg) {
    Args* args = reinterpret_cast<Args*>(arg);
    // auto info = args->pf->GetOutputInformation(0);
    auto output_tmp = args->pf->GetOutput();
    auto input_tmp = args->pf->GetInput();
    vtkStructuredGrid* output = dynamic_cast<vtkStructuredGrid*>(output_tmp);
    vtkStructuredGrid* input = dynamic_cast<vtkStructuredGrid*>(input_tmp);
    output->ShallowCopy(input);
    output->SetExtent(args->local_extent);
}

template<typename Grid>
VTKWriter<Grid>::VTKWriter(mpi::ExecutionManager* ex_man, double _time, int _step)
    : exec_man(ex_man), time(_time), step(_step) {
}

template <typename Grid>
template <typename... PairLike>
bool VTKWriter<Grid>::Write(const std::string& base_path,
                            const std::string& parallel_folder_name,
                            const PairLike&... data) {
    static_assert(sizeof...(PairLike) > 0, "at least one data set has to be provided for writing!");
    using LO = typename Grid::LocalOrdinalType;
    using vtkOrdinal = typename VTKOptions<Grid>::vtkOrdinal;
    std::string parallel_data_path = details::VTKGetParallelOutputPath(*exec_man,
                                                                       base_path,
                                                                       parallel_folder_name);

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
        // get a representation for this grid
        std::string grid_name = grid.first;
        const typename Grid::Representation* grep{representations[grid_name]};

        // some data for the root file writing
        // std::list<VTKXMLPStructuredGridCData> c_data;

        // allocate vtkGrid
        vtkNew<GridType> vtkDataSet;
        vtkDataSet = VTKOptions<Grid>::GetGrid(*grep);
        vtkOrdinal num_cells_glob = vtkDataSet->GetNumberOfCells();
        if (num_cells_glob != grep->GetNumberGlobalCells()) {
            ERROR << "Number of cells are incompatible! VTK computed "
                  << num_cells_glob << " vs the number of global cells: "
                  << grep->GetNumberGlobalCells() << ERROR_CLOSE;
        }
        LO num_cells_loc = grep->GetNumberLocalCells();

        // add time stamp, if time is negative this will be skipped
        AddTimeStamp(vtkDataSet);

        // add data per data set
        for (auto num_instance : grid.second) {
            // access vtk output type
            vtkNew<vtkDoubleArray> data_array;
            int num_components{0};
            std::string data_name;
            std::vector<std::string> cnames;
            VTKDataAgglomerateType otype{VTKDataAgglomerateType::SCALARS};
            // allocate field according to the specified type
            auto SetInformation = [&](std::size_t pos, auto instance) {
                if (pos == num_instance) {
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
                        vtkOrdinal vtk_id = VTKOptions<Grid>::Map(*grep, cell_id);
                        // tuple
                        for (int n{0}; n < num_components; n++) {
                            tuple_like[n] = instance.second->At(cell_id, n);
                        }
                        data_array->SetTuple(vtk_id, tuple_like.data());
                    }
                }
            };
            LoopThroughData<0>(SetData, data_tuple);
            vtkDataSet->GetCellData()->AddArray(data_array);

            // // and add to data set
            switch (otype) {
            case VTKDataAgglomerateType::SCALARS:
                // vtkDataSet->GetCellData()->SetScalars(data_array);
                break;
            case VTKDataAgglomerateType::VECTORS:
                // vtkDataSet->GetCellData()->SetVectors(data_array);
                break;
            default:
                ERROR << "Unsupported output type provided!" << ERROR_CLOSE;
            }

            // VTKXMLPStructuredGridCData c_data_loc;
            // c_data_loc.name = data_name;
            // c_data_loc.number_components = cnames.size();
            // c_data_loc.output_type = VTKOutputType::CELL_DATA;
            // c_data_loc.data_type = otype;
            // c_data.push_back(c_data_loc);
        }  // end loop through instances

        // write to file
        // vtkNew<WriterType> writer;
        // std::string fname = details::VTKGetOutputFileName(
        //                         exec_man, parallel_data_path,
        //                         grep->GetName(), step, writer->GetDefaultFileExtension());
        // writer->SetFileName(fname.c_str());
        // writer->SetInputData(vtkDataSet);
        // if (!writer->Write()) {
        //     ERROR << "VTK writer returned with an error!" << ERROR_CLOSE;
        //     return false;
        // }
        // write root file --- IMPROVE THIS!
        VTKPXMLStructuredGridWriter pwriter(exec_man);
        std::string root_path = details::VTKGetParallelOutputFileName(exec_man, base_path,
                                                                      grep->GetName(),
                                                                      step,
                                                                      pwriter.GetDefaultFileExtension());
        // std::string piece_path = details::VTKGetOutputFileName(exec_man,
        //                                                        details::VTKGetParallelOutputPath(*exec_man,
        //                                                        "", parallel_folder_name),
        //                                                        grep->GetName(),
        //                                                        step,
        //                                                        writer->GetDefaultFileExtension());

        // VTKExtent extent_pglobal = VTKOptions<Grid>::GetPExtentGlobal(*grep);
        VTKExtent extent_plocal  = VTKOptions<Grid>::GetPExtentLocal(*grep);
        // pwriter.SetFileName(root_path);
        // pwriter.SetGhostLevel(grep->GetNumberGhostCells());
        // pwriter.SetGlobalExtent(extent_pglobal);
        // pwriter.SetPieceExtent(extent_plocal);
        // pwriter.SetPieceFileName(piece_path);
        // pwriter.SetTime(time);
        // for (auto it = c_data.begin(); it != c_data.end(); it++) {
        //     pwriter.AddComponents(*it);
        // }
        // pwriter.Write();
        // Create a vtkProgrammableFilter
        vtkNew<vtkProgrammableFilter> pf;
        vtkNew<vtkMPIController> contr;
        // contr->Initialize(&argc, &argv, 1);
        // contr->Initialize();
        // Initialize an instance of Args
        Args args;
        args.pf = pf;
        for (int i = 0; i < 6; ++i)
            args.local_extent[i] = extent_plocal[i];

        pf->SetExecuteMethod(execute, &args);

        // Create a structured grid and assign point data and cell data to it
        // auto structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        // structuredGrid->SetExtent(global_extent);
        pf->SetInputData(vtkDataSet);
        // structuredGrid->SetPoints(points);
        // structuredGrid->GetCellData()->AddArray(density);

        // Create the parallel writer and call some functions
        auto parallel_writer = vtkSmartPointer<vtkXMLPStructuredGridWriter>::New();
        parallel_writer->SetInputConnection(pf->GetOutputPort());
        parallel_writer->SetController(contr);
        parallel_writer->SetFileName(root_path.c_str());
        parallel_writer->SetNumberOfPieces(exec_man->GetNumberProcesses());
        parallel_writer->SetStartPiece(exec_man->GetRank());
        parallel_writer->SetEndPiece(exec_man->GetRank());
        parallel_writer->SetDataModeToBinary();
        parallel_writer->SetGhostLevel(grep->GetNumberGhostCells());
        parallel_writer->Update();
        parallel_writer->Write();
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
    std::string pfolder_name = std::to_string(exec_man->GetRank());
    if (!fman.CreateDirectory(fman.GetOutputPath().native() + '/' + pfolder_name, true)) {
        ERROR << "An error occured while creating the folder for paralle output!" << ERROR_CLOSE;
        return false;
    }
    return Write(fman.GetOutputPath().native(), pfolder_name, data...);
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
void VTKWriter<Grid>::AddTimeStamp(vtkNew<GridType>& data_set) {    // NOLINT
    if (time < 0)
        return;

    vtkNew<vtkDoubleArray> time_array;
    time_array->SetName("TimeValue");
    time_array->SetNumberOfTuples(1);
    time_array->SetTuple1(0, time);
    data_set->GetFieldData()->AddArray(time_array);
}

}  // end namespace dare::io
