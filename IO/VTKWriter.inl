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

namespace dare::io {

template<typename Grid>
VTKWriter<Grid>::VTKWriter(mpi::ExecutionManager* ex_man)
    : exec_man(ex_man) {
}

template <typename Grid>
template <typename... PairLike>
bool VTKWriter<Grid>::Write(const std::string& base_path,
                            const std::string& parallel_folder_name,
                            const PairLike&... data) {
    static_assert(sizeof...(PairLike) > 0, "at least one data set has to be provided for writing!");
    using LO = typename Grid::LocalOrdinalType;
    const std::string proc_id = std::to_string(exec_man->GetRank());
    std::string parallel_data_path(base_path);
    if (!parallel_data_path.empty()) {
        if (base_path.back() != '/')
            parallel_data_path += '/';
    }
    parallel_data_path += parallel_folder_name;
    if (parallel_data_path.back() != '/')
        parallel_data_path += '/';

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

        // allocate vtkGrid
        vtkNew<GridType> vtkDataSet;
        vtkDataSet = VTKOptions<Grid>::GetGrid(*grep);
        int num_cells = vtkDataSet->GetNumberOfCells();

        if (num_cells != grep->GetNumberLocalCells()) {
            ERROR << "Number of cells are incompatible! VTK computed "
                  << num_cells << " vs the number of local cells: "
                  << grep->GetNumberLocalCells() << ERROR_CLOSE;
        }

        // add data per data set
        for (auto num_instance : grid.second) {
            // access vtk output type
            vtkNew<vtkDoubleArray> data_array;
            int num_components{0};
            std::vector<std::string> cnames;
            VTKOutputType otype{VTKOutputType::SCALAR_DATA};
            // allocate field according to the specified type
            auto SetInformation = [&](std::size_t pos, auto instance) {
                if (pos == num_instance) {
                    data_array->SetName(instance.second->GetName().c_str());
                    num_components = static_cast<int>(instance.second->GetNumComponents());
                    cnames.resize(num_components);
                    for (int n{0}; n < num_components; n++) {
                        cnames[n] = instance.second->GetComponentName(n);
                    }
                }
            };
            LoopThroughData<0>(SetInformation, data_tuple);
            data_array->SetNumberOfComponents(num_components);
            data_array->SetNumberOfTuples(vtkDataSet->GetNumberOfCells());
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
                    for (LO cell_id{0}; cell_id < num_cells; cell_id++) {
                        int vtk_id = VTKOptions<Grid>::Map(*grep, cell_id);
                        // tuple
                        for (int n{0}; n < num_components; n++) {
                            tuple_like[n] = instance.second->At(cell_id, n);
                        }
                        data_array->SetTuple(vtk_id, tuple_like.data());
                    }
                }
            };
            LoopThroughData<0>(SetData, data_tuple);

            // and add to data set
            switch (otype) {
            case VTKOutputType::SCALAR_DATA:
                vtkDataSet->GetCellData()->SetScalars(data_array);
                break;
            case VTKOutputType::VECTORS:
                vtkDataSet->GetCellData()->SetVectors(data_array);
                break;
            default:
                ERROR << "Unsupported output type provided!" << ERROR_CLOSE;
            }
        }  // end loop through instances

        // write to file
        vtkNew<WriterType> writer;
        std::ostringstream os;

        os << parallel_data_path
             << grep->GetName()
             << "." << writer->GetDefaultFileExtension();
        exec_man->Print(dare::mpi::Verbosity::Medium) << "Writing to file "
                << parallel_data_path << "*" << grep->GetName()
                << "." << writer->GetDefaultFileExtension() << std::endl;
        writer->SetFileName(os.str().c_str());
        writer->SetInputData(vtkDataSet);
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

}  // end namespace dare::io
