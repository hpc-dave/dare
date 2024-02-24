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

#ifndef IO_VTKPXMLSTRUCTUREDGRIDWRITER_H_
#define IO_VTKPXMLSTRUCTUREDGRIDWRITER_H_

#include <vtkStructuredGrid.h>

#include <iostream>

#include <vector>
#include <string>

#include "MPI/ExecutionManager.h"
#include "Utilities/Vector.h"
#include "VTKOptions.h"

namespace dare::io {

struct VTKXMLPStructuredGridCData {
    std::string name;
    int number_components{0};
    VTKOutputType output_type{VTKOutputType::CELL_DATA};
    VTKDataAgglomerateType data_type{VTKDataAgglomerateType::SCALARS};

    /*!
     * @brief fast check if the dataset is valid
     */
    bool IsValid() const {
        bool valid{true};
        valid &= !name.empty();
        valid &= number_components > 0;
        return valid;
    }
};

class VTKPXMLStructuredGridWriter {
public:
    explicit VTKPXMLStructuredGridWriter(mpi::ExecutionManager* exec_man);

    void AddComponents(const VTKXMLPStructuredGridCData& data);
    void SetTime(double time);
    void SetGhostLevel(int ghost_level);
    void SetFileName(const std::string& fname);
    void SetGlobalExtent(const VTKExtent& extent);
    void SetPieceFileName(const std::string& fname);
    void SetPieceExtent(const VTKExtent& extent);
    bool Write();

    const std::string& GetDefaultFileExtension() const;

private:
    void GatherPieceInformation();

    void WriteHeader(std::ostringstream& os);

    void AddComponents(std::ostringstream& os);

    void AddPieces(std::ostringstream& os);

    void AddTimeStamp(std::ostringstream& os);

    void AddClosing(std::ostringstream& os);

    // global properties
    mpi::ExecutionManager* exec_man;
    bool is_root;
    std::string file_name;
    std::string path_piece;
    const std::string file_extension;
    VTKExtent extent_global;
    VTKExtent extent_piece;
    int ghost_levels;
    double time;
    const std::string in;

    // properties of the pieces
    std::vector<int> extent_pieces;
    std::vector<std::string> path_pieces;

    // properties of the components
    std::vector<VTKXMLPStructuredGridCData> component_data;
};

}  // end namespace dare::io

#endif  // IO_VTKPXMLSTRUCTUREDGRIDWRITER_H_
