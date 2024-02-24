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

#include <bit>  // for c++20 std::endian, once available
#include <fstream>
#include <iomanip>
#include <list>
#include <type_traits>

namespace dare::io {

VTKPXMLStructuredGridWriter::VTKPXMLStructuredGridWriter(mpi::ExecutionManager* ex_man)
    : exec_man(ex_man),
      is_root(true),
      file_name("output.pvts"), file_extension("pvts"),
      extent_global(-1, -1, -1, -1, -1, -1), ghost_levels(-1), time(-1.), in("  ") {
    if (!exec_man) {
        ERROR << "No ExecutionManager was provided, cannot collect output data!" << ERROR_CLOSE;
    } else {
        is_root = exec_man->AmIRoot();
    }
}

void VTKPXMLStructuredGridWriter::AddComponents(const VTKXMLPStructuredGridCData& data) {
    if (!is_root)
        return;
    if (!data.IsValid()) {
        ERROR << "Provided component data is invalid, won't add it to output!" << ERROR_CLOSE;
        return;
    }
    component_data.push_back(data);
}
void VTKPXMLStructuredGridWriter::SetTime(double t) {
    if (!is_root)
        return;
    time = t;
}

void VTKPXMLStructuredGridWriter::SetGhostLevel(int g_level) {
    if (!is_root)
        return;
    if (g_level < 0) {
        ERROR << "Provided ghost level (" << g_level << ") is invalid, should be >=0" << ERROR_CLOSE;
    }
    ghost_levels = g_level;
}

void VTKPXMLStructuredGridWriter::SetFileName(const std::string& fname) {
    if (!is_root)
        return;
    if (fname.empty()) {
        ERROR << "Empty file name provided! I will not update current file name " << file_name << ERROR_CLOSE;
        return;
    }
    file_name = fname;
}

void VTKPXMLStructuredGridWriter::SetGlobalExtent(const VTKExtent& extent) {
    extent_global = extent;
}

void VTKPXMLStructuredGridWriter::SetPieceFileName(const std::string& fname) {
    if (fname.empty()) {
        ERROR << "Name of piece is empty!" << ERROR_CLOSE;
    }
    path_piece = fname;
}

void VTKPXMLStructuredGridWriter::SetPieceExtent(const VTKExtent& extent) {
    extent_piece = extent;
}

bool VTKPXMLStructuredGridWriter::Write() {
    /*
     * Here the most important stuff happens. First the root process gathers the required information
     * of each piece of the other processes.
     * Then the different parts of the file are combined to a entity
     * and finally written to file
     */

    // Collect information of pieces
    GatherPieceInformation();

    // buffer information
    std::ostringstream os;

    // Start with header
    WriteHeader(os);

    // Add components
    AddComponents(os);

    // Add location of each piece
    AddPieces(os);

    // Add time stamp
    AddTimeStamp(os);

    // Add closure
    AddClosing(os);

    // Write to File
    std::ofstream ofs(file_name, std::ios::out);
    if (!ofs) {
        ERROR << "Could not open file: " << file_name << ERROR_CLOSE;
        return false;
    }
    ofs << os.str();
    return true;
}

const std::string& VTKPXMLStructuredGridWriter::GetDefaultFileExtension() const {
    return file_extension;
}

void VTKPXMLStructuredGridWriter::GatherPieceInformation() {
    const int tag_name{6001};
    const int tag_extent{6001};
    const int num_requests{2};
    int num_pieces = exec_man->GetNumberProcesses();
    int rank_root = exec_man->GetRankRoot();
    path_pieces.resize(num_pieces);
    extent_pieces.resize(num_pieces * 6);   // extents as one single array for use with MPI_Gather
    MPI_Request requests[num_requests];     // NOLINT

    // Collect all file names at root process
    exec_man->Isend(path_piece.c_str(), path_piece.size(), rank_root, tag_name, &requests[0]);
    if (is_root) {
        for (int n{0}; n < num_pieces; n++) {
            MPI_Status status;
            exec_man->Probe(n, tag_name, &status);
            int size_path{0};
            MPI_Get_count(&status, MPI_CHAR, &size_path);
            path_pieces[n].resize(size_path);
            exec_man->Recv(path_pieces[n].data(), size_path, n, tag_name);
        }
    }

    // Collect all extents at root process
    exec_man->Isend(extent_piece.data(), 6, rank_root, tag_extent, &requests[1]);
    if (is_root) {
        for (int n{0}; n < num_pieces; n++) {
            exec_man->Recv(&extent_pieces[n * 6], 6, n, tag_extent);
        }
    }

    exec_man->Waitall(num_requests, requests);
}

void VTKPXMLStructuredGridWriter::WriteHeader(std::ostringstream& os) {
    if (!is_root)
        return;

    std::string endianness = "LittleEndian";
    // substitute with the following for C++20
    // (std::endian::native == std::endian::little ? "LittleEndian" : "BigEndian");
    std::string extent_str = std::to_string(extent_global[0]);
    for (std::size_t n{1}; n < 6; n++) {
        extent_str += " ";
        extent_str += std::to_string(extent_global[n]);
    }

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"" << endianness << "\">\n";
    os << in << "<PStructuredGrid WholeExtent=\"" << extent_str << "\" GhostLevel=\"" << std::to_string(ghost_levels) << "\">\n";  // NOLINT
    os << in << in << "<PPoints>\n";
    os << in << in << in << "<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n";  // NOLINT
    os << in << in << in << "</PDataArray>\n";
    os << in << in << "</PPoints>\n";
}

void VTKPXMLStructuredGridWriter::AddComponents(std::ostringstream& os) {
    if (!is_root)
        return;

    enum vtkTypeAndAgglomerates {
        POINT_AND_SCALAR,
        POINT_AND_VECTOR,
        CELL_AND_SCALAR,
        CELL_AND_VECTOR,
        NUM_TYPES
    };
    std::array<std::string, 2> OutputStrings[NUM_TYPES];
    OutputStrings[POINT_AND_SCALAR] = {"PointData", "scalars"};
    OutputStrings[POINT_AND_VECTOR] = {"PointData", "vectors"};
    OutputStrings[CELL_AND_SCALAR] = {"CellData", "scalars"};
    OutputStrings[CELL_AND_VECTOR] = {"CellData", "vectors"};
    auto TestTypeAndForm = [](const VTKXMLPStructuredGridCData& c,
                              VTKOutputType otype, VTKDataAgglomerateType dtype) {
        return (c.output_type == otype) && (c.data_type == dtype);
    };

    std::list<VTKXMLPStructuredGridCData> datasets[NUM_TYPES];
    for (auto it = component_data.begin(); it != component_data.end(); it++) {
        if (TestTypeAndForm(*it, VTKOutputType::POINT_DATA, VTKDataAgglomerateType::SCALARS)) {
            datasets[POINT_AND_SCALAR].push_back(*it);
        } else if (TestTypeAndForm(*it, VTKOutputType::POINT_DATA, VTKDataAgglomerateType::VECTORS)) {
            datasets[POINT_AND_VECTOR].push_back(*it);
        } else if (TestTypeAndForm(*it, VTKOutputType::CELL_DATA, VTKDataAgglomerateType::SCALARS)) {
            datasets[CELL_AND_SCALAR].push_back(*it);
        } else if (TestTypeAndForm(*it, VTKOutputType::CELL_DATA, VTKDataAgglomerateType::VECTORS)) {
            datasets[CELL_AND_VECTOR].push_back(*it);
        } else {
            ERROR << "dataset with cannot be printed via VTK!" << ERROR_CLOSE;
        }
    }
    for (std::size_t n{0}; n < NUM_TYPES; n++) {
        if (datasets[n].empty())
            continue;

        os << in << in << "<P" << OutputStrings[n][0] << " Scalars=\"" << OutputStrings[n][1] << "\">\n";  // check that one
        for (auto it = datasets[n].begin(); it != datasets[n].end(); it++) {
            os << in << in << in << "<PDataArray type=\"Float64\" Name=\"" << it->name << "\" NumberOfComponents=\"" << std::to_string(it->number_components) << "\" format=\"appended\" offset=\"0\">\n";  // NOLINT
            os << in << in << in << "</PDataArray>\n";
        }
        os << in << in << "</P" << OutputStrings[n][0] << ">\n";
    }
}

void VTKPXMLStructuredGridWriter::AddPieces(std::ostringstream& os) {
    if (!is_root)
        return;
    if (static_cast<int>(path_pieces.size()) != exec_man->GetNumberProcesses()) {
        ERROR << "Number of provided paths is inconsistent with number of processes" << ERROR_CLOSE;
    }
    if (static_cast<int>(extent_pieces.size()) != (exec_man->GetNumberProcesses() * 6)) {
        ERROR << "Size of extents is inconsistent with number of processes" << ERROR_CLOSE;
    }

    for (int n{0}; n < static_cast<int>(path_pieces.size()); n++) {
        int offset = n * 6;
        std::string extent_str = std::to_string(extent_pieces[offset]);
        for (int n{1}; n < 6; n++) {
            extent_str += " ";
            extent_str += std::to_string(extent_pieces[offset + n]);
        }
        os << in << in << "<Piece Extent=\"" << extent_str << "\" Source=\"" << path_pieces[n] << "\"/>\n";
    }
}

void VTKPXMLStructuredGridWriter::AddTimeStamp(std::ostringstream& os) {
    if (!is_root || (time < 0.))
        return;
    os << in << in << "<FieldData>\n";
    os << in << in << in << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> " << std::setprecision(6) << time << '\n';  // NOLINT
    os << in << in << in << "</DataArray>\n";
    os << in << in << "</FieldData>\n";
}

void VTKPXMLStructuredGridWriter::AddClosing(std::ostringstream& os) {
    if (!is_root)
        return;
    os << in << "</PStructuredGrid>\n";
    os << "</VTKFile>\n";
}

}  // namespace dare::io
