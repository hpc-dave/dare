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

#include "VTKWriter.h"

namespace dare::io::details {
std::string VTKGetParallelOutputPath(const dare::mpi::ExecutionManager& exman,
                                     const std::string& base_path,
                                     const std::string& pfolder_name) {
    const std::string proc_id = std::to_string(exman.GetRank());
    std::string parallel_data_path(base_path);
    if (!parallel_data_path.empty()) {
        if (base_path.back() != '/')
            parallel_data_path += '/';
    }
    parallel_data_path += pfolder_name;
    if (parallel_data_path.back() != '/')
        parallel_data_path += '/';
    return parallel_data_path;
}

std::string VTKGetOutputFileName(dare::mpi::ExecutionManager* exman,
                                         const std::string& parallel_data_path,
                                         const std::string& grid_name,
                                         int step,
                                         const std::string& ext) {
    std::ostringstream os;
    os << parallel_data_path;
    if (!parallel_data_path.empty() && (parallel_data_path.back() != '/'))
        os << '/';
    os << grid_name
       << '_' << std::to_string(step)
       << "." << ext;
    exman->Print(dare::mpi::Verbosity::Medium) << "Writing to file "
                                                  << parallel_data_path << "*" << grid_name
                                                  << '_' << std::to_string(step)
                                                  << "." << ext
                                                  << std::endl;
    return os.str();
}
std::string VTKGetParallelOutputFileName(dare::mpi::ExecutionManager* exman,
                                 const std::string& output_path,
                                 const std::string& grid_name,
                                 int step,
                                 const std::string& ext) {
    std::ostringstream os;
    os << output_path;
    if (!output_path.empty() && (output_path.back() != '/'))
        os << '/';
    os << grid_name
       << '_' << std::to_string(step)
       << "." << ext;
    exman->Print(dare::mpi::Verbosity::Medium) << "Writing to file "
                                               << output_path << grid_name
                                               << '_' << std::to_string(step)
                                               << "." << ext
                                               << std::endl;
    return os.str();
}

bool VTKWritePXMLRootFile(const std::string& path,
                          const VTKExtent& whole_extent,
                          int ghost_level,
                          const std::list<std::string>& path_distributed_files,
                          const std::list<VTKXMLPStructuredGridComponentData>& comp_data,
                          const std::list<VTKExtent>& extent_subdomains,
                          double time) {
    enum vtkTypeAndAgglomerates {
        POINT_AND_SCALAR,
        POINT_AND_VECTOR,
        CELL_AND_SCALAR,
        CELL_AND_VECTOR,
        NUM_TYPES
    };
    auto TestTypeAndForm = [](const VTKXMLPStructuredGridComponentData& c, std::string type, std::string form) {
        return boost::iequals(c.OutputType, type) && boost::iequals(c.DataAgglomerateType, form);
    };
    // prepare some strings
    std::string endianness = "LittleEndian";
   // substitute with the following for C++20
   // (std::endian::native == std::endian::little ? "LittleEndian" : "BigEndian");
    std::string whole_extent_str = std::to_string(whole_extent[0]);
    for (std::size_t n{1}; n < 6; n++) {
        whole_extent_str += " ";
        whole_extent_str += std::to_string(whole_extent[n]);
    }
    std::list<VTKXMLPStructuredGridComponentData> datasets[NUM_TYPES];
    for (auto it = comp_data.begin(); it != comp_data.end(); it++) {
        if (TestTypeAndForm(*it, "Point", "Scalar")) {
            datasets[POINT_AND_SCALAR].push_back(*it);
        } else if (TestTypeAndForm(*it, "Point", "Vector")) {
            datasets[POINT_AND_VECTOR].push_back(*it);
        } else if (TestTypeAndForm(*it, "Cell", "Scalar")) {
            datasets[CELL_AND_SCALAR].push_back(*it);
        } else if (TestTypeAndForm(*it, "Cell", "Vector")) {
            datasets[CELL_AND_VECTOR].push_back(*it);
        } else {
            ERROR << "dataset with " << it->OutputType << " and "
                  << it->DataAgglomerateType << " cannot be printed via VTK!" << ERROR_CLOSE;
        }
    }
    // Assemble file body
    std::ostringstream os;
    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"" << endianness << "\">\n";
    os << "<PStructuredGrid WholeExtent=\"" << whole_extent_str << "\" GhostLevel=\"" << std::to_string(ghost_level) << "\">\n";    // NOLINT
    os << "<PPoints>\n";
    os << "<PDataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\">\n";         // NOLINT
    os << "</PDataArray>\n";
    os << "</PPoints>\n";
    if (!datasets[POINT_AND_SCALAR].empty()) {
        os << "<PPointData Scalars=\"scalars\">\n";
        for (auto it = datasets[POINT_AND_SCALAR].begin(); it != datasets[POINT_AND_SCALAR].end(); it++) {
            os << "<PDataArray type=\"Float64\" Name=\"" << it->Name << "\" NumberOfComponents=\"" << it->NumberOfComponents << "\" format=\"appended\" offset=\"0\">\n";  // NOLINT
            os << "</PDataArray>\n";
        }
        os << "</PPointData>\n";
    }
    if (!datasets[POINT_AND_VECTOR].empty()) {
        os << "<PPointData Scalars=\"vectors\">\n";  // check that one
        for (auto it = datasets[POINT_AND_VECTOR].begin(); it != datasets[POINT_AND_VECTOR].end(); it++) {
            os << "<PDataArray type=\"Float64\" Name=\"" << it->Name << "\" NumberOfComponents=\"" << it->NumberOfComponents << "\" format=\"appended\" offset=\"0\">\n";  // NOLINT
            os << "</PDataArray>\n";
        }
        os << "</PPointData>\n";
    }
    if (!datasets[CELL_AND_SCALAR].empty()) {
        os << "<PCellData Scalars=\"scalars\">\n";
        for (auto it = datasets[CELL_AND_SCALAR].begin(); it != datasets[CELL_AND_SCALAR].end(); it++) {
            os << "<PDataArray type=\"Float64\" Name=\"" << it->Name << "\" NumberOfComponents=\"" << it->NumberOfComponents << "\" format=\"appended\" offset=\"0\">\n";  // NOLINT
            os << "</PDataArray>\n";
        }
        os << "</PCellData>\n";
    }
    if (!datasets[CELL_AND_VECTOR].empty()) {
        os << "<PCellData Scalars=\"vectors\">\n";
        for (auto it = datasets[CELL_AND_VECTOR].begin(); it != datasets[CELL_AND_VECTOR].end(); it++) {
            os << "<PDataArray type=\"Float64\" Name=\"" << it->Name << "\" NumberOfComponents=\"" << it->NumberOfComponents << "\" format=\"appended\" offset=\"0\">\n";  // NOLINT
            os << "</PDataArray>\n";
        }
        os << "</PCellData>\n";
    }
    if (path_distributed_files.size() != extent_subdomains.size()) {
        ERROR << "number of subdomains are inconsistent between paths and extents!" << ERROR_CLOSE;
        return false;
    }
    auto it_path = path_distributed_files.begin();
    auto it_extent = extent_subdomains.begin();
    for (std::size_t n{0}; n < path_distributed_files.size(); n++) {
        std::string extent_str = std::to_string((*it_extent)[0]);
        for (std::size_t n{1}; n < 6; n++) {
            extent_str += " ";
            extent_str += std::to_string((*it_extent)[n]);
        }
        os << "<Piece Extent=\"" << extent_str << "\" Source=\"" << *it_path << "\"/>\n";
        ++it_path;
        ++it_extent;
    }
    if (time >= 0.) {
        os << "<FieldData>\n";
        os << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> " << std::setprecision(6) << time << '\n'; //NOLINT
        os << "</DataArray>\n";
        os << "</FieldData>\n";
        os << "</PStructuredGrid>\n";
        os << "</VTKFile>\n";
    }
    std::ofstream ofs(path, std::ios::out);
    if (!ofs) {
        ERROR << "Could not open file: " << path << ERROR_CLOSE;
        return false;
    }
    ofs << os.str();
    return true;
}
}  // end namespace dare::io::details
