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

#include <vtkIndent.h>

#include <fstream>
#include <iomanip>
#include <type_traits>

namespace dare::io {

vtkStandardNewMacro(VTKPXMLStructuredGridWriter);

void VTKPXMLStructuredGridWriter::SetPPieceExtent(const VTKExtent& local_extent, dare::mpi::ExecutionManager* exman) {
    extent_array.resize(exman->GetNumberProcesses() * 6);
    if ((exman->Allgather(local_extent.data(), 6, extent_array.data(), 6) != MPI_SUCCESS) && exman->AmIRoot())
        ERROR << "A not further specified problem occured during communication!" << ERROR_CLOSE;
}

void VTKPXMLStructuredGridWriter::WritePPieceAttributes(int index) {
    int offset = index * 6;
    std::string extent_s = std::to_string(extent_array[offset]);
    for (int n{1}; n < 6; n++) {
        extent_s += " ";
        extent_s += std::to_string(extent_array[offset + n]);
    }
    std::string piecename(this->CreatePieceFileName(index));
    this->WriteStringAttribute("Source", piecename.c_str());
    this->WriteStringAttribute("Extent", extent_s.c_str());
}

}  // namespace dare::io
