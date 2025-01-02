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

#include <string>

namespace dare::io::details {

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

}  // end namespace dare::io::details
