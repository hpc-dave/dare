/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder

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

#include <mpi.h>
#include <string>
#include <iostream>

#include "IO/TerminalOutput.h"
#include "Utilities/Errors.h"
#include "PMTerminalOutputDefault.h"

namespace dare {

PMTerminalOutputDefault::PMTerminalOutputDefault(int dim, int acc, int width)
: dimension{dim}, accuracy{acc}, print_width{width} {
}

PMTerminalOutputDefault::~PMTerminalOutputDefault() {}

void PMTerminalOutputDefault::PrintHeader(uint64_t tstep, double time) {
    dare::Print(dare::Verbosity::Low) << "step: " << std::to_string(tstep) << " - time: " << time << std::endl;
    dare::Print(dare::Verbosity::Low) << "| ";
}

void PMTerminalOutputDefault::PrintMomentum(int dim, int iter, bool success) {
    if (dim >= dimension) {
        ERROR << "provided dimension is larger than specified (" << dim << " >= " << dimension << ")" << ERROR_CLOSE;
    }
    std::string f_id;
    switch (dim) {
    case 0:
        f_id = "X";
        break;
    case 1:
        f_id = "Y";
        break;
    case 2:
        f_id = "Z";
        break;
    default:
        ERROR << "Invalid dimension for output: " << std::to_string(dim) << ERROR_CLOSE;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    Print(dare::Verbosity::Low) << f_id << ": " << std::to_string(iter) << " it "
            << (success? "": "FAILURE ") << "| " << std::flush;
}

void PMTerminalOutputDefault::PrintInitialDefectInternal(double defect) {
    Print(dare::Verbosity::Low) << "defect: " << defect << '\n'
                                << "Continuity:\n"
                                << std::flush;
}


void PMTerminalOutputDefault::PrintContinuity(int loop, int iter, bool success, double defect, bool is_last) {
    Print(dare::Verbosity::Low)
        << std::to_string(loop) << " -> "
        << (is_last? "MAX" : std::to_string(iter))
        << " it - defect: " << defect << (success? "" : "FAILURE") <<  '\n' << std::flush;
}

}  // namespace dare
