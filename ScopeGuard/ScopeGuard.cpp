/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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

#include "ScopeGuard.h"

namespace dare {

ScopeGuard::ScopeGuard(int* _argc, char*** _argv, bool suppress_output) : argc(_argc), argv(_argv) {
    int num_proc = 1;
    int my_rank = 0;

    int flag_mpi_init = 0;
    MPI_Initialized(&flag_mpi_init);

    manage_mpi = !static_cast<bool>(flag_mpi_init);

    if (manage_mpi)
        MPI_Init(_argc, _argv);

    tpetra_scope = Teuchos::rcp(new Tpetra::ScopeGuard(_argc, _argv));

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    is_root = my_rank == 0;
    if (AmIRoot() && !suppress_output)
        std::cout << "\nRunning with " << num_proc << " procs.\n"
                  << std::endl;

    std::string option;
    if (HasArgument("-D", &option) || HasArgument("--debug", &option)) {
        bool skip_by_error{false};
        int proc{-1};
        try {
            proc = std::stoi(option);
        } catch (std::exception& ex) {
            if (AmIRoot())
                std::cout << "Cannot interpret argument of -D/--debug, wrong format! Found value: " << option;
            skip_by_error = true;
        }
        if (!skip_by_error && (proc == my_rank)) {
            int stop{0};
            while (stop == 0) {
                // if you end up here, then change the value of stop in your debugger
                // and continue the debugging
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (HasArgument("-H") || HasArgument("--help")) {
        if (AmIRoot())
            PrintHelp();
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(0);
    }
}

ScopeGuard::~ScopeGuard() {
    tpetra_scope = Teuchos::ENull::null;
    if (manage_mpi)
        MPI_Finalize();
}

bool ScopeGuard::HasArgument(std::string arg, std::string* options) const {
    for (int n = 1; n < *argc; ++n) {
        if (boost::iequals(arg, (*argv)[n])) {
            if (options != nullptr && n < (*argc - 1))
                *options = (*argv)[n + 1];
            return true;
        }
    }
    return false;
}

bool ScopeGuard::AmIRoot() const {
    return is_root;
}

void ScopeGuard::Terminate(std::string message, int error_code) const {
    if (is_root)
        std::cerr << "Terminating the execution! " << message << "!" << std::endl;

    MPI_Abort(MPI_COMM_WORLD, error_code);
}

void ScopeGuard::PrintHelp() {
    std::cout << "Use following command line options with this executable:\n";
    std::cout << "    -D    or --debug <rank>          The process of <rank> will be caught in an infinite loop,\n"
              << "                                     which can be escaped from with a debugger.\n";
    std::cout << "    -H    or --help                  Displays this message\n";
    std::cout << std::endl;
}
}  // end namespace dare
