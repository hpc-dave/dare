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

#include <gtest/gtest.h>
#include <mpi.h>

#include <boost/algorithm/string/predicate.hpp>

#include "ScopeGuard/ScopeGuard.h"

static std::string xml_report_arg = "";  // NOLINT

int main(int argc, char* argv[]) {
    // initialize MPI environment
    dare::ScopeGuard scope_guard(&argc, &argv);
    {
        int rank{-1};
        int num_proc{-1};
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

        // if output should be generated, we adapt the output files name to include the process ID
        for (int n{0}; n < argc; n++) {
            std::string arg = argv[n];
            if (arg.size() > 15) {
                if (boost::iequals(arg.substr(0, 15), "--gtest_output=")) {
                    xml_report_arg = arg;
                    std::string proc_id = "_";
                    proc_id += std::to_string(rank);
                    xml_report_arg.insert(xml_report_arg.size() - 4, proc_id, 0, proc_id.size());
                    argv[n] = xml_report_arg.data();
                }
            }
        }

        int ret = -1;
        {
            // initializing google test
            testing::InitGoogleTest(&argc, argv);
            // delete the default printers to avoid confusing output on COUT
            testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
            if (!scope_guard.AmIRoot()) {
                // remove all output to cout from non-root processes
                delete listeners.Release(listeners.default_result_printer());
                // remove report writers
                // delete listeners.Release(listeners.default_xml_generator());
            }
            // run the tests
            ret = RUN_ALL_TESTS();
        }

        bool result = ret == MPI_SUCCESS;
        int did_fail = ret == 0 ? 0 : 1;
        bool result_global{false};
        int n_fail{0};
        MPI_Allreduce(&result, &result_global, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
        MPI_Allreduce(&did_fail, &n_fail, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (result_global && scope_guard.AmIRoot()) {
            std::cout << "All tests were successful over all processes" << std::endl;
        } else if (scope_guard.AmIRoot()) {
            std::cout << "On " << n_fail << " out of " << num_proc << " process" << (num_proc > 1 ? "es" : "")
                      << " the tests failed!" << std::endl;
        }

        return ret;
    }
}
