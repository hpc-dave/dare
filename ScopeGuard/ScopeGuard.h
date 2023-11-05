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

#ifndef SCOPEGUARD_SCOPEGUARD_H_
#define SCOPEGUARD_SCOPEGUARD_H_

#include <mpi.h>
#include <omp.h>

#include <iostream>
#include <string>
#include <boost/algorithm/string/predicate.hpp>


namespace dare {

/*! \class ScopeGuard
 * \brief helper to take care of environment settings
 *
 * This helper takes care of settings in OpenMP and MPI environments.
 * The ScopeGuard should be the <b> first </b> object to be constructed,
 * so it will also be the <b> last </b> object to go out of scope.
 *
 * If the MPI_Init is called prior to the ScopeGuard, you are responsible to
 * call MPI_Finalize yourself!
 *
 * Application should be as follows:
 * int main(int argc, char* argv[]){
 *
 *   FoxBerry::ScopeGuard scope(argc, argv);
 *
 *   // execute your code here
 *
 * }
 */

class ScopeGuard {
public:
    /*!
     * \brief constructor
     * @param argc number of arguments
     * @param argv pointer to the array with string arguments
     * @param suppress_output identifier, if anything should be printed
     *
     * The constructor will initialize the environments according to the provided
     * terminal arguments
     */
    ScopeGuard(int argc, char* argv[], bool suppress_output = false);

    /*!
     * \brief default destructor
     * Any environment which was activated in the constructor is close here
     */
    virtual ~ScopeGuard();

    /*!
     * \brief returns true, if the terminal constains a certain arguments
     * @param arg argument of test for
     * @param options will be filled with the following argument, if it exists
     *
     * Example:
     *
     * In terminal: ./exec -c 1
     *
     * In code:
     * int main(int argc, char** argv){
     *
     * FoxBerry::ScopeGuard scope(argc, argv);
     * {
     *   if(scope.HasArgument("-c")){
     *    // will execute, since "-c" exists in the list or arguments
     *   }
     *
     *   if(scope.HasArgument("c"){
     *    // won't execute, since "c" is only part of an argument, not a full argument!
     *   }
     * }
     */
    bool HasArgument(std::string arg, std::string* options = nullptr) const;

    /*!
     * \brief true, if this is the root process
     */
    bool AmIRoot() const;

    /*!
     * \brief determines the maximum number of threads for the process
     * \note only relevant if running with OpenMP
     */
    int GetMaximumThreadsForProcess() const;

    /*!
     * \brief terminates the execution and returns an error code
     */
    void Terminate(std::string message, int error_code = 1) const;

private:
    /*!
     * @brief prints help message
     */
    void PrintHelp();

    int argc;         // number of arguments
    char** argv;      // argument array

    bool is_root;     // identifier, if this process is root
    bool manage_mpi;  // identifier, if this instance is responsible for managing MPI
};

}  // end namespace dare

#endif  // SCOPEGUARD_SCOPEGUARD_H_
