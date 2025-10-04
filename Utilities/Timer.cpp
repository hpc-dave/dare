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

#include "Utilities/Errors.h"
#include "Timer.h"

namespace dare {

Timer::Timer(bool blocking, MPI_Comm comm) : is_blocking{blocking}, communicator{comm} {}
Timer::~Timer() {}

void Timer::Tic(std::string ident) {
    Synchronize();
    _data[ident].time_begin = MPI_Wtime();
}

Timer::ValueType Timer::Toc(std::string ident) {
    auto it = _data.find(ident);
#ifndef DARE_NDEBUG
    if (it == _data.end()) {
        ERROR << "Could not find " << ident << " in the timer instance" << ERROR_CLOSE;
        return 0.;
    }
#endif
    Synchronize();
    it->second.time_elapsed = MPI_Wtime() - it->second.time_begin;
    it->second.time_elapsed_total += it->second.time_elapsed;
    return it->second.time_elapsed;
}

Timer::ValueType Timer::GetElapsedTime(std::string ident) const {
    auto it = _data.find(ident);
#ifndef DARE_NDEBUG
    if (it == _data.end()) {
        ERROR << "Could not find " << ident << " in the timer instance" << ERROR_CLOSE;
        return 0.;
    }
#endif
    return it->second.time_elapsed;
}

Timer::ValueType Timer::GetTotalElapsedTime(std::string ident) const {
    auto it = _data.find(ident);
#ifndef DARE_NDEBUG
    if (it == _data.end()) {
        ERROR << "Could not find " << ident << " in the timer instance" << ERROR_CLOSE;
        return 0.;
    }
#endif
    return it->second.time_elapsed_total;
}

void Timer::Synchronize() {
    if (is_blocking)
        MPI_Barrier(communicator);
}

}  // namespace dare
