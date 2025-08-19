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

#ifndef UTILITIES_TIMER_H_
#define UTILITIES_TIMER_H_

#include <mpi.h>
#include <string>
#include <iostream>
#include <map>

namespace dare {
class Timer {
public:
    using ValueType = double;
    struct ValueTriplet {
        explicit ValueTriplet(ValueType tic = 0., ValueType t_el = 0., ValueType t_tot = 0.)
            : time_begin{tic}, time_elapsed{t_el}, time_elapsed_total{t_tot} {}
        ValueType time_begin;
        ValueType time_elapsed;
        ValueType time_elapsed_total;
    };
    using Container = std::map<std::string, ValueTriplet>;

    explicit Timer(bool is_blocking = true, MPI_Comm comm = MPI_COMM_WORLD);
    virtual ~Timer();

    void Tic(std::string ident);
    ValueType Toc(std::string ident);

    ValueType GetElapsedTime(std::string ident) const;
    ValueType GetTotalElapsedTime(std::string ident) const;

    template <typename OS>
    friend inline OS& operator<<(OS& os, const Timer& timer);

private:
    void Synchronize();

    bool is_blocking;
    MPI_Comm communicator;
    Container _data;
};

template <typename OS>
OS& operator<<(OS& os, const dare::Timer& timer) {
    for (const auto& e : timer._data) {
        os << e.first << ": " << e.second.time_elapsed << " s - total: "
           << e.second.time_elapsed_total << " s" << std::endl;
    }
}

}  // namespace dare


#endif  // UTILITIES_TIMER_H_
