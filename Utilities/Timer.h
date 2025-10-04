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
/*! \class Timer
 * @brief allows timing by means of an indentifier
 */
class Timer {
public:
    using ValueType = double;       //!< Time value

    /*!
     * @brief a helper structure to store the time data
     */
    struct ValueTriplet {
        explicit ValueTriplet(ValueType tic = 0., ValueType t_el = 0., ValueType t_tot = 0.)
            : time_begin{tic}, time_elapsed{t_el}, time_elapsed_total{t_tot} {}
        ValueType time_begin;           //!< time at which profiling starts (Tic)
        ValueType time_elapsed;         //!< elapsed time when profiling stops (Toc)
        ValueType time_elapsed_total;   //!< summed up time of all tic-toc pairs
    };
    using Container = std::map<std::string, ValueTriplet>;  //!< container for storing all members

    /*!
     * @brief constructor with some customization options
     * @param is_blocking flag, if the timer should synchronize all processes before timing
     * @param comm communicator to block the processes on
     */
    explicit Timer(bool is_blocking = true, MPI_Comm comm = MPI_COMM_WORLD);

    /*!
     * @brief default destructor
     */
    virtual ~Timer();

    /*!
     * @brief starts profiling
     * @param ident unique identifier
     * \note if the time is set as blocking, an MPI-Barrier is called
     */
    void Tic(std::string ident);

    /*!
     * @brief determines elapsed time since last Tic
     * @param ident unique identifier
     * @return elapsed time since last Tic in s
     * \note if the time is set as blocking, an MPI-Barrier is called
     * if no Tic on the identifier was called before, a 0 is returned
     */
    ValueType Toc(std::string ident);

    /*!
     * @brief provides the latest elapsed time for the identifier
     * @param ident identifier
     * @return elapsed time for the identifier in s
     */
    ValueType GetElapsedTime(std::string ident) const;

    /*!
     * @brief provides the sum of all elapsed times for an indentifier
     * @param ident identifier
     * @return elapsed time in s
     */
    ValueType GetTotalElapsedTime(std::string ident) const;

    /*!
     * @brief convenient printing option
     * @tparam OS ostream type of object
     * @param os ostream
     * @param timer the timer instance
     * @return the ostream after addition of the times
     */
    template <typename OS>
    friend inline OS& operator<<(OS& os, const Timer& timer);

private:
    /*!
     * @brief if the timer is set to blocking, all processes will be synchronized at this point
     */
    void Synchronize();

    bool is_blocking;       //!< flag if the timer is blocking
    MPI_Comm communicator;  //!< communicator on which to block
    Container _data;        //!< the actual data
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
