/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder
 *
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

#ifndef IO_TERMINALOUTPUT_H_
#define IO_TERMINALOUTPUT_H_

#include "BlackHoleOStream.h"

namespace dare {
enum class Verbosity : uint8_t {
    None,    //!< All output is swallowed
    Low,     //!< minimum production output
    Medium,  //!< a bit more information for rough debugging
    High     //!< huge amount of output, only for detailed debuggin
};

namespace detail {
static Verbosity verbosity_level{Verbosity::Low};
static BlackHoleOStream black_hole_ostream;
static bool is_root_for_print{true};
inline void SetRootForPrint(bool v) { is_root_for_print = v; }
inline bool IsRootForPrint() { return is_root_for_print; }
}  // namespace detail

inline void SetVerbosity(Verbosity level) {
    detail::verbosity_level = level;
}

/*!
 * \brief allows printing with output control
 * @param level level below which the output will be swallowed
 * Messages provided via this function will only be printed
 * by the root processor and if the \p level is lower than
 * the internally specified output verbosity
 * \note not threadsafe, as all output to terminal
 * \note must not be called before allocating ScopeGuard
 */
inline std::ostream& Print(Verbosity level) {
    if (level > detail::verbosity_level || !detail::IsRootForPrint())
        return detail::black_hole_ostream;
    else
        return std::cout;
}

/*!
 * \brief allows printing with output control
 * @param level level below which the output will be swallowed
 * Messages provided via this function will only be printed
 * by the root processor and if the \p level is lower than
 * the internally specified output verbosity
 * \note not threadsafe, as all output to terminal
 * \note must not be called before allocating ScopeGuard
 */
inline std::ostream& PrintAll(Verbosity level) {
    if (level > detail::verbosity_level)
        return detail::black_hole_ostream;
    else
        return std::cout;
}

}  // namespace dare

#endif  // IO_TERMINALOUTPUT_H_
