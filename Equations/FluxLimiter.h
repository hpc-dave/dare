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

#ifndef EQUATIONS_FLUXLIMITER_H_
#define EQUATIONS_FLUXLIMITER_H_

#include <algorithm>
#include <cmath>

#include "Utilities/Vector.h"

namespace dare::Matrix {

/*!
 * \brief central differencing limiter for TVD schemes
 * \f[
 *  \Psi(r) = 1
 * \f]
 */
struct CDS {
    /*!
     * @brief computes flux-limiter for multiple values
     * @tparam SC type of scalars
     * @tparam N number of components
     * @param r face-gradient
     */
    template <std::size_t N, typename SC>
    static dare::utils::Vector<N, SC> GetValue(const dare::utils::Vector<N, SC>& r) {
        dare::utils::Vector<N, SC> ret;
        for (std::size_t n{0}; n < N; n++)
            ret[n] = GetValue(r[n]);
        return ret;
    }

    /*!
     * @brief computes flux-limiter for single value
     * @tparam SC type of scalar
     * @param r face-gradient
     */
    template <typename SC>
    static SC GetValue(const SC& r) {
        return 1.;
    }
};

/*!
 * \brief upwind limiter for TVD schemes
 * \f[
 *  \Psi(r) = 0
 * \f]
 */
struct UPWIND {
    /*!
     * @brief computes flux-limiter for multiple values
     * @tparam SC type of scalars
     * @tparam N number of components
     * @param r face-gradient
     */
    template<std::size_t N, typename SC>
    static dare::utils::Vector<N, SC> GetValue(const dare::utils::Vector<N, SC>& r) {
        dare::utils::Vector<N, SC> ret;
        for (std::size_t n{0}; n < N; n++)
            ret[n] = GetValue(r[n]);
        return ret;
    }

    /*!
     * @brief computes flux-limiter for single value
     * @tparam SC type of scalar
     * @param r face-gradient
     */
    template <typename SC>
    static SC GetValue(const SC& r) {
        return 0.;
    }
};

/*!
 * \brief Van Albada's limiter for TVD schemes
 * \f[
 *  \Psi(r) = \frac{r+r^2}{1+r^2}
 * \f]
 */
struct VANALBADA {
    /*!
     * @brief computes flux-limiter for multiple values
     * @tparam SC type of scalars
     * @tparam N number of components
     * @param r face-gradient
     */
    template <std::size_t N, typename SC>
    static dare::utils::Vector<N, SC> GetValue(const dare::utils::Vector<N, SC>& r) {
        dare::utils::Vector<N, SC> ret;
        for (std::size_t n{0}; n < N; n++)
            ret[n] = GetValue(r[n]);
        return ret;
    }

    /*!
     * @brief computes flux-limiter for single value
     * @tparam SC type of scalar
     * @param r face-gradient
     * In the case of r == NaN, the value 0 will be returned
     */
    template <typename SC>
    static SC GetValue(SC r) {
        const SC ZERO{0};
        r = std::isnan(r) ? ZERO : r;
        const SC r_sqr{r * r};
        return (r + r_sqr) / (1. + r_sqr);
    }
};

/*!
 * \brief MinMod limiter for TVD schemes
 * \f[
 *  \Psi(r) = \begin{cases}
 *              0 &\quad r\leq 0 \\
 *              min(r,1) & \quad r > 0
 *            \end{cases}
 * \f]
 * \note sometimes referred to as second order Barton scheme
 */
struct MINMOD {
    /*!
     * @brief computes flux-limiter for multiple values
     * @tparam SC type of scalars
     * @tparam N number of components
     * @param r face-gradient
     */
    template <std::size_t N, typename SC>
    static dare::utils::Vector<N, SC> GetValue(const dare::utils::Vector<N, SC>& r) {
        dare::utils::Vector<N, SC> ret;
        for (std::size_t n{0}; n < N; n++)
            ret[n] = GetValue(r[n]);
        return ret;
    }

    /*!
     * @brief computes flux-limiter for single value
     * @tparam SC type of scalar
     * @param r face-gradient
     * If the face-gradient is NaN, the limit value is 0
     */
    template <typename SC>
    static SC GetValue(SC r) {
        const SC ONE{1};
        const SC ZERO{0};
        r = std::isnan(r) ? ZERO : r;
        return std::max(ZERO, std::min(r, ONE));
    }
};
}  // namespace dare::Matrix

#endif  // EQUATIONS_FLUXLIMITER_H_
