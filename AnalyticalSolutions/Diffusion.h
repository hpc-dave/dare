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

#ifndef ANALYTICALSOLUTIONS_DIFFUSION_H_
#define ANALYTICALSOLUTIONS_DIFFUSION_H_

#include <math.h>
#include <cassert>
#include <cstdlib>

namespace dare::analytical {

/*!
 * @brief computes analytical solution for a transient diffusion problem in a slab
 * @tparam SC scalar type
 * @param tau dimensionless time
 * @param zeta dimensionless distance
 * @param n_modes number of modes to compute
 * @return analytical value
 * Here the problem is solved with following boundary and initial conditions:
 * phi(0, zeta) = 0
 * phi(tau > 0, 1) = 1
 * dphi/dzeta (tau > 0, 0) = 0
 */
template<typename SC>
SC Diffusion_N0D1(SC tau, SC zeta, int n_modes = 100) {
    // infinite solution for this problem is simply 1
    SC phi_inf{1.};
    SC phi_tr{0.};
    // now we need to add the transient component
    for (int n{0}; n < n_modes; ++n) {
        SC npi = (n + 0.5) * M_PI;
        phi_tr += 2 * std::cos(n * M_PI) / npi * std::exp(-npi * npi * tau) * std::cos(npi * zeta);
    }
    SC phi = phi_inf - phi_tr;
    return phi;
}

}  // namespace dare::analytical
#endif  // ANALYTICALSOLUTIONS_DIFFUSION_H_
