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
#include "Grid/Cartesian.h"

namespace dare::algorithm {

template <typename PM, std::size_t Dim, typename... Args>
void free_pm_initialize(PM* pm, const dare::Grid::Cartesian<Dim>& grid, Args&&... bc_args) {
    static_assert(PM::dimension == Dim, "The projection method and grid do not have the same dimension!");  // NOLINT
    static_assert(PM::Dimension < 3, "Not equipped for higher dimensions");
    std::string m_names[] = {"u", "v", "w"};
    typename PM::GridType::Options opt;
    for (auto& o : opt)
        o = 0.;

    for (std::size_t d{0}; d < Dim; d++) {
        auto opt_loc = opt;
        opt_loc[d] = 1;
        pm->InitializeMomentum(d,
                               m_names[d],
                               grid.GetRepresentation(opt_loc),
                               grid.GetExecutionManager(),
                               PM::num_tsteps_momentum,
                               bc_args...);
    }
    pm->GetContinuity()->Initialize("pressure",
                                    grid,
                                    grid.GetRepresentation(opt),
                                    grid.GetExecutionManager(),
                                    2,
                                    bc_args...);
}

template <typename PM>
    requires(std::is_same_v<typename PM::GridType, dare::Grid::Cartesian<PM::dimension>>)
void free_pm_solve(PM* pm) {
    static_assert(dare::always_false<PM>, "not yet implemented");
}

}  // namespace dare::algorithm
