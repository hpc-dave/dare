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

namespace dare::algorithm {

template <typename PM>
void free_compile_time_check(PM) {
    // by default this is empty and accepts everything, can be overloaded for certain properties
}

template <typename PM, typename Grid, typename... Args>
void free_pm_initialize(PM* pm, const Grid& grid, Args&&... args) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method and grid");   // NOLINT
}

template <typename PM, std::size_t dir>
void free_pm_build_momentum(PM* pm) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method");  // NOLINT
    // using GridType = typename PM::GridType;
    // using LO = typename GridType::LocalOrdinalType;
    // using DensityType = typename PM::DensityVariableType;
    // using ViscosityType = typename PM::ViscosityVariableType;
    // using PorosityType = typename PM::ViscosityVariableType;
    // using IndexLocal = typename GridType::Index;
    // using DDT = dare::Matrix::DDT<GridType>;

    // auto BuildStrategy = [= pm](auto mblock) {
    //     auto g_r{mblock->GetRepresentation()};
    //     LO o_loc{mblock->GetLocalOrdinal()};
    //     IndexLocal ind{mblock->GetIndex()};

    //     DDT ddt(*g_r, o_loc, pm->GetTimeStepSize());

    //     DensityType rho{pm->GetDensity()};
    //     ViscosityType mu{pm->GetViscosity()};
    //     PorosityType epsilon{pm->GetPorosity()};
    //     // accumulation
    //     (*mblock) = ddt(epsilon, rho, *pm->GetMomentum(dir));
    // };

    // momentum[dir]->Build(BuildStrategy);
}

template <typename PM>
std::pair<bool, int> free_pm_solve_momentum(PM* pm, std::size_t dir) {
    using IterType = typename PM::MomentumIterationType;
    std::pair<bool, int> = std::make_pair(false, static_cast<int>(-1));
    if constexpr (uses_fixed_point_iterations_v<IterType>) {
        ret = continuity.Solve(dare::Matrix::UpdateFieldCopy{});
    } else if constexpr (uses_newton_iterations_v<IterType>) {
        ret = continuity.Solve(dare::Matrix::UpdateFieldAddInto{});
    } else {
        static_assert(dare::always_false<IterType>, "Solving the momentum equations is not implemented for the specified algorithm type");  // NOLINT
    }
    return ret;
}

template <typename PM>
void free_pm_build_continuity(PM* pm, int iteration) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method and grid");  // NOLINT
}

template <typename PM>
std::pair<bool, int> free_pm_solve_continuity(PM* pm, int iteration) {
    using IterType = typename PM::ContinuityIterationType;
    std::pair<bool, int> = std::make_pair(false, static_cast<int>(-1));
    if constexpr (uses_newton_iterations_v<IterType>) {
        auto UpdateStrategy = [](const typename Continuity::MatrixSystemType& m, typename PM::ContinuityType* c) {
            m.CopyTo(&c->GetdP()->GetDataVector());
        };
        ret = continuity.Solve(UpdateStrategy);
    } else {
        static_assert(dare::always_false<IterType>, "Updating the pressure is not implemented for anything except Newton iterations");  // NOLINT
    }
    return ret;
}

template <typename PM>
void free_pm_update_pressure(PM* pm, int iteration) {
    if constexpr (uses_newton_iterations_v<typename PM::ContinuityIterationType>) {
        auto UpdateStrategy = [](const typename Continuity::MatrixSystemType& m, typename PM::ContinuityType* c) {
            m.AddTo(&c->GetField()->GetDataVector());
        };
        continuity.Solve(UpdateStrategy);
    } else {
        static_assert(dare::always_false<PM>, "Updating the pressure is not implemented for anything except Newton iterations");  // NOLINT
    }
}

template <typename PM>
void free_pm_update_velocity(PM* pm, int iteration) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method and grid");  // NOLINT
}

template <typename PM>
bool free_pm_continuity_convergence(PM* pm, int iteration) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method and grid");  // NOLINT
    return false;
}

}  // namespace dare::algorithm
