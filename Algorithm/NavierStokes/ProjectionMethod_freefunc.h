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

#ifndef ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_FREEFUNC_H_
#define ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_FREEFUNC_H_
#include <utility>

#include "Utilities/Errors.h"

namespace dare {

template <typename PM>
void free_compile_time_check(PM*) {
    // by default this is empty and accepts everything, can be overloaded for certain properties
}

template <typename PM, typename Grid, typename... Args>
void free_pm_initialize(PM* pm, Grid* grid, Args&&... args) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the initialization using the specified types of the projection method and grid");   // NOLINT
}

template <typename PM, dare::NaturalNumber Direction>
void free_pm_build_momentum(PM* pm, Direction) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method");  // NOLINT
}

template <typename PM, typename Normalizer, typename MatrixBlock>
void free_pm_apply_normalizer(PM* pm, Normalizer normalizer, MatrixBlock* mb) {
    if constexpr (!dare::is_none_v<Normalizer>)
        static_assert(dare::always_false<Normalizer>, "Unknown normalizer treatment");
}

template<typename PM, typename Normalizer, typename MatrixBlock>
    requires std::is_arithmetic_v<Normalizer>
void free_pm_apply_normalizer(PM* pm, Normalizer normalizer, MatrixBlock* mb) {
    (*mb) *= normalizer;
}

template <typename PM, typename Normalizer, typename DataVector>
void free_pm_revert_normalizer(PM* pm, Normalizer normalizer, DataVector* vec) {
    if constexpr (!dare::is_none_v<Normalizer>)
        static_assert(dare::always_false<Normalizer>, "Unknown normalizer treatment");
}

template <typename PM, typename Normalizer, typename DataVector>
    requires std::is_arithmetic_v<Normalizer>
void free_pm_revert_normalizer(PM* pm, Normalizer normalizer, DataVector* vec) {
    Normalizer r_norm{1. / normalizer};
    for (std::size_t n{0}; n < vec->GetSize(); n++)
        (*vec)[n] *= r_norm;
}

template <typename PM, typename MatrixBlock>
void free_pm_momentum_initialguess(PM* pm,
                                   typename PM::SC center_coef,
                                   MatrixBlock* mb) {
    mb->GetInitialGuess(0) += mb->GetRhs(0) / center_coef;
}

template <typename PM, typename MatrixBlock>
void free_pm_momentum_apply_boundary_conditions(PM* pm,
                                                const typename PM::MomentumType& momentum,
                                                MatrixBlock* mb) {
    momentum.GetBoundaryStrategy().Apply(mb);
}

template <typename PM, dare::NaturalNumber Direction>
std::pair<bool, int> free_pm_solve_momentum(PM* pm, Direction dir) {
    using IterType = typename PM::MomentumIterationType;
    std::pair<bool, int> ret = std::make_pair(false, static_cast<int>(-1));
    if constexpr (uses_fixed_point_iterations_v<IterType>) {
        ret = pm->GetMomentum(dir)->Solve(dare::UpdateFieldCopy{});
    } else {
        static_assert(dare::always_false<IterType>, "Solving the momentum equations is not implemented for the specified algorithm type");  // NOLINT
    }

    // Reverse normalization
    free_pm_revert_normalizer(pm,
                              pm->GetMomentum(dir)->GetCustomMember()->normalizer,
                              pm->GetMomentum()->GetField()->GetDataVector());
    return ret;
}

template <typename PM>
void free_pm_build_continuity(PM* pm, int iteration) {
    static_assert(dare::always_false<PM>, "Could not find the specialization for the specified types of the projection method and grid");  // NOLINT
}

template <typename PM>
std::pair<bool, int> free_pm_solve_continuity(PM* pm, int iteration) {
    using IterType = typename PM::ContinuityIterationType;
    using ContType = typename PM::ContinuityType;
    using ContMSystem = typename ContType::MatrixSystemType;
    std::pair<bool, int> ret = std::make_pair(false, static_cast<int>(-1));
    if constexpr (uses_newton_iterations_v<IterType>) {
        auto UpdateStrategy = [=](const ContMSystem& m, ContType* c) {
            m.CopyTo(&pm->GetContinuity()->GetdP()->GetDataVector());
        };
        ret = pm->GetContinuity()->Solve(UpdateStrategy);
    } else {
        static_assert(dare::always_false<IterType>, "Updating the pressure is not implemented for anything except Newton iterations");  // NOLINT
    }
    return ret;
}

template <typename PM>
void free_pm_update_pressure(PM* pm, int iteration) {
    if constexpr (uses_newton_iterations_v<typename PM::ContinuityIterationType>) {
        *(pm->GetContinuity()->GetPressure()) += *(pm->GetContinuity()->GetdP());
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

}  // namespace dare

#endif  // ALGORITHM_NAVIERSTOKES_PROJECTIONMETHOD_FREEFUNC_H_
