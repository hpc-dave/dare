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

#ifndef ALGORITHM_NAVIERSTOKES_PM_MOMENTUM_H_
#define ALGORITHM_NAVIERSTOKES_PM_MOMENTUM_H_

#include <string>
#include <set>
#include <functional>
#include <utility>

#include "Data/Field.h"
#include "MPI/ExecutionManager.h"
#include "MatrixSystem/Trilinos.h"
#include "MatrixSystem/TrilinosSolver.h"

namespace dare::algorithm {

// consider upgrading this class to a general momentum class
template<typename Grid, typename BoundaryStrategy, typename CustomMember>
class PMMomentum {
public:
    using GridType = Grid;
    using GridRepresentation = typename GridType::Representation;
    using SC = typename GridType::ScalarType;
    using FieldType = dare::Data::Field<GridType, SC, 1>;
    using BoundaryStrategyType = BoundaryStrategy;
    using CustomMemberType = CustomMember;
    using MatrixSystemType = dare::Matrix::Trilinos<SC>;
    using MatrixSolverType = dare::Matrix::TrilinosSolver<SC>;
    using PreconditionerPropertyType = typename MatrixSolverType::PropertyType;
    using SolverPropertyType = typename MatrixSolverType::PropertyType;
    using SelfType = PMMomentum<Grid, BoundaryStrategy, CustomMember>;

    PMMomentum(const std::string& name,
               GridRepresentation grid,
               dare::mpi::ExecutionManager* ex_man,
               std::size_t num_tsteps,
               BoundaryStrategy bc_strat);

    explicit PMMomentum(const SelfType&) = delete;
    SelfType& operator=(const SelfType&) = delete;

    void PreStep();

    template <typename BuildStrategy>
    void Build(BuildStrategy build);

    std::pair<bool, int> Solve(SolverPropertyType sprop, PreconditionerPropertyType mprop);

    void UpdateBoundaries();

    void PostStep();

    BoundaryStrategyType* GetBoundaryStrategy();
    const BoundaryStrategyType& GetBoundaryStrategy() const;

    FieldType* GetField();
    const FieldType& GetField() const;

    CustomMemberType* GetCustomMember();
    const CustomMemberType& GetCustomMember() const;

    void AddPreStepStrategy(std::function<void(SelfType*)> f);
    void SetPreStepStrategy(std::function<void(SelfType*)> f);
    void ClearPreStepStrategy();
    void AddPostStepStrategy(std::function<void(SelfType*)> f);
    void SetPostStepStrategy(std::function<void(SelfType*)> f);
    void ClearPostStepStrategy();

private:
    GridRepresentation grep;
    dare::mpi::ExecutionManager* exec_man;
    FieldType data;
    BoundaryStrategyType boundary_strategy;
    CustomMemberType custom_member;
    std::set<std::function<void(SelfType*)>> pre_step_strategy;
    std::set<std::function<void(SelfType*)>> post_step_strategy;
    MatrixSystemType matrix_system;
};

}  // namespace dare::algorithm

#include "PM_Momentum.inl"

#endif  // ALGORITHM_NAVIERSTOKES_PM_MOMENTUM_H_
