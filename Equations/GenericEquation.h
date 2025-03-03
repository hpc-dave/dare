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

#ifndef EQUATIONS_GENERICEQUATION_H_
#define EQUATIONS_GENERICEQUATION_H_

#include <functional>
#include <set>
#include <string>
#include <utility>

#include "Data/Field.h"
#include "MPI/ExecutionManager.h"
#include "MatrixSystem/Trilinos.h"
#include "MatrixSystem/TrilinosSolver.h"
#include "Utilities/Observer.h"
#include "Algorithm/AlgorithmTraits.h"

namespace dare {

struct UpdateFieldCopy {
    template<typename EQ, typename MS>
    void operator()(EQ& eq, const MS& ms) const {   // NOLINT
        ms.CopyTo(eq.GetField()->GetDataVector());
    }
};

struct UpdateFieldAddInto {
    template <typename EQ, typename MS>
    void operator()(EQ& eq, const MS& ms) const {  // NOLINT
        ms.AddTo(eq.GetField()->GetDataVector());
    }
};

template <typename Grid, typename BoundaryStrategy, typename CustomMember>
class GenericEquation {
public:
    using GridType = Grid;
    using GridRepresentation = typename GridType::Representation;
    using BoundaryStrategyType = BoundaryStrategy;
    using SC = typename GridType::ScalarType;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using FieldType = Field<GridType, SC, 1>;
    using CustomMemberType = CustomMember;
    using MatrixSystemType = dare::Trilinos<SC>;            // could be made a template parameter
    using MatrixSolverType = dare::TrilinosSolver<SC>;      // could be made a template parameter
    using SolverNumericalPropertiesType = typename MatrixSolverType::NumericalPropertiesType;
    using SelfType = GenericEquation<Grid, BoundaryStrategy, CustomMember>;

    GenericEquation(const std::string& name,
                    GridRepresentation grid,
                    std::size_t num_tsteps,
                    BoundaryStrategy bc_strat);

    explicit GenericEquation(const SelfType&) = delete;
    SelfType& operator=(const SelfType&) = delete;

    void PreStep();

    template <typename BuildStrategy>
    void Build(BuildStrategy build);

    template <typename BuildStrategy>
    void UpdateRhs(BuildStrategy build);

    template <typename Lambda>
    std::pair<bool, int> Solve(Lambda UpdateStrategy =
                                   UpdateFieldCopy{});

    void UpdateBoundaries();

    void PostStep();

    GridRepresentation* GetGridRepresentation();
    const GridRepresentation& GetGridRepresentation() const;

    BoundaryStrategyType* GetBoundaryStrategy();
    const BoundaryStrategyType& GetBoundaryStrategy() const;

    FieldType* GetField();
    const FieldType& GetField() const;

    CustomMemberType* GetCustomMember();
    const CustomMemberType& GetCustomMember() const;

    void SetSolverNumericalProperties(const SolverNumericalPropertiesType& prop);
    const SolverNumericalPropertiesType& GetSolverNumericalProperties() const;

    void AddPreStepStrategy(std::function<void(SelfType*)> f);
    void SetPreStepStrategy(std::function<void(SelfType*)> f);
    void ClearPreStepStrategy();
    void AddPostStepStrategy(std::function<void(SelfType*)> f);
    void SetPostStepStrategy(std::function<void(SelfType*)> f);
    void ClearPostStepStrategy();

private:
    GridRepresentation grep;
    dare::ExecutionManager* exec_man;
    FieldType field;
    BoundaryStrategyType boundary_strategy;
    CustomMemberType custom_member;
    std::set<std::function<void(SelfType*)>> pre_step_strategy;
    std::set<std::function<void(SelfType*)>> post_step_strategy;
    MatrixSystemType matrix_system;
    SolverNumericalPropertiesType solver_prop;
};

}  // namespace dare

#include "GenericEquation.inl"

#endif  // EQUATIONS_GENERICEQUATION_H_
