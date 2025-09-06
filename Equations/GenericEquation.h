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

/*! \struct UpdateFieldCopy
 * @brief updating strategy for copying data from the matrix system to the field
 */
struct UpdateFieldCopy {
    template<typename EQ, typename MS>
    void operator()(EQ& eq, const MS& ms) const {   // NOLINT
        ms.CopyTo(&eq.GetField()->GetDataVector());
    }
};

/*! \struct UpdateFieldAddInto
 * @brief updating strategy for adding data from the matrix system to the field
 */
struct UpdateFieldAddInto {
    template <typename EQ, typename MS>
    void operator()(EQ& eq, const MS& ms) const {  // NOLINT
        ms.AddTo(&eq.GetField()->GetDataVector());
    }
};

/*!
 * @brief a generic base class which provides a standard interface for implementing equations
 * @tparam Grid type of the utilized grid
 * @tparam BoundaryStrategy object type of the boundary strategy
 * @tparam CustomMember a custom object, that can be added to the class
 */
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

    /*!
     * @brief only constructor
     * @tparam BC type of the input argument that serves as input to the boundary strategy
     * @param name name of the equation, used for creating output and debugging
     * @param grid a representation object of the grid
     * @param num_tsteps number of timesteps that are required for solving the equation
     * @param bc_strat the boundary strategy, which serves as input for the BoundaryStrategy object
     */
    template<typename BC>
    GenericEquation(const std::string& name,
                    GridRepresentation grid,
                    std::size_t num_tsteps,
                    BC bc_strat);

    /*!
     * @brief deleted copy constructor
     */
    explicit GenericEquation(const SelfType&) = delete;

    /*!
     * @brief deleted copy assignment operator
     */
    SelfType& operator=(const SelfType&) = delete;

    /*!
     * @brief a function that is usually executed before solving the equation
     */
    void PreStep();

    /*!
     * @brief builds the matrix system according to the build strategy
     * @tparam BuildStrategy strategy object type which takes a matrix block as parameter
     * @param build the actual strategy
     */
    template <typename BuildStrategy>
    void Build(BuildStrategy build);

    /*!
     * @brief only updates the right hand side (B-vector) of the equation system and leaves the matrix 'as is'
     * @tparam BuildStrategy strategy object type, takes a matrix block as parameter
     * @param build the actual strategy
     */
    template <typename BuildStrategy>
    void UpdateRhs(BuildStrategy build);

    /*!
     * @brief calls the solver with prior set properties and solves the equation
     * @tparam Lambda update strategy that defines, how values are transferred from the matrix system to the field
     * @param UpdateStrategy the update strategy
     * @param build_preconditioner if true, the preconditioner is build prior to solving
     * @return a pair with a boolean and integer, indicating convergence and the number of iterations to reach it
     */
    template <typename Lambda>
    std::pair<bool, int> Solve(Lambda UpdateStrategy =
                                   UpdateFieldCopy{},
                                bool build_preconditioner = true);

    /*!
     * @brief updates the boundaries of the underlying field and exchanges halo cell information
     */
    void UpdateBoundaries();

    /*!
     * @brief additional operations that may be executed after the solving step
     */
    void PostStep();

    /*!
     * @brief provides access to the grid representation
     */
    GridRepresentation* GetGridRepresentation();
    const GridRepresentation& GetGridRepresentation() const;

    /*!
     * @brief provides access to the underlying boundary strategy
     */
    BoundaryStrategyType* GetBoundaryStrategy();
    const BoundaryStrategyType& GetBoundaryStrategy() const;

    /*!
     * @brief provides access to the field
     */
    FieldType* GetField();
    const FieldType& GetField() const;

    /*!
     * @brief provides access to the defined custom member
     */
    CustomMemberType* GetCustomMember();
    const CustomMemberType& GetCustomMember() const;

    /*!
     * @brief allows setting solver properties, e.g. convergence tolerance
     * @param prop the associated property type
     */
    void SetSolverNumericalProperties(const SolverNumericalPropertiesType& prop);
    const SolverNumericalPropertiesType& GetSolverNumericalProperties() const;

    /*!
     * @brief provides access to the matrix system, mainly useful for debugging
     */
    MatrixSystemType* GetMatrixSystem();
    const MatrixSystemType& GetMatrixSystem() const;

    /*!
     * @brief adds a prestep strategy to the existing ones
     * @param f a prestep strategy as function pointer
     */
    void AddPreStepStrategy(std::function<void(SelfType*)> f);

    /*!
     * @brief sets a prestep strategy and overwrites existing ones
     * @param f a prestep strategy as function pointer
     */
    void SetPreStepStrategy(std::function<void(SelfType*)> f);

    /*!
     * @brief removes all prestep strategies
     */
    void ClearPreStepStrategy();

    /*!
     * @brief adds a poststep strategy to the existing ones
     * @param f a poststep strategy as function pointer
     */
    void AddPostStepStrategy(std::function<void(SelfType*)> f);

    /*!
     * @brief sets a poststep strategy and overwrites existing ones
     * @param f a poststep strategy
     */
    void SetPostStepStrategy(std::function<void(SelfType*)> f);

    /*!
     * @brief removes all poststep strategies
     */
    void ClearPostStepStrategy();

private:
    GridRepresentation grep;                                      //!< grid representation
    dare::ExecutionManager* exec_man;                             //!< reference to execution manager
    FieldType field;                                              //!< the field with data
    BoundaryStrategyType boundary_strategy;                       //!< strategy for dealing with boundaries
    CustomMemberType custom_member;                               //!< a custom member
    std::set<std::function<void(SelfType*)>> pre_step_strategy;   //!< set of prestep strategies
    std::set<std::function<void(SelfType*)>> post_step_strategy;  //!< set of postsetp strategies
    MatrixSystemType matrix_system;                               //!< the matrix system Ax=b
    SolverNumericalPropertiesType solver_prop;                    //!< dedicated solver properties
};

}  // namespace dare

#include "GenericEquation.inl"

#endif  // EQUATIONS_GENERICEQUATION_H_
