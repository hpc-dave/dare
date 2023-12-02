/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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

#ifndef MATRIXSYSTEM_BICGSTAB2_H_
#define MATRIXSYSTEM_BICGSTAB2_H_

#include "BelosMultiVecTraits.hpp"
#include "BelosSolverManager.hpp"

namespace dare::Matrix {

/*!
 * \class BiCGSTAB2
 * \brief BiConjugate Gradient Stabilized (2) method
 * \tparam SC scalar type
 * \tparam MV multivector type
 * \tparam OperatorType
 *
 * BiConjugate Gradients Stabilized (2) method works with both symmetric and non-symmetric matrices.
 * This method is stable and is an improved version of BiCGSTAB.
 *
 * Was proposed by Sleijpen G, Fokkema D., "BICGStab(L) for linear equations involving unsymmetric matrices with complex spectrum",
 * 1993, ETNA, Vol. 1, pp. 11-32.
 *
 * The algorithm has been taken from H. van der Vorst, "Iterative Krylov Methods for Large Linear Systems",
 * 2003, Cambridge Monographs on Applied and Computational Mathematics, vol. 13.
 *
 * \note In general it is just an unrolled algorithm from Sleijpen and Fokkema
 *
 * \note Currently, the implementation is just an updated version for Tpetra of the Epetra based hybrid implementation
 * in the BumbleBeeMPI submodule, implemented by Maxim Masterov
 *
 * To developers:
 * Tpetra supports multivector operations, however this implementation is not yet adapted to such operations. If you
 * target to adapt the solver, have a look at the functions provided in <b>BelosMultiVecTraits.hpp</b>
 */
template <typename SC, typename MV, typename OP>
class BiCGStab2 : public Belos::SolverManager<ST, MV, OP> {
    using ScalarType = SC;
    using MultiVectorType = MV;
    using OperatorType = OP;
    using LocalOrdinal = typename MV::local_ordinal_type;
    using GlobalOrdinal = typename MV::global_ordinal_type;
    using Node = typename MV::node_type;
    using STS = Teuchos::ScalarTraits<ST>;
    using Magnitude = typename STS::magnitudeType;
    using Problem = Belos::LinearProblem<ST, MV, OP>;
    using MapType = typename MultiVectorType::map_type;

    enum ConvergenceCriteria {
        RNORM,   //!< L2-norm of residual vec
        RBNORM,  //!< L2-norm of residual vec normalized by L2-norm of right hand side
        RWNORM,  //!< Weighted L2-norm of residual vec
        RRNORM,  //!< L2-norm of residual vec normalized by L2-norm of initial residual
        INTERN   //!< Default criteria determined for each solver separately (see implementation)
    };

    //! @name Constructors/Destructor
    //@{
    explicit BiCGStab2(const Teuchos::RCP<Teuchos::ParameterList>& params);
    virtual ~BiCGStab2();

    /// \brief clone the solver manager.
    ///
    /// Implements the DII inversion and injection pattern
    Teuchos::RCP<Belos::SolverManager<ST, MV, OP> > clone() const override;
    //@}

    //! @name Accessor methods
    //@{
    //! Return a reference to the linear problem being solved by this solver manager.
    const Belos::LinearProblem<ScalarType, MV, OP>& getProblem() const override;

    //! Return the valid parameters for this solver manager.
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

    //! Return the current parameters being used for this solver manager.
    Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const override;

    /// \brief Tolerance achieved by the last \c solve() invocation.
    ///
    /// This is the maximum over all right-hand sides' achieved
    /// convergence tolerances, and is set whether or not the solve
    /// actually managed to achieve the desired convergence tolerance.
    ///
    /// The default implementation throws std::runtime_error.  This is
    /// in case the idea of a single convergence tolerance doesn't make
    /// sense for some solvers.  It also serves as a gradual upgrade
    /// path (since this method is a later addition to the \c
    /// SolverManager interface).
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType achievedTol() const;

    //! Get the iteration count for the most recent call to \c solve().
    int getNumIters() const override;

    /*! \brief Returns whether a loss of accuracy was detected in the solver.
     *  \note This method is normally applicable to GMRES-type solvers.
     */
    bool isLOADetected() const override;

    //@}

    //! @name Set methods
    //@{

    //! Set the linear problem that needs to be solved.
    void setProblem(const Teuchos::RCP<Belos::LinearProblem<ScalarType, MV, OP> >& problem) override;

    /// \brief Set the parameters to use when solving the linear problem.
    ///
    /// \param params [in/out] List of parameters to use when solving
    ///   the linear problem.  This list will be modified as necessary
    ///   to include default parameters that need not be provided.  If
    ///   params is null, then this method uses default parameters.
    ///
    /// \note The ParameterList returned by \c getValidParameters() has
    ///   all the parameters that the solver understands, possibly
    ///   including human-readable documentation and validators.
    void setParameters(const Teuchos::RCP<Teuchos::ParameterList>& params) override;

    //! Set user-defined convergence status test.
    virtual void setUserConvStatusTest(
        const Teuchos::RCP<Belos::StatusTest<ST, MV, OP> >& userConvStatusTest,
        const typename Belos::StatusTestCombo<ST, MV, OP>::ComboType& comboType =
            Belos::StatusTestCombo<ST, MV, OP>::SEQ);

    //! Set user-defined debug status test.
    virtual void setDebugStatusTest(const Teuchos::RCP<Belos::StatusTest<ScalarType, MV, OP> >& debugStatusTest);

    //@}

    //! @name Reset methods
    //@{

    /// \brief Reset the solver manager.
    ///
    /// Reset the solver manager in a way specified by the \c
    /// ResetType parameter.  This informs the solver manager that the
    /// solver should prepare for the next call to solve by resetting
    /// certain elements of the iterative solver strategy.
    void reset(const Belos::ResetType type) override;
    //@}

    //! @name Solver application methods
    //@{

    /// \brief Iterate until the status test tells us to stop.
    //
    /// This method performs possibly repeated calls to the underlying
    /// linear solver's iterate() routine, until the problem has been
    /// solved (as decided by the solver manager via the status
    /// test(s)), or the solver manager decides to quit.
    ///
    /// \return A \c Belos::ReturnType enum specifying:
    ///   - Belos::Converged: the linear problem was solved to the
    ///     specification required by the solver manager.
    ///   - Belos::Unconverged: the linear problem was not solved to the
    ///     specification desired by the solver manager.
    Belos::ReturnType solve() override;
    //@}

private:
    /*!
     * \brief solves the matrix system
     * @param A matrix
     * @param x solution vector (also initial guess)
     * @param B right hand side vector
     */
    bool Solve(const OP& A, MV& x, const MV& B);  // NOLINT

    Teuchos::RCP<const OP> A;       //!< pointer to matrix
    Teuchos::RCP<const MV> B;       //!< pointer to rhs
    Teuchos::RCP<MV> x;             //!< pointer to solution vector
    Teuchos::RCP<Problem> problem;  //!< pointer to linear porblem formulation

    Magnitude tol;      //!< convergence tolerance
    int max_iteration;  //!< maximum number of iterations
    int verbosity;      //!< level of verbosity

    Magnitude norm_residual;  //!< residual of the latest solving
    int num_iteration;        //!< number of iterations of latest solving
    int stop_criteria;        //!< identifier of stop criteria
    bool converged;           //!< identifier, if converged
};
}  // end namespace dare::Matrix

#include "BiCGStab2.inl"

#endif  // MATRIXSYSTEM_BICGSTAB2_H_
