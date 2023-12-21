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

namespace dare::Matrix {

template <class ST, class MV, class OP>
BiCGStab2<ST, MV, OP>::BiCGStab2(const Teuchos::RCP<Teuchos::ParameterList>& params)
    : tol(1e-8),
      max_iteration(100),
      verbosity(0),
      norm_residual(1e-14),
      num_iteration(0),
      stop_criteria(RNORM),
      converged(false) {
    setParameters(params);
}

template <class ST, class MV, class OP>
BiCGStab2<ST, MV, OP>::~BiCGStab2() {
}

template <class ST, class MV, class OP>
Teuchos::RCP<Belos::SolverManager<ST, MV, OP>> BiCGStab2<ST, MV, OP>::clone() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "clone() not implemented");
    Teuchos::RCP<Teuchos::ParameterList> params;
    params->set("Convergence Tolerance", tol);
    params->set("Maximum Iterations", max_iteration);
    params->set("Verbosity", verbosity);
    return Teuchos::rcp(new BiCGStab2<ST, MV, OP>(params));
}

template <class ST, class MV, class OP>
const Belos::LinearProblem<ST, MV, OP>& BiCGStab2<ST, MV, OP>::getProblem() const {
    return *problem;
}

template <class ST, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList> BiCGStab2<ST, MV, OP>::getValidParameters() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "getValidParameters() not implemented");
    return rcp(new Teuchos::ParameterList());
}

template <class ST, class MV, class OP>
Teuchos::RCP<const Teuchos::ParameterList> BiCGStab2<ST, MV, OP>::getCurrentParameters() const {
    Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
    params->set("Convergence Tolerance", tol);
    params->set("Maximum Iterations", max_iteration);
    params->set("Verbosity", verbosity);
    return rcp(new Teuchos::ParameterList());
}

template <class ST, class MV, class OP>
typename Teuchos::ScalarTraits<ST>::magnitudeType BiCGStab2<ST, MV, OP>::achievedTol() const {
    return tol;
}

template <class ST, class MV, class OP>
int BiCGStab2<ST, MV, OP>::getNumIters() const {
    return num_iteration;
}

template <class ST, class MV, class OP>
bool BiCGStab2<ST, MV, OP>::isLOADetected() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "isLOADDetected() not implemented");
    return false;
}

template <class ST, class MV, class OP>
void BiCGStab2<ST, MV, OP>::setProblem(const Teuchos::RCP<Belos::LinearProblem<ScalarType, MV, OP>>& _problem) {
    problem = _problem;
}

template <class ST, class MV, class OP>
void BiCGStab2<ST, MV, OP>::setParameters(const Teuchos::RCP<Teuchos::ParameterList>& params) {
    tol = Teuchos::get<Magnitude>(*params, "Convergence Tolerance");
    max_iteration = Teuchos::get<int>(*params, "Maximum Iterations");
    try {
        std::string verb = Teuchos::get<std::string>(*params, "Verbosity");
        if (verb.compare("none") == 0) {
            verbosity = 0;
        } else if (verb.compare("low") == 0) {
            verbosity = 1;
        } else if (verb.compare("medium") == 0) {
            verbosity = 2;
        } else if (verb.compare("high") == 0) {
            verbosity = 3;
        }
    } catch (Teuchos::Exceptions::InvalidParameterName& b) {
    }
}

template <class ST, class MV, class OP>
void BiCGStab2<ST, MV, OP>::setUserConvStatusTest(
    const Teuchos::RCP<Belos::StatusTest<ST, MV, OP>>& userConvStatusTest,
    const typename Belos::StatusTestCombo<ST, MV, OP>::ComboType& comboType) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                "Error, the function setUserConvStatusTest() has not been"
                                    << " overridden for the class" << this->description() << " yet!");
}

template <class ST, class MV, class OP>
void BiCGStab2<ST, MV, OP>::setDebugStatusTest(
    const Teuchos::RCP<Belos::StatusTest<ST, MV, OP>>& debugStatusTest) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true,
                                "Error, the function setDebugStatusTest() has not been"
                                    << " overridden for the class"
                                    << this->description() << " yet!");
}

template <class ST, class MV, class OP>
void BiCGStab2<ST, MV, OP>::reset(const Belos::ResetType type) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "reset() not implemented");
}

template <class ST, class MV, class OP>
Belos::ReturnType BiCGStab2<ST, MV, OP>::solve() {
    using Teuchos::RCP;
    TEUCHOS_TEST_FOR_EXCEPTION(problem.is_null(), std::runtime_error,
                               "The linear problem has not "
                               "yet been set.  Please call setProblem with a nonnull argument before "
                               "calling this method.");

    B = problem->getRHS();
    TEUCHOS_TEST_FOR_EXCEPTION(B.is_null(), std::runtime_error,
                               "The linear problem's right-hand "
                               "side(s) B has/have not yet been set.  Please call setProblem with "
                               "a nonnull argument before calling this method.");

    x = problem->getLHS();
    TEUCHOS_TEST_FOR_EXCEPTION(x.is_null(), std::runtime_error,
                               "The linear problem's left-hand "
                               "side(s) X has/have not yet been set.  Please call setProblem with "
                               "a nonnull argument before calling this method.");

    A = problem->getOperator();
    TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::runtime_error,
                               "The linear problem's Matrix A "
                               "has/have not yet been set.  Please call setProblem with "
                               "a nonnull argument before calling this method.");
    converged = Solve(*A, *x, *B);

    return converged ? Belos::ReturnType::Converged : Belos::ReturnType::Unconverged;
}

template <class ST, class MV, class OP>
bool BiCGStab2<ST, MV, OP>::Solve(const OP& A, MV& x, const MV& B) {  // NOLINT
    /*
     * For developers: Have a look at "BelosMultiVectorTraits.hpp"
     * The functions there provide quite convenient possibilities for
     * introducing multivectors instead of the singular implementation
     * below
     */
    TEUCHOS_TEST_FOR_EXCEPTION(x.getNumVectors() != B.getNumVectors(),
                               std::runtime_error,
                               "x.getNumVectors() = " << x.getNumVectors()
                                                      << " != B.getNumVectors() = " << B.getNumVectors() << ".");

    using Teuchos::RCP;

    int k = 0;                   // iteration number
    Magnitude alpha = 0.;        // part of the method
    Magnitude rho[2] = {0.};     // part of the method
    Magnitude gamma = 0.;        // part of the method
    Magnitude beta = 0.;         // part of the method
    long double omega_1 = 0.0L;  // part of the method, stored as a long to prevent overflow
    long double omega_2 = 0.0L;  // part of the method, stored as a long to prevent overflow
    long double mu = 0.0L;       // part of the method, stored as a long to prevent overflow
    long double nu = 0.0L;       // part of the method, stored as a long to prevent overflow
    long double tau = 0.0L;      // part of the method, stored as a long to prevent overflow

    Magnitude r_norm_0 = 0.;           // Preconditioned norm
    Magnitude convergence_check = 0.;  // keeps new residual
    Magnitude normalizer = 1.;         // normalizer for the residual
    const Magnitude ONE = STS::one();

    const bool is_left_prec{problem->isLeftPrec()};    // identifier, if the problem is left conditioned
    const bool is_right_prec{problem->isRightPrec()};  // identifier, if the problem is right conditioned

    /*
     * MPI communicators
     */
    RCP<const MapType> map = x.getMap();
    RCP<const Teuchos::Comm<int>> comm = map->getComm();

    const int my_rank = comm->getRank();
    const bool is_root{my_rank == 0};  // identifier, if the process is root
    if (!std::is_same<double, Magnitude>()) {
        if (is_root)
            std::cerr << "Warning, scalar type of system provided to BiCGStab2 is NOT a 'double'!"
                      << " Make sure the scalar type is compatible!\n";
    }

    // Anonymous function for error printing
    auto PrintError = [&](int step, std::string what, std::map<std::string, double> values) {
        if (my_rank != 0)
            return;

        if (verbosity > 0)
            std::cout << "BiCGSTAB(2) has been interrupted due to " << what
                      << " @ step " << std::to_string(step) << std::endl;

        if (verbosity > 1)
            for (const auto& [key, value] : values)
                std::cout << " " << key << std::scientific << value << std::endl;
    };

    const std::size_t num_vec = x.getNumVectors();  // number of vectors (columns)

    Kokkos::View<Magnitude*, Kokkos::HostSpace> norms("norms", num_vec);                  // Kokkos view for norms
    Kokkos::View<typename MV::dot_type*, Kokkos::HostSpace> dot_product("dot", num_vec);  // Kokkos view for dot prod

    if (num_vec != 1) {
        if (is_root)
            std::cerr << "Warning, BiCGStab2 has been tested for single vectors only!"
                      << " Use with multiple vectors may lead to unexpected results\n";
    }

    // FIXME: Possibility for optimization? Storing those as part of the class?
    MultiVectorType x0(map, num_vec);
    MultiVectorType r(map, num_vec);
    MultiVectorType r_hat_0(map, num_vec);
    MultiVectorType u(map, num_vec);
    MultiVectorType v(map, num_vec);
    MultiVectorType s(map, num_vec);
    MultiVectorType w(map, num_vec);
    MultiVectorType t(map, num_vec);
    MultiVectorType u_hat(map, num_vec);
    MultiVectorType r_hat(map, num_vec);
    MultiVectorType v_hat(map, num_vec);
    MultiVectorType s_hat(map, num_vec);
    MultiVectorType tmp(map, num_vec);

    x0.assign(x);
    r.assign(B);

    // Right preconditioner
    if (is_right_prec) {
        problem->applyRightPrec(x0, tmp);
        x0.assign(tmp);
    }

    //! (0) \f$ r = \hat{r}_0 = b - A * x_0\f$
    //    r = (b - Matrix * x0);
    A.apply(x0, v);
    r.update(-ONE, v, ONE);
    r_hat_0.assign(r);

    //! (1) \f$ u = 0 \f$, \f$ w = 0 \f$, \f$ v = 0 \f$, \f$ \alpha = \rho[0] = \omega_1 = \omega_2 = 1\f$
    u.putScalar(0.);
    v.putScalar(0.);
    w.putScalar(0.);
    rho[0] = omega_1 = omega_2 = ONE;

    //! (2) Solve \f$ M y = r \f$, set \f$ r = y \f$
    if (is_left_prec) {
        // Case of left preconditioner
        problem->applyLeftPrec(r, tmp);
        r.assign(tmp);
    }

    r.norm2(norms);
    r_norm_0 = norms[0];  // TODO(dave): fix for multiple vectors

    /*!
     * Prepare stop criteria
     */
    switch (stop_criteria) {
    case RNORM:
        normalizer = 1.;
        break;
    case RRNORM:
    case INTERN:
        normalizer = r_norm_0;
        break;
    case RBNORM:
        B.norm2(norms);
        normalizer = norms[0];  // TODO(dave): adapt for multiple vectors
        break;
    case RWNORM:
        normalizer = 1.;
        break;
    default:
        normalizer = 1.;
        break;
    }

    /*
     * Check residual. Stop if initial guess satisfies convergence criteria.
     */
    convergence_check = r_norm_0 / normalizer;
    if (convergence_check < tol) {
        num_iteration = k;
        norm_residual = convergence_check;
        x.assign(x0);
        return true;
    }

    ++k;

    //! Start iterative loop
    while (1) {
        if (k > max_iteration) {
            break;
        }

        //! (3) \f$ \rho[0] = - \omega_2 \rho[0] \ \f$
        rho[0] = -omega_2 * rho[0];
        if (rho[0] == 0.0) {
            if (is_root)
                PrintError(3, "rho[0] == 0.0", {{"omega2: ", omega_2}});
            break;
        }

        /*!
         * Even Bi-CG step
         */
        //! (4) \f$ \rho[1] = <\hat{r}_0, r>\f$, \f$ \beta = \alpha \rho[1] / \rho[0] \f$, \f$ \rho[0] = \rho[1] \f$
        r_hat_0.dot(r, dot_product);
        rho[1] = dot_product[0];  // FIXME: adapt for mulitvectors
        beta = alpha * rho[1] / rho[0];
        rho[0] = rho[1];
        if (rho[0] == 0.0) {
            if (is_root)
                PrintError(4, "rho[0] == 0.0",
                          {{"rho[0]: ", rho[0]}, {"rho[1]: ", rho[1]}, {"alpha:  ", alpha}, {"beta:   ", beta}});
            break;
        }

        //! (5) \f$ u = r - \beta u \f$
        u.update(ONE, r, -beta);

        if (is_right_prec) {
            // Case of right preconditioner
            //! (6) \f$ v = A M^{-1} u \f$
            problem->applyRightPrec(u, tmp);
            A.apply(tmp, v);
        } else if (is_left_prec) {
            // Case of left preconditioner
            //! (6) \f$ v = M^{-1} A u \f$
            A.apply(u, tmp);
            problem->applyLeftPrec(tmp, v);
        } else {
            A.apply(u, v);
        }

        //! (7) \f$ \gamma = <v, \hat{r}_0> \f$, \f$ \alpha = \rho[0] / \gamma \f$
        v.dot(r_hat_0, dot_product);
        gamma = dot_product[0];  // FIXME: adapt for multivectors

        // Check for breakdown (probably may occur)
        if (gamma == 0.0) {
            v.norm2(norms);
            Magnitude norm_v = norms[0];
            r_hat_0.norm2(norms);
            Magnitude norm_r_hat_0 = norms[0];
            if (is_root)
                PrintError(7, "gamma == 0.0", {{"norm v:       ", norm_v}, {"norm r_hat_0: ", norm_r_hat_0}});
            break;
        }

        alpha = rho[0] / gamma;

        //! (8) \f$ r = r - \alpha v \f$
        r.update(-alpha, v, ONE);

        if (is_right_prec) {
            // Case of right preconditioner
            //! (9) \f$ s = A M^{-1} r \f$
            problem->applyRightPrec(r, tmp);
            A.apply(tmp, s);
        } else if (is_left_prec) {
            // Case of left preconditioner
            //! (9) \f$ s = M^{-1} A r \f$
            A.apply(r, tmp);
            problem->applyLeftPrec(tmp, s);
        } else {
            A.apply(r, s);
        }

        //! (10) \f$ x = x + \alpha u \f$
        x.update(alpha, u, 1.);

        /*!
         * Odd Bi-CG step
         */
        //! (11) \f$ \rho[1] = <\hat{r}_0, s>\f$, \f$ \beta = \alpha \rho[1] / \rho[0] \f$, \f$ \rho[0] = \rho[1] \f$
        r_hat_0.dot(s, dot_product);
        rho[1] = dot_product[0];
        beta = alpha * rho[1] / rho[0];
        rho[0] = rho[1];

        //! (12) \f$ v = s - \beta v \f$
        v.update(ONE, s, -beta);

        if (is_right_prec) {
            // Case of right preconditioner
            //! (13) \f$ w = A M^{-1} v \f$
            problem->applyRightPrec(v, tmp);
            A.apply(tmp, w);
        } else if (is_left_prec) {
            // Case of left preconditioner
            //! (13) \f$ w = M^{-1} A v \f$
            A.apply(v, tmp);
            problem->applyLeftPrec(tmp, w);
        } else {
            A.apply(v, w);
        }

        //! (14) \f$ \gamma = <w, \hat{r}_0> \f$, \f$ \alpha = \rho[0] / \gamma \f$
        w.dot(r_hat_0, dot_product);
        gamma = dot_product[0];

        // Check for breakdown (probably may occur)
        if (gamma == 0.0) {
            u.norm2(norms);
            double norms_u = norms[0];
            s.norm2(norms);
            double norms_s = norms[0];
            v.norm2(norms);
            double norms_v = norms[0];
            w.norm2(norms);
            double norms_w = norms[0];
            r_hat_0.norm2(norms);
            double norms_r_hat_0 = norms[0];
            if (is_root)
                PrintError(14, "gamma == 0.0",
                            {{"alpha:        ", alpha},
                             {"beta:         ", beta},
                             {"norm u:       ", norms_u},
                             {"norm s:       ", norms_s},
                             {"norm v:       ", norms_v},
                             {"norm w:       ", norms_w},
                             {"norm r_hat_0: ", norms_r_hat_0}});
            break;
        }

        alpha = rho[0] / gamma;

        //! (15) \f$ u = r - \beta u \f$
        u.update(1., r, -beta);

        //! (16) \f$ r = r - \alpha v \f$
        r.update(-alpha, v, ONE);

        //! (17) \f$ s = s - \alpha w\f$
        s.update(-alpha, w, ONE);

        if (is_right_prec) {
            // Case of right preconditioner
            //! (18) \f$ t = A M^{-1} s\f$
            problem->applyRightPrec(s, tmp);
            A.apply(tmp, t);
        } else if (is_left_prec) {
            // Case of left preconditioner
            //! (18) \f$ t = M^{-1} A s\f$
            A.apply(s, tmp);
            problem->applyLeftPrec(tmp, t);
        } else {
            A.apply(s, t);
        }

        /*
         * GCR(2)-part
         */
        //! (19) \f$ \omega_1 = <r, s> \f$, \f$ \mu = <s, s> \f$, \f$ \nu = <s, t> \f$, \f$ \tau = <t, t> \f$
        //! (20) \f$ \omega_2 = <r, t> \f$
        r.dot(s, dot_product);
        omega_1 = dot_product[0];

        s.dot(s, dot_product);
        mu = dot_product[0];

        s.dot(t, dot_product);
        nu = dot_product[0];

        t.dot(t, dot_product);
        tau = dot_product[0];

        r.dot(t, dot_product);
        omega_2 = dot_product[0];

        if (mu == 0.) {
            s.norm2(norms);
            double norm_s = norms[0];
            if (is_root)
                PrintError(20, "mu == 0.0", {{"norm s: ", norm_s}});
            break;
        }

        //! (21) \f$ \tau = \tau - \nu^2 / \mu \f$
        tau -= nu * nu / mu;

        if (tau == 0.0) {
            // it may happen, that the norm is already reached by luck and tau is subsequently 0
            if (is_right_prec) {
                problem->applyRightPrec(x, tmp);
                x.assign(tmp);
            }
            r.assign(B);
            A.apply(x, v);
            r.update(-1., v, 1.);
            r.norm2(norms);
            convergence_check = norms[0] / normalizer;
            if (convergence_check <= tol)
                return true;

            s.norm2(norms);
            double norm_s = norms[0];
            t.norm2(norms);
            double norm_t = norms[0];

            if (is_root)
                PrintError(21, "tau == 0.0", {{"nu:     ", nu}, {"norm s: ", norm_s}, {"norm t: ", norm_t}});
            break;
        }

        //! (22) \f$ \omega_2 = (\omega_2 - \nu \omega_1 / \mu) / \tau \f$
        omega_2 = (omega_2 - (nu * omega_1) / mu) / tau;

        //! (23) \f$ \omega_1 = (\omega_1 - \nu \omega_2) / \mu \f$
        omega_1 = (omega_1 - nu * omega_2) / mu;

        //! (24) \f$ x = x + \omega_1 r + \omega_2 s + \alpha u \f$
        x.update(static_cast<Magnitude>(omega_1), r, static_cast<Magnitude>(omega_2), s, ONE);
        x.update(alpha, u, ONE);

        //! (25) \f$ r = r - \omega_1 s - \omega_2 t \f$
        r.update(static_cast<Magnitude>(-omega_1), s, static_cast<Magnitude>(-omega_2), t, ONE);

        /*!
         * Check convergence
         */
        if (is_left_prec) {
            // Case of left preconditioner
            problem->applyLeftPrec(r, tmp);
            tmp.norm2(norms);
        } else {
            r.norm2(norms);
        }
        convergence_check = norms[0] / normalizer;

        if (convergence_check <= tol && k > 1) {
            break;
        }

        //! (25) \f$ u = u - \omega_1 v - \omega_2 w \f$
        u.update(static_cast<Magnitude>(-omega_1), v, static_cast<Magnitude>(-omega_2), w, ONE);

        ++k;
    }

    if (is_right_prec) {
        problem->applyRightPrec(x, tmp);
        x.assign(tmp);
    }

    num_iteration = k;
    norm_residual = convergence_check;

    return convergence_check <= tol;
}

}  // end namespace dare::Matrix
