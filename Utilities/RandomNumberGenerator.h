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

#ifndef UTILITIES_RANDOMNUMBERGENERATOR_H_
#define UTILITIES_RANDOMNUMBERGENERATOR_H_

#include <random>


namespace dare {

template<typename Target, typename Generator, typename Distribution>
class PseudoRandomNumberGenerator {
public:
    using GeneratorType = Generator;
    using DistributionType = Distribution;
    using TargetType = Target;
    using ParameterType = typename DistributionType::param_type;
    using ResultType = typename DistributionType::result_type;

    /*!
     * @brief construcor for the Pseudo random numer generator
     * @tparam Args types for initializing the distribution type
     * @param args arguments for the distribution
     */
    template<typename... Args>
    explicit PseudoRandomNumberGenerator(Args&&... args): factor(1), distribution(args...) {
    }

    void SetPreFactor(TargetType f) {
        factor = f;
    }

    TargetType Generate() {
        return factor * distribution(generator);
    }

private:
    TargetType factor;
    GeneratorType generator;
    DistributionType distribution;
};


template <typename Target>
using PseudoRandomTGenerator = PseudoRandomNumberGenerator<Target,
                                        std::default_random_engine,
                                        std::uniform_int_distribution<int>>;

using PseudoRandomIntGenerator = PseudoRandomTGenerator<int>;
using PseudoRandomDoubleGenerator = PseudoRandomTGenerator<double>;
}  // namespace dare

#endif  // UTILITIES_RANDOMNUMBERGENERATOR_H_
