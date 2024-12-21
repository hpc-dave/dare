/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder
 *
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

#ifndef EQUATIONS_TIMEDISCRETIZATIONSCHEMES_H_
#define EQUATIONS_TIMEDISCRETIZATIONSCHEMES_H_

namespace dare::Matrix {
struct EULER_BACKWARD {
    static const std::size_t NUM_TIMESTEPS{0};

    template <typename SC>
    static constexpr dare::utils::Vector<NUM_TIMESTEPS + 1, SC> GetWeights() {
        return dare::utils::Vector<NUM_TIMESTEPS + 1, SC>(static_cast<SC>(1.));
    }
};

struct EULER_FORWARD {
    static const std::size_t NUM_TIMESTEPS{1};

    template <typename SC>
    static constexpr dare::utils::Vector<NUM_TIMESTEPS + 1, SC> GetWeights() {
        return dare::utils::Vector<NUM_TIMESTEPS + 1, SC>(static_cast<SC>(0.), static_cast<SC>(1.));
    }
};

struct CRANK_NICHOLSON {
    static const std::size_t NUM_TIMESTEPS{1};

    template <typename SC>
    static constexpr dare::utils::Vector<NUM_TIMESTEPS + 1, SC> GetWeights() {
        return dare::utils::Vector<NUM_TIMESTEPS + 1, SC>(static_cast<SC>(0.5), static_cast<SC>(0.5));
    }
};

struct ADAMS_BASHFORT {
    static const std::size_t NUM_TIMESTEPS{2};

    template <typename SC>
    static constexpr dare::utils::Vector<NUM_TIMESTEPS + 1, SC> GetWeights() {
        return dare::utils::Vector<NUM_TIMESTEPS + 1, SC>(static_cast<SC>(0.),
                                                          static_cast<SC>(1.5),
                                                          static_cast<SC>(-0.5));
    }
};

struct ADAMS_MOULTON {
    static const std::size_t NUM_TIMESTEPS{2};

    template <typename SC>
    constexpr dare::utils::Vector<NUM_TIMESTEPS + 1, SC> GetWeights() {
        return dare::utils::Vector<NUM_TIMESTEPS + 1, SC>(static_cast<SC>(5. / 12.),
                                                          static_cast<SC>(8. / 12.),
                                                          static_cast<SC>(-1. / 12.));
    }
};

}  // namespace dare::Matrix

#endif  // EQUATIONS_TIMEDISCRETIZATIONSCHEMES_H_
