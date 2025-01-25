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

#ifndef ALGORITHM_CONSTANTTIMESTEP_H_
#define ALGORITHM_CONSTANTTIMESTEP_H_

#include "Data/DefaultTypes.h"
#include "Algorithm/AlgorithmTraits.h"
#include "Utilities/Observer.h"
namespace dare {
template<typename T = dare::defaults::ScalarType>
class ConstantTimeStep {
public:
    enum class StateChange {
        Update
    };
    using ValueType = T;
    using ObserverType = dare::Observer<ConstantTimeStep<T>, StateChange>;

    explicit ConstantTimeStep(ValueType t) : dt(t) {
    }

    virtual ~ConstantTimeStep() {}

    ValueType GetTimeStepSize() const {
        return dt;
    }

    bool Attach(ObserverType* o) {
        // empty, because the time step is constant
        return true;
    }

    bool Detach(ObserverType* o) {
        // empty, because the time step is constant
        return true;
    }

    void Notify() {
        // empty, nothing to notify
    }

private:
    ValueType dt;
};
}  // namespace dare

// static test
static_assert(dare::TimeStepper<dare::ConstantTimeStep<double>>, "does not fullfill requirements");

#endif  // ALGORITHM_CONSTANTTIMESTEP_H_
