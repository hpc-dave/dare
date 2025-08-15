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

#ifndef ALGORITHM_ADAPTIVETIMESTEP_H_
#define ALGORITHM_ADAPTIVETIMESTEP_H_

#include <set>
#include "Data/DefaultTypes.h"
#include "Algorithm/AlgorithmTraits.h"
#include "Utilities/Observer.h"
namespace dare {

/*! \class AdaptiveTimeStep
 * @tparam T type of the time values
 *
 * This class controls timestepping and allows adaption of the timestep
 */
template<typename T = dare::defaults::ScalarType>
class AdaptiveTimeStep {
public:
    enum class StateChange {
        AdvanceTimeStep,
        UpdateTimeStepSize
    };
    using ValueType = T;    //!< type used for keeping track of the time
    using CounterType = TimeStepCounter;    //!< integer value for the current time step number
    using ObserverType = dare::Observer<AdaptiveTimeStep<T>, StateChange>;  //!< the type of the observer

    /*!
     * @brief main constructor
     * @param _dt time step size
     * @param current_time the time at the specified time step
     * @param counter_begin the time step to start with
     */
    explicit AdaptiveTimeStep(ValueType _dt, ValueType current_time = 0., CounterType counter_begin)
        : dt{_dt}, time{current_time}, counter{counter_begin} {}

    /*!
     * @brief returns the time step size
     */
    ValueType GetTimeStepSize() const {
        return dt;
    }

    /*!
     * @brief returns the current simulation time
     */
    ValueType GetTime() const {
        return time;
    }

    /*!
     * @brief returns the number of conducted timesteps
     */
    TimeStepCounter GetTimeStepCounter() const {
        return counter;
    }

    /*!
     * @brief advances to the next timestep and notfies all observers
     */
    void AdvanceTimeStep() {
        time += dt;
        ++counter;
        Notify(StateChange::AdvancetimeStep);
    }

    /*!
     * @brief attaches an observer
     * @param o pointer to the observer
     * @return true if successful
     */
    bool Attach(ObserverType* o) {
        auto [pos, success] = observers.emplace(o);
        return success;
    }

    /*!
     * @brief detaches an observer
     * @param o pointer to the observer to detach
     * @return true if the observer was found and detached
     */
    bool Detach(ObserverType* o) {
        return (observers.erase(o) > 0U);
    }

    /*!
     * @brief notifies all observers of a state change
     * @param property the state that has changed
     */
    void Notify(StateChange property) {
        for (auto iter = observers.begin(); iter != observers.end();) {
            auto const pos = iter++;
            (*pos)->Update(*this, property);
        }
    }

    /*!
     * @brief adapts the timestep size and notifies the attached observers
     * @param _dt value to adapt the timestep to
     */
    void AdaptTimeStepSize(ValueType _dt) {
        dt = _dt;
        Notify(StateChange::Update);
    }

    /*!
     * @brief a conversion function which allows using the observer as if it were of ValueType
     */
    operator ValueType() { return dt; }

private:
    ValueType dt;                       //!< the time step size
    ValueType time;                     //!< current simulation time
    CounterType counter;                //!< time step counter
    std::set<ObserverType*> observers;  //!< attached observers
};
}  // namespace dare

// static test
static_assert(dare::TimeStepper<AdaptiveTimeStep<double>>, "does not fullfill requirements");

#endif  // ALGORITHM_ADAPTIVETIMESTEP_H_
