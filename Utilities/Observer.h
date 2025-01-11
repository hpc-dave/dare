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

#ifndef UTILITIES_OBSERVER_H_
#define UTILITIES_OBSERVER_H_

#include <functional>
#include <utility>

namespace dare::utils {

/*!
 * \brief a generic observer for implementation of the Observer pattern
 * @tparam Subject type of object which is observed
 * @tparam StateTag a tag which allows differentiating between different callbacks
 */
template <typename Subject, typename StateTag>
class Observer {
public:
    using OnUpdate = std::function<void(Subject const&, StateTag)>;

    /*!
     * @brief constructor
     * @param onUpdate function pointer to update function
     * 
     * The signature of the update function has to be void(const Subject&, StateTag)
     */
    explicit Observer(OnUpdate onUpdate)
        : onUpdate_{std::move(onUpdate)} {
        // Possibly respond on an invalid/empty std::function instance
    }

    /*!
     * @brief update function called by the subject to notify a change
     * @param subject reference to the subject
     * @param property the state that was changed
     */
    void update(const Subject& subject, StateTag property) {
        onUpdate_(subject, property);
    }

private:
    OnUpdate onUpdate_;     //!< function pointer with the update function
};

}  // namespace dare::utils

#endif  // UTILITIES_OBSERVER_H_
