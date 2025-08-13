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
#include <memory>

namespace dare {

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
    void Update(const Subject& subject, StateTag property) {
        onUpdate_(subject, property);
    }

private:
    OnUpdate onUpdate_;     //!< function pointer with the update function
};

template <typename T>
concept Observable =
    requires {
        typename T::StateChange;
    } && requires(T t, Observer<T, typename T::StateChange>* o, typename T::StateChange s) {
        { t.Attach(o) } -> std::same_as<bool>;
        { t.Detach(o) } -> std::same_as<bool>;
        t.Notify(s);
    };  // NOLINT

namespace detail {
/*!
 * @brief a concept class for external polymorphism to handle observers
 *
 * The class is empty, because we are only interested in the
 * automatically generated virtual destructor. Otherwise the
 * observer is not called by the owning object, only
 * by the observed object, and that one knows the type!
 */
class ObserverHandleConcept {
};

}  // end namespace detail

using UniqueObserverHandle = std::unique_ptr<detail::ObserverHandleConcept>;

/*!
 * @brief an observer handle for managing the lifetime of observers
 * @tparam T an observer type
 * No functions except the constructor are required, since the owning
 * instance of an observer does not require any further direct access.
 */
template <typename Subject, typename StateTag>
class ObserverHandleModel : public detail::ObserverHandleConcept {
public:
    using ObserverType = Observer<Subject, StateTag>;
    explicit ObserverHandleModel(ObserverType&& observer) : obs(std::move(observer)) {}

    /*!
     * @brief a dedicated function for attaching the observer to the observable object
     * @return address of the observer
     */
    ObserverType* GetObserver() { return &obs; }

private:
    ObserverType obs;  //!< actual instance of the observer
};

/*!
 * @brief A convenience function to get an observer handle
 * @tparam Lambda Update function type
 * @tparam T the observable type
 * @param on_update actual update function
 *
 * First an observer is instantiated, which then is wrapped in an observer model.
 * The address of the observer is attached to the observable object and finally
 * the model returned as basic concept (or handle)
 */
template<Observable T, typename Lambda>
UniqueObserverHandle make_observer_handle(T* observable, Lambda on_update) {
    using ObserverType = T::ObserverType;
    using TObsModel = dare::ObserverHandleModel<T, typename T::StateChange>;

    ObserverType obs(on_update);
    auto obs_model = std::make_unique<TObsModel>(std::move(obs));
    observable->Attach(obs_model->GetObserver());
    return obs_model;
}

}  // namespace dare

#endif  // UTILITIES_OBSERVER_H_
