/*
 * MIT License
 *
 * Copyright (c) 2022 David Rieder

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

#ifndef UTILITIES_INITIALIZATIONTRACKER_H_
#define UTILITIES_INITIALIZATIONTRACKER_H_

namespace dare::utils {

/*! \class InitializationTracker
 * @brief little helper to a keep track of initialization
 */
class InitializationTracker {
public:
    /*!
     * @brief constructor
     * @param init identifier if inizialized at construction
     */
    explicit InitializationTracker(bool init = false) : is_initialized(init) {}

    /*!
     * @brief default constructor
     */
    virtual ~InitializationTracker() {}

    /*!
     * @brief default copy constructor
     * @param other 
     */
    InitializationTracker(const InitializationTracker& other) = default;

    /*!
     * @brief default copy assignment operator
     * @param other 
     * @return 
     */
    InitializationTracker& operator=(const InitializationTracker& other) = default;

    /*!
     * @brief provides status
     * @return true if initialized
     */
    bool IsInitialized() const {
        return is_initialized;
    }

protected:
    /*!
     * @brief sets tracker to true
     */
    void Initialize() {
        is_initialized = true;
    }

private:
    bool is_initialized;    //!< identifier if initialized
};

}  // namespace dare::utils

#endif  // UTILITIES_INITIALIZATIONTRACKER_H_
