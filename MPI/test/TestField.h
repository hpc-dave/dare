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


#ifndef MPI_TEST_TESTFIELD_H_
#define MPI_TEST_TESTFIELD_H_

#include <vector>

namespace dare::test {
class TestField {
public:
    void ResizeByGridSize(int grid_size) {
        data.resize(grid_size * GetNumEquations());
    }
    double* GetData() {
        return data.data();
    }
    std::size_t GetNumEquations() const {
        return num_eq;
    }

    double& at(std::size_t pos) {
        return data[pos];
    }
    double at(std::size_t pos) const {
        return data[pos];
    }

private:
    const std::size_t num_eq = 5;
    std::vector<double> data;
};
}  // namespace dare::test

#endif  // MPI_TEST_TESTFIELD_H_
