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

#ifndef ALGORITHM_NAVIERSTOKES_PMTERMINALOUTPUTDEFAULT_H_
#define ALGORITHM_NAVIERSTOKES_PMTERMINALOUTPUTDEFAULT_H_

namespace dare {
class PMTerminalOutputDefault{
public:
    explicit PMTerminalOutputDefault(int dim, int acc = 2, int width = 15);

    virtual ~PMTerminalOutputDefault();

    template<typename PM>
    void PrintHeader(const PM& pm);
    void PrintHeader(uint64_t tstep, double time);

    template <typename PM>
    void PrintMomentum(int dim, int iter, bool success, const PM& pm);
    void PrintMomentum(int dim, int iter, bool success);

    template <typename PM>
    void PrintInitialDefect(const PM& pm);
    void PrintInitialDefectInternal(double defect);

    template<typename PM>
    void PrintContinuity(int loop, int iter, bool success, const PM& pm);
    void PrintContinuity(int loop, int iter, bool success, double defect, bool is_last);

private:
    int dimension;
    int accuracy;
    int print_width;
};

template <typename PM>
void PMTerminalOutputDefault::PrintHeader(const PM& pm) {
    PrintHeader(pm.GetTimeStepCounter(), pm.GetTime());
}


template <typename PM>
void PMTerminalOutputDefault::PrintMomentum(int dim, int iter, bool success, const PM& pm) {
    PrintMomentum(dim, iter, success);
}

template<typename PM>
void PMTerminalOutputDefault::PrintInitialDefect(const PM& pm) {
    PrintInitialDefectInternal(pm.GetMaxContinuityDefect());
}

template <typename PM>
void PMTerminalOutputDefault::PrintContinuity(int loop, int iter, bool success, const PM& pm) {
    PrintContinuity(loop, iter, success, pm.GetMaxContinuityDefect(), iter == pm.GetMaxLoopIterations() - 1);
}

}  // namespace dare
#endif  // ALGORITHM_NAVIERSTOKES_PMTERMINALOUTPUTDEFAULT_H_
