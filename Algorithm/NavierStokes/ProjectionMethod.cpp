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

#include "ProjectionMethod.h"
#include "Utilities/PropertyInformation.h"

namespace dare {
namespace detail {
ParameterList GetDefaultParameterListPM() {
    ParameterList l;
    l.set("continuity: Newton iterations max", 20,
            "The maximum number of Newton iterations for minimizing the continuity defect");
    l.set("continuity: solver properties", dare::None{},
            "A list of solver properties, depends on the solver backend");
    l.set("continuity: Jacobian", "mutable", "An optimization option for avoiding recomputation of the Jacobian");
    l.set("momentum: Newton iterations max", 1,
          "The default maximum number of Newton iterations for minimizing the defect of the momentum equation (currently not in use)");  // NOLINT
    l.set("momentum[0]: Newton iterations max", 1, "maximum Newton iterations in x-direction");
    l.set("momentum[1]: Newton iterations max", 1, "maximum Newton iterations in y-direction");
    l.set("momentum[2]: Newton iterations max", 1, "maximum Newton iterations in z-direction");
    l.set("momentum: solver properties", dare::None{},
            "A list of default solver properties, depends on the solver backend");
    l.set("momentum[0]: solver properties", dare::None{}, "solver properties for the momentum equation x-direction");
    l.set("momentum[1]: solver properties", dare::None{}, "solver properties for the momentum equation y-direction");
    l.set("momentum[2]: solver properties", dare::None{}, "solver properties for the momentum equation z-direction");
    l.set("momentum: Jacobian", "mutable", "An optimization option for avoiding recomputation of the Jacobian");
    l.set("momentum[0]: Jacobian", "mutable",
            "An optimization option for avoiding recomputation of the Jacobian in x-direction");
    l.set("momentum[1]: Jacobian", "mutable",
            "An optimization option for avoiding recomputation of the Jacobian in y-direction");
    l.set("momentum[2]: Jacobian", "mutable",
            "An optimization option for avoiding recomputation of the Jacobian in z-direction");

    return l;
}
}  // namespace detail
}  // namespace dare
