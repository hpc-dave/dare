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

#ifndef EQUATIONS_DDT_H_
#define EQUATIONS_DDT_H_

#include <tuple>
#include <type_traits>

#include "Utilities/Vector.h"
#include "Utilities/Errors.h"
#include "Data/Stencil.h"
#include "Data/Field.h"
#include "TimeDiscretizationSchemes.h"

namespace dare::Matrix {

namespace timeschemes {
struct IMPLICIT {
    static const std::size_t NUM_TIMESTEPS{1};

    template <typename SC>
    static constexpr dare::utils::Vector<NUM_TIMESTEPS + 1, SC> GetWeights() {
        return dare::utils::Vector<NUM_TIMESTEPS + 1, SC>(static_cast<SC>(1.), static_cast<SC>(1.));
    }
};
}  // namespace timeschemes

/*!
 * \brief computes time derivative
 * @tparam Grid type of grid
 * @tparam TimeDiscretization discretization scheme for the time derivative
 */
template<typename Grid, typename TimeDiscretization = dare::Matrix::timeschemes::IMPLICIT>
class DDT {
public:
    using LO = typename Grid::LocalOrdinalType;
    using SC = typename Grid::ScalarType;
    using Index = typename Grid::Index;
    using GridRepresentation = typename Grid::Representation;
    static const std::size_t NUM_TIMESTEPS = TimeDiscretization::NUM_TIMESTEPS;
    static const std::size_t NUM_TFIELDS{NUM_TIMESTEPS + 1};

    DDT(const GridRepresentation& grep,
        LO local_ordinal,
        SC dt);

    virtual ~DDT();

    template <typename... Args>
    auto operator()(const Args&... args);

private:
    template <int I, std::size_t NUM_COMPONENTS, typename... Args>
    void Iterate(dare::utils::Vector<NUM_COMPONENTS, dare::utils::Vector<NUM_TFIELDS, SC>>& v,  // NOLINT
                 const std::tuple<const Args&...>& args);

    SC dt;              //!< discrete timestep
    LO ordinal;         //!< internal ordinal
    LO local_ordinal;   //!< local ordinal
    SC volume;          //!< discrete volume
};

}  // end namespace dare::Matrix

#include "DDT.inl"

#endif  // EQUATIONS_DDT_H_
