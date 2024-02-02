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

namespace dare::Matrix {

template <typename Grid, typename TimeDiscretization>
DDT<Grid, TimeDiscretization>::DDT(const GridRepresentation& grep,
         LO lo,
         SC _dt)
         : dt(_dt), ordinal(lo), volume(grep.GetCellVolume(lo)) {
}

template <typename Grid, typename TimeDiscretization>
DDT<Grid, TimeDiscretization>::~DDT() {
}

template <typename Grid, typename TimeDiscretization>
template <typename... Args>
auto DDT<Grid, TimeDiscretization>::operator()(const Args&... args) {
    const int LAST_POS{sizeof...(args) - 1};
    static_assert(LAST_POS >= 0, "minimally one argument (the relevant field) is required!");
    auto tuple_val = std::forward_as_tuple(args...);
    using LastType = std::remove_reference_t<decltype(std::get<LAST_POS>(tuple_val))>;
    static_assert(dare::is_field_v<LastType>, "The last element needs to be a field!");
    const std::size_t NUM_COMPONENTS{LastType::NUM_COMPONENTS};
    using TransientTerms = dare::utils::Vector<NUM_TFIELDS, SC>;
    using ComponentVector = dare::utils::Vector<NUM_COMPONENTS, TransientTerms>;

#ifndef DARE_NDEBUG
    std::size_t num_last_timesteps = std::get<LAST_POS>(std::forward_as_tuple(args...)).GetNumberTimesteps();
    if (num_last_timesteps < NUM_TFIELDS) {
        ERROR << "The provided field does not store enough time-steps! Minimally required for "
        << typeid(TimeDiscretization).name() << " " << (NUM_TFIELDS > 1?"are": "is") << " " << NUM_TFIELDS
        << ", but found " << num_last_timesteps << ERROR_CLOSE;
    }
#endif
    ComponentVector transient_terms;
    // set initial values
    transient_terms[0] = volume * TimeDiscretization::template GetWeights<SC>() / dt;
    for (std::size_t n{1}; n < NUM_COMPONENTS; n++) {
        transient_terms[n] = transient_terms[n - 1];
    }

    // multiply with the remaining arguments
    Iterate<LAST_POS - 1>(transient_terms, std::forward_as_tuple(args...));

    // add the information to the matrix stencil
    using Stencil = Data::CenterMatrixStencil<Grid, SC, NUM_COMPONENTS>;
    Stencil stencil;
    for (std::size_t n{0}; n < NUM_COMPONENTS; n++) {
        stencil.Center(n) = transient_terms[n][0];
        for (std::size_t t{1}; t < NUM_TFIELDS; t++) {
            stencil.GetRHS(n) -= transient_terms[n][t];
        }
    }
    return stencil;
}

template <typename Grid, typename TimeDiscretization>
template <int I, std::size_t NUM_COMPONENTS, typename... Args>
void DDT<Grid, TimeDiscretization>::Iterate(
    dare::utils::Vector<NUM_COMPONENTS, dare::utils::Vector<NUM_TFIELDS, SC>>& transient_terms, // NOLINT
    const std::tuple<const Args&...>& args) {
    if constexpr (I == -1) {
        // end of iteration, do nothing
        // Note, this is a convenient choice here, so we don't need to make a decision in the initial call of
        // this recursive function
    } else {
        using Type = std::remove_reference_t<decltype(std::get<I>(args))>;
        static_assert(dare::is_field_v<Type>, "the provided arguments need to be of type Field!");

        auto field = std::get<I>(args);
        // apply values
        for (std::size_t n{0}; n < NUM_COMPONENTS; n++) {
            for (std::size_t t{0}; t < NUM_TFIELDS; t++) {
                transient_terms[n][t] *= field.GetDataVector(t).At(ordinal, n);
            }
        }

        // iterate one element further
        Iterate<I - 1>(transient_terms, args);
    }
}
}  // end namespace dare::Matrix
