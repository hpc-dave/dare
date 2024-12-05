/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder

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
#ifndef MATRIXSYSTEM_OPERATORS_CARTESIAN_H_
#define MATRIXSYSTEM_OPERATORS_CARTESIAN_H_

namespace dare::Matrix {

template <std::size_t Dim, typename TimeDiscretization>
Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::Divergence(
    const GridRepresentation& grid, LO ordinal_internal)
    : A(grid.GetFaceArea()),
      ind(grid.MapInternalToLocal(grid.MapOrdinalToIndexLocalInternal(ordinal_internal))),
      grep(&grid) {
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename... Args>
auto Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::operator()(const Args&... args) {
    const int LAST_POS{sizeof...(args) - 1};
    static_assert(LAST_POS >= 0, "minimally one argument (the relevant field) is required!");
    auto tuple_val = std::forward_as_tuple(args...);
    using LastType = std::remove_cv_t<std::remove_reference_t<decltype(std::get<LAST_POS>(tuple_val))>>;
    const std::size_t NUM_COMPONENTS{LastType::NUM_COMPONENTS};
    using SC = typename LastType::ScalarType;

#ifndef DARE_NDEBUG
    if constexpr(dare::is_field_v<LastType>) {
        std::size_t num_last_timesteps = std::get<LAST_POS>(std::forward_as_tuple(args...)).GetNumberTimesteps();
        if (num_last_timesteps < NUM_TFIELDS) {
            ERROR << "The provided field does not store enough time-steps! Minimally required for "
                << typeid(TimeDiscretization).name() << " " << (NUM_TFIELDS > 1 ? "are" : "is") << " " << NUM_TFIELDS
                << ", but found " << num_last_timesteps << ERROR_CLOSE;
        }
    }
#endif

    if constexpr (dare::is_field_v<LastType> || dare::is_face_matrix_stencil_v<LastType>) {
        TFaceMatrixStencil<SC, NUM_COMPONENTS> f;
        f[0] = GetFaceMatrixStencil(std::get<LAST_POS>(tuple_val));

        for (std::size_t t{1}; t < NUM_TFIELDS; t++)
            f[t] = f[0];

        // mutiply with all remaining tuple elements
        MultiplyAll(&f, args...);

        // apply divergence
        return ApplyDivergence(f);
    } if constexpr (dare::is_face_value_stencil_v<LastType> || dare::is_gridvector_v<LastType>) {
        dare::Data::FaceValueStencil<GridType, SC, NUM_COMPONENTS> f;
        if constexpr (dare::is_gridvector_v<LastType>) {
            f = PopulateFaceValueFromField(std::get<LAST_POS>(tuple_val));
        } else {
            f = std::get<LAST_POS>(tuple_val);
        }

        // mutiply with all remaining tuple elements
        MultiplyAll(&f, args...);

        // Apply divergence
        return ApplyDivergence(f);
    } else {
        // bit of a workaround, since some compilers have issues if there is a 'false'
        // provided directly
        // static_assert(dare::always_false<>, "No implementation provided for this type");
    }
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
void Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::Multiply(
    SC value, TFaceMatrixStencil<SC, N>* s) const {
    for (auto& i : *s)
        i *= value;
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
void Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::Multiply(
    const dare::Data::FaceValueStencil<GridType, SC, N>& f,
    TFaceMatrixStencil<SC, N>* s) const {
    for (auto& i : *s)
        i *= f;
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
void Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::Multiply(
    const TFaceValueStencil<SC, N>& f,
    TFaceMatrixStencil<SC, N>* s) const {
#ifndef DARE_NDEBUG
    if (f.size() < NUM_TFIELDS) {
        ERROR << "The number of timesteps provided in f is not enough ("
              << f.size() << " < " << NUM_TFIELDS << ")" << ERROR_CLOSE;
    }
#endif
    for (std::size_t n{0}; n < NUM_TFIELDS; n++)
        (*s)[n] *= f[n];
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
void Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::Multiply(
    const dare::Data::GridVector<GridType, SC, N>& f, TFaceMatrixStencil<SC, N>* s) const {
    Multiply(PopulateFaceValueFromField(f), s);
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
void Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::Multiply(
    const dare::Data::Field<GridType, SC, N>& f, TFaceMatrixStencil<SC, N>* s) const {
    for (std::size_t n{0}; n < NUM_TIMESTEPS; n++) {
        s->At(n) *= PopulateFaceValueFromField(f.GetDataVector(n));
    }
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
dare::Data::CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::ApplyDivergence(
    const TFaceMatrixStencil<SC, N>& s) const {
    dare::Data::CenterMatrixStencil<GridType, SC, N> s_c;
    // at timestep 0, we deal with the implicit components
    for (std::size_t n{0}; n < N; n++) {
        // Divergence in X
        SC coef_c = s[0].GetValueCenter(Positions::EAST, n);
        SC coef_f = s[0].GetValueNeighbor(Positions::EAST, n);
        s_c.GetValue(Positions::EAST, n) = A[0] * coef_f;
        s_c.GetValue(Positions::CENTER, n) = A[0] * coef_c;

        coef_c = s[0].GetValueCenter(Positions::WEST, n);
        coef_f = s[0].GetValueNeighbor(Positions::WEST, n);
        s_c.GetValue(Positions::WEST, n) = -A[0] * coef_f;
        s_c.GetValue(Positions::CENTER, n) -= A[0] * coef_c;

        // Divergence in Y
        if constexpr (Dim > 1) {
            SC coef_c = s[0].GetValueCenter(Positions::NORTH, n);
            SC coef_f = s[0].GetValueNeighbor(Positions::NORTH, n);
            s_c.GetValue(Positions::NORTH, n) = A[1] * coef_f;
            s_c.GetValue(Positions::CENTER, n) += A[1] * coef_c;

            coef_c = s[0].GetValueCenter(Positions::SOUTH, n);
            coef_f = s[0].GetValueNeighbor(Positions::SOUTH, n);
            s_c.GetValue(Positions::SOUTH, n) = -A[1] * coef_f;
            s_c.GetValue(Positions::CENTER, n) -= A[1] * coef_c;
        }

        // Divergence in Z
        if constexpr (Dim > 2) {
            SC coef_c = s[0].GetValueCenter(Positions::TOP, n);
            SC coef_f = s[0].GetValueNeighbor(Positions::TOP, n);
            s_c.GetValue(Positions::TOP, n) = A[2] * coef_f;
            s_c.GetValue(Positions::CENTER, n) += A[2] * coef_c;

            coef_c = s[0].GetValueCenter(Positions::BOTTOM, n);
            coef_f = s[0].GetValueNeighbor(Positions::BOTTOM, n);
            s_c.GetValue(Positions::BOTTOM, n) = -A[2] * coef_f;
            s_c.GetValue(Positions::CENTER, n) -= A[2] * coef_c;
        }
    }
    const SC coef_td = TimeDiscretization::template GetWeights<SC>()[0];
    s_c *= coef_td;  // Apply scheme weight

    // Add explicit components (e.g. from TVD Schemes)
    for (std::size_t n{0}; n < N; n++) {
        s_c.GetRHS(n) = A[0] * (s[0].GetRHS(Positions::EAST, n) - s[0].GetRHS(Positions::WEST, n));
        if constexpr(Dim > 1)
            s_c.GetRHS(n) += A[1] * (s[0].GetRHS(Positions::NORTH, n) - s[0].GetRHS(Positions::SOUTH, n));
        if constexpr(Dim > 2)
            s_c.GetRHS(n) += A[2] * (s[0].GetRHS(Positions::TOP, n) - s[0].GetRHS(Positions::BOTTOM, n));
    }

    // the remaining components are explicit
    for (std::size_t t{1}; t < NUM_TFIELDS; t++) {
        for (std::size_t n{0}; n < N; n++) {
            const SC coef_td = TimeDiscretization::template GetWeights<SC>()[t];
            s_c.GetRHS(n) += coef_td * A[0] * (s[t].GetRHS(Positions::EAST, n) - s[t].GetRHS(Positions::WEST, n));
            if constexpr (Dim > 1)
                s_c.GetRHS(n) += coef_td * A[1] * (s[t].GetRHS(Positions::NORTH, n) - s[t].GetRHS(Positions::SOUTH, n));
            if constexpr (Dim > 2)
                s_c.GetRHS(n) += coef_td * A[2] * (s[t].GetRHS(Positions::TOP, n) - s[t].GetRHS(Positions::BOTTOM, n));
            }
    }
    return s_c;
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
dare::utils::Vector<N, SC>
Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::ApplyDivergence(
    const dare::Data::FaceValueStencil<GridType, SC, N>& s) const {
    // Note that here we don't care about the time discretization, since
    // we have no temporal information
    dare::utils::Vector<N, SC> div_v;
    for (std::size_t n{0}; n < N; n++) {
        // Divergence in X
        div_v[n] = A[0] * s.GetValue(Positions::EAST, n);
        div_v[n] -= A[0] * s.GetValue(Positions::WEST, n);

        // Divergence in Y
        if constexpr (Dim > 1) {
            div_v[n] += A[1] * s.GetValue(Positions::NORTH, n);
            div_v[n] -= A[1] * s.GetValue(Positions::SOUTH, n);
        }

        // Divergence in Z
        if constexpr (Dim > 2) {
            div_v[n] += A[2] * s.GetValue(Positions::TOP, n);
            div_v[n] -= A[2] * s.GetValue(Positions::BOTTOM, n);
        }
    }
    return div_v;
}

// template <std::size_t Dim, typename TimeDiscretization>
// template <typename SC, std::size_t N>
// dare::Data::CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
// Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::Apply(
//     const dare::Data::FaceMatrixStencil<GridType, SC, N>& s) const {
//     dare::Data::CenterMatrixStencil<GridType, SC, N> s_c;
//     for (std::size_t n{0}; n < N; n++) {
//         // Divergence in X
//         SC coef_c = s.GetValueCenter(Positions::EAST, n);
//         SC coef_f = s.GetValueNeighbor(Positions::EAST, n);
//         s_c.GetValue(Positions::EAST, n) = A[0] * coef_f;
//         s_c.GetValue(Positions::CENTER, n) = A[0] * coef_c;
//         s_c.GetRHS(n) = A[0] * s.GetRHS(Positions::EAST, n);

//         coef_c = s.GetValueCenter(Positions::WEST, n);
//         coef_f = s.GetValueNeighbor(Positions::WEST, n);
//         s_c.GetValue(Positions::WEST, n) = -A[0] * coef_f;
//         s_c.GetValue(Positions::CENTER, n) -= A[0] * coef_c;
//         s_c.GetRHS(n) -= A[0] * s.GetRHS(Positions::WEST, n);

//         // Divergence in Y
//         if constexpr (Dim > 1) {
//             SC coef_c = s.GetValueCenter(Positions::NORTH, n);
//             SC coef_f = s.GetValueNeighbor(Positions::NORTH, n);
//             s_c.GetValue(Positions::NORTH, n) = A[1] * coef_f;
//             s_c.GetValue(Positions::CENTER, n) += A[1] * coef_c;
//             s_c.GetRHS(n) += A[1] * s.GetRHS(Positions::NORTH, n);

//             coef_c = s.GetValueCenter(Positions::SOUTH, n);
//             coef_f = s.GetValueNeighbor(Positions::SOUTH, n);
//             s_c.GetValue(Positions::SOUTH, n) = -A[1] * coef_f;
//             s_c.GetValue(Positions::CENTER, n) -= A[1] * coef_c;
//             s_c.GetRHS(n) -= A[1] * s.GetRHS(Positions::SOUTH, n);
//         }

//         // Divergence in Z
//         if constexpr (Dim > 2) {
//             SC coef_c = s.GetValueCenter(Positions::TOP, n);
//             SC coef_f = s.GetValueNeighbor(Positions::TOP, n);
//             s_c.GetValue(Positions::TOP, n) = A[2] * coef_f;
//             s_c.GetValue(Positions::CENTER, n) += A[2] * coef_c;
//             s_c.GetRHS(n) += A[2] * s.GetRHS(Positions::TOP, n);

//             coef_c = s.GetValueCenter(Positions::BOTTOM, n);
//             coef_f = s.GetValueNeighbor(Positions::BOTTOM, n);
//             s_c.GetValue(Positions::BOTTOM, n) = -A[2] * coef_f;
//             s_c.GetValue(Positions::CENTER, n) -= A[2] * coef_c;
//             s_c.GetRHS(n) -= A[2] * s.GetRHS(Positions::BOTTOM, n);
//         }
//     }
//     return s_c;
// }

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::PopulateFaceValueFromField(
    const dare::Data::GridVector<GridType, SC, N>& f) const {
    static_assert(dare::always_false<>, "need to interpolated those values!");
    using Pos = typename dare::Grid::Cartesian<Dim>::NeighborID;
    Index ind_nb{ind};
    dare::Data::FaceValueStencil<GridType, SC, N> s;

    ind_nb.i()--;
    auto v = dare::math::InterpolateToFace(*grep, ind, Pos::WEST, f);
    for (std::size_t n{0}; n < N; n++)
        s.SetValue(Pos::WEST, v[n]);
    ind_nb.i() += 2;
    for (std::size_t n{0}; n < N; n++)
        s.SetValue(Pos::EAST, f.At(ind_nb, n));
    if constexpr (Dim > 1) {
        ind_nb.i()--;
        ind_nb.j()--;
        for (std::size_t n{0}; n < N; n++)
            s.SetValue(Pos::SOUTH, f.At(ind_nb, n));
        ind_nb.j() += 2;
        for (std::size_t n{0}; n < N; n++)
            s.SetValue(Pos::NORTH, f.At(ind_nb, n));
    }
    if constexpr (Dim > 2) {
        ind_nb.j()--;
        ind_nb.k()--;
        for (std::size_t n{0}; n < N; n++)
            s.SetValue(Pos::BOTTOM, f.At(ind_nb, n));
        ind_nb.k() += 2;
        for (std::size_t n{0}; n < N; n++)
            s.SetValue(Pos::TOP, f.At(ind_nb, n));
    }
    return s;
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
typename Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::TFaceValueStencil<SC, N>
Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::PopulateFaceValueFromField(
    const dare::Data::Field<GridType, SC, N>& f) const {
    static_assert(dare::always_false<>, "need to include temporal information");
    return PopulateFaceValueFromField(f.GetDataVector(0));
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename Stencil, typename Arg, typename... Args>
void Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::MultiplyAll(
    Stencil* f, const Arg& arg, const Args&... args) {
    if constexpr (sizeof...(args) > 0) {
        Multiply(arg, f);
        MultiplyAll(f, args...);
    }
}


template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
dare::Data::FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::GetFaceMatrixStencil(
    const dare::Data::Field<GridType, SC, N>& field) const {
    dare::Data::FaceMatrixStencil<GridType, SC, N> f;
    f.SetValues(Positions::WEST, 0, -1., 1.);
    f.SetValues(Positions::EAST, 0, 1., -1.);
    if constexpr (Dim > 1) {
        f.SetValues(Positions::SOUTH, 0, -1., 1.);
        f.SetValues(Positions::NORTH, 0, 1., -1.);
    }
    if constexpr (Dim > 2) {
        f.SetValues(Positions::BOTTOM, 0, -1., 1.);
        f.SetValues(Positions::TOP, 0, 1., -1.);
    }

    for (std::size_t n{1}; n < N; n++) {
        f.GetDataCenter()[n] = f.GetDataCenter()[n - 1];
        f.GetDataNeighbor()[n] = f.GetDataNeighbor()[n - 1];
    }
    return f;
}

template <std::size_t Dim, typename TimeDiscretization>
template <typename SC, std::size_t N>
dare::Data::FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
Divergence<dare::Grid::Cartesian<Dim>, TimeDiscretization>::GetFaceMatrixStencil(
    const dare::Data::FaceMatrixStencil<GridType, SC, N>& s) const {
    return s;
}

}  // end namespace dare::Matrix

#endif  // MATRIXSYSTEM_OPERATORS_CARTESIAN_H_
