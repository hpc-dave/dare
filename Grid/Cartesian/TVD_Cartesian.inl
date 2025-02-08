/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder
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

namespace dare {
template <std::size_t Dim, typename SC, typename FluxLimiter>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::TVD(
    const GridRepresentation& grid,
    LO ordinal_internal,
    dare::Vector<Dim, const dare::GridVector<GridType, SC, 1>*> v)
    : ind(grid.MapOrdinalToIndexLocal(grid.MapInternalToLocal(ordinal_internal))), grep(&grid) {
    for (std::size_t id{0}; id < (Dim * 2); id++) {
        const CartesianNeighbor cnb = ToCartesianNeighbor(id + 1);
        const SC value = InterpolateToFace(grid, ind, cnb, *v[id / 2], 0);
        velocity.SetValue(cnb, 0, value);
        upwind[id] = value >= static_cast<SC>(0.);
    }
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::TVD(const GridRepresentation& grid,
                                                      LO ordinal_internal,
                                                      const dare::Vector<Dim, SC>& v)
    : ind(grid.MapOrdinalToIndexLocal(grid.MapInternalToLocal(ordinal_internal))), grep(&grid) {
    for (std::size_t dim{0}; dim < Dim; dim++) {
        const CartesianNeighbor face_low = ToCartesianNeighbor(dim * 2 + 1);
        const CartesianNeighbor face_up = ToCartesianNeighbor(dim * 2 + 2);
        const SC value = v[dim];
        velocity.SetValue(face_low, 0, value);
        velocity.SetValue(face_up, 0, value);
        upwind[dim * 2] = upwind[dim * 2 + 1] = (value >= static_cast<SC>(0.));
    }
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::~TVD() {
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N>
dare::FaceValueStencil<dare::Cartesian<Dim>, SC, N>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    const dare::ExtendedStencil<dare::CenterValueStencil<GridType, SC, N>>& s) const {
    const dare::CenterValueStencil<GridType, SC, N>& s_close = s.first;
    const dare::CenterValueStencil<GridType, SC, N>& s_far = s.second;
    dare::FaceValueStencil<GridType, SC, N> face_values;
    dare::Vector<N, SC> phi_UU, phi_U, phi_D, r_f, flux_lim;
    Index ind_UU, ind_U, ind_D;
    for (std::size_t d{0}; d < Dim; d++) {
        /*
         * at lower face
         * 1) Determine upwind direction
         * 2) Get Values from stencils
         * 3) Compute face gradient and insert into flux limiter
         * 4) Determine face values
         */
        // 1) determine directions
        SC vel_pos{static_cast<SC>(upwind[d * 2])};
        SC vel_neg{static_cast<SC>(!upwind[d * 2])};
        const CartesianNeighbor center = CartesianNeighbor::CENTER;
        CartesianNeighbor face_low = ToCartesianNeighbor(d * 2 + 1);
        CartesianNeighbor face_up = ToCartesianNeighbor(d * 2 + 2);

        // 2) Get values from stencil
        phi_UU = vel_pos * s_far.GetValues(face_low) + vel_neg * s_close.GetValues(face_up);
        phi_U = vel_pos * s_close.GetValues(face_low) + vel_neg * s_close.GetValues(center);
        phi_D = vel_pos * s_close.GetValues(center) + vel_neg * s_close.GetValues(face_low);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Determine face values
        auto phi_face = phi_U + 0.5 * flux_lim * (phi_D - phi_U);
        face_values.SetValues(ToCartesianNeighbor(d * 2 + 1), phi_face);

        /*
         * The same is happening at the upper face, just one cell further
         */
        // 1) determine directions
        vel_pos = static_cast<SC>(upwind[d * 2 + 1]);
        vel_neg = static_cast<SC>(!upwind[d * 2 + 1]);

        // 2) Get values from stencil
        phi_UU = vel_pos * s_close.GetValues(face_low) + vel_neg * s_far.GetValues(face_up);
        phi_U = vel_pos * s_close.GetValues(center) + vel_neg * s_close.GetValues(face_up);
        phi_D = vel_pos * s_close.GetValues(face_up) + vel_neg * s_close.GetValues(center);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Determine face values
        phi_face = phi_U + 0.5 * flux_lim * (phi_D - phi_U);
        face_values.SetValues(ToCartesianNeighbor(d * 2 + 2), phi_face);
    }

    return face_values;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N>
dare::FaceValueStencil<dare::Cartesian<Dim>, SC, N>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    const dare::GridVector<GridType, SC, N>& field) const {
    return Interpolate(ComputeExtendedValueStencil(field, field));
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N>
dare::FaceValueStencil<dare::Cartesian<Dim>, SC, N>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    const dare::GridVector<GridType, SC, N>* field) const {
    return Interpolate(*field);
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template<std::size_t N>
dare::FaceValueStencil<dare::Cartesian<Dim>, SC, N>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    const dare::Vector<N, SC>& values) const {
    dare::FaceValueStencil<dare::Cartesian<Dim>, SC, N> f;
    f.SetValues(values);
    return f;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
dare::FaceValueStencil<dare::Cartesian<Dim>, SC, 1>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    SC value) const {
    dare::FaceValueStencil<dare::Cartesian<Dim>, SC, 1> f;
    f.SetAll(value);
    return f;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
dare::FaceValueStencil<dare::Cartesian<Dim>, SC, 1>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    dare::None value) const {
    return Interpolate(1.);
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N>
dare::FaceMatrixStencil<dare::Cartesian<Dim>, SC, N>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::operator*(
    const dare::GridVector<GridType, SC, N>& field) const {
    return Apply(field);
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
const dare::FaceValueStencil<dare::Cartesian<Dim>, SC, 1>&
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::GetVelocityStencil() const {
    return velocity;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
const typename dare::Cartesian<Dim>::Index&
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::GetIndex() const {
    return ind;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <typename... Args>
auto TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::operator()(const Args&... args) const {
    return Apply(args...);
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <typename... Args>
auto TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::Apply(const Args&... args) const {
    auto tuple_val = std::forward_as_tuple(args...);
    using LastType = std::remove_cvref_t<decltype(std::get<sizeof...(args) - 1>(tuple_val))>;
    static const std::size_t N = dare::detail::free_tvd_cartesian_extract_num_components<LastType>::GetValue();
    dare::ExtendedStencil<dare::CenterValueStencil<GridType, SC, N>>
        phi = ComputeExtendedValueStencil(dare::convert_to_ref(std::get<sizeof...(args) - 1>(tuple_val)), args...);

    // for readability, compiler will optimize out (hopefully :) )
    const dare::CenterValueStencil<GridType, SC, N>& phi_close = phi.first;
    const dare::CenterValueStencil<GridType, SC, N>& phi_far = phi.second;

    dare::FaceMatrixStencil<GridType, SC, N> s;
    dare::Vector<N, SC> ONES, phi_UU, phi_U, phi_D, r_f, flux_lim;
    ONES.SetAllValues(static_cast<SC>(1.));     // a little helper for setting the upwind scheme later

    const CartesianNeighbor center = CartesianNeighbor::CENTER;
    for (std::size_t d{0}; d < Dim; d++) {
        /*
         * at lower face
         * 1) Determine upwind direction
         * 2) Get Values from stencils
         * 3) Compute face gradient and insert into flux limiter
         * 4) Determine face values
         */
        // 1) determine directions
        SC vel_pos{static_cast<SC>(upwind[d * 2])};
        SC vel_neg{static_cast<SC>(!upwind[d * 2])};
        // CartesianNeighbor face = ToCartesianNeighbor(d * 2 + 1);
        CartesianNeighbor face_low = ToCartesianNeighbor(d * 2 + 1);
        CartesianNeighbor face_up = ToCartesianNeighbor(d * 2 + 2);

        // 2) Get values from stencil
        phi_UU = vel_pos * phi_far.GetValues(face_low) + vel_neg * phi_close.GetValues(face_up);
        phi_U = vel_pos * phi_close.GetValues(face_low) + vel_neg * phi_close.GetValues(center);
        phi_D = vel_pos * phi_close.GetValues(center) + vel_neg * phi_close.GetValues(face_low);
        SC vel{velocity.GetValue(face_low, 0)};

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Set upwind scheme
        s.SetValueNeighbor(face_low, vel_pos * vel * ONES);
        s.SetValueCenter(face_low, vel_neg * vel * ONES);

        // 5) add deferred correction, note the negative sign to account for rhs
        auto phi_explicit = -0.5 * flux_lim * (phi_D - phi_U) * vel;
        s.SetRHS(face_low, phi_explicit);

        /*
         * The same is happening at the upper face, just one cell further
         */
        // 1) determine directions
        vel_pos = static_cast<SC>(upwind[d * 2 + 1]);
        vel_neg = static_cast<SC>(!upwind[d * 2 + 1]);

        // 2) Get values from stencil
        phi_UU = vel_pos * phi_close.GetValues(face_low) + vel_neg * phi_far.GetValues(face_up);
        phi_U = vel_pos * phi_close.GetValues(center) + vel_neg * phi_close.GetValues(face_up);
        phi_D = vel_pos * phi_close.GetValues(face_up) + vel_neg * phi_close.GetValues(center);
        vel = velocity.GetValue(face_up, 0);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Set Upwind scheme
        s.SetValueNeighbor(face_up, vel_neg * vel * ONES);
        s.SetValueCenter(face_up, vel_pos * vel * ONES);

        // deferred correction, note the negative sign to account for rhs
        phi_explicit = -0.5 * flux_lim * (phi_D - phi_U) * vel;
        s.SetRHS(face_up, phi_explicit);
    }
    return s;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N, typename... Args>
dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>
TVD<dare::Cartesian<Dim>, SC, FluxLimiter>::ComputeExtendedValueStencil(
    const dare::GridVector<GridType, SC, N>& f,
    const Args&... values) const {
    const std::size_t NUM_VALUES = sizeof...(values);
    static_assert(NUM_VALUES > 0, "No values were provided!");

    dare::ExtendedStencil<dare::CenterValueStencil<GridType, SC, N>> s;
    s.first.SetAll(1.);
    s.second.SetAll(1.);
    detail::free_tvd_cartesian_get_extended_stencil(*grep, ind, &s, values...);
    return s;
}

}  // end namespace dare
