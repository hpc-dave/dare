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
template <std::size_t Dim, typename SC, typename FluxLimiter>
TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter>::TVD(
    const GridRepresentation& grid,
    LO ordinal_internal,
    dare::utils::Vector<Dim, const dare::Data::GridVector<GridType, SC, 1>*> v)
    : ind(grid.MapOrdinalToIndexLocal(grid.MapInternalToLocal(ordinal_internal))) {
    for (std::size_t id{0}; id < (Dim * 2); id++) {
        const Grid::CartesianNeighbor cnb = Grid::ToCartesianNeighbor(id + 1);
        const SC value = math::InterpolateToFace(grid, ind, cnb, *v[id / 2]);
        velocity.SetValue(cnb, 0, value);
        upwind[id] = value >= static_cast<SC>(0.);
    }
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter>::TVD(const GridRepresentation& grid,
                                                      LO ordinal_internal,
                                                      const dare::utils::Vector<Dim, SC>& v)
    : ind(grid.MapOrdinalToIndexLocal(grid.MapInternalToLocal(ordinal_internal))) {
    for (std::size_t dim{0}; dim < Dim; dim++) {
        const Grid::CartesianNeighbor cnb_low = Grid::ToCartesianNeighbor(dim * 2 + 1);
        const Grid::CartesianNeighbor cnb_up = Grid::ToCartesianNeighbor(dim * 2 + 2);
        const SC value = v[dim];
        velocity.SetValue(cnb_low, 0, value);
        velocity.SetValue(cnb_up, 0, value);
        upwind[dim * 2] = upwind[dim * 2 + 1] = (value >= static_cast<SC>(0.));
    }
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter>::~TVD() {
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N>
dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    const dare::Data::CenterValueStencil<GridType, SC, N>& s_close,
    const dare::Data::CenterValueStencil<GridType, SC, N>& s_far) const {
    dare::Data::FaceValueStencil<GridType, SC, N> face_values;
    dare::utils::Vector<N, SC> phi_UU, phi_U, phi_D, r_f, flux_lim;
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
        SC vel_pos{upwind[d * 2]};
        SC vel_neg{!upwind[d * 2]};
        const Grid::CartesianNeighbor center = Grid::CartesianNeighbor::CENTER;
        Grid::CartesianNeighbor cnb_low = Grid::ToCartesianNeighbor(d * 2 + 1);
        Grid::CartesianNeighbor cnb_up = Grid::ToCartesianNeighbor(d * 2 + 1);

        // 2) Get values from stencil
        phi_UU = vel_pos * s_far.GetValues(cnb_low) + vel_neg * s_close.GetValues(cnb_up);
        phi_U = vel_pos * s_close.GetValues(cnb_low) + vel_neg * s_close.GetValues(center);
        phi_D = vel_pos * s_close.GetValues(center) + vel_neg * s_close.GetValues(cnb_low);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Determine face values
        auto phi_face = phi_U + 0.5 * flux_lim * (phi_D - phi_U);
        face_values.SetValues(Grid::ToCartesianNeighbor(d * 2 + 1), phi_face);

        /*
         * The same is happening at the upper face, just one cell further
         */
        // 1) determine directions
        vel_pos = upwind[d * 2 + 1];
        vel_neg = !upwind[d * 2 + 1];

        // 2) Get values from stencil
        phi_UU = vel_pos * s_close.GetValues(cnb_low) + vel_neg * s_far.GetValues(cnb_up);
        phi_U = vel_pos * s_close.GetValues(center) + vel_neg * s_close.GetValues(cnb_up);
        phi_D = vel_pos * s_close.GetValues(cnb_up) + vel_neg * s_close.GetValues(center);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Determine face values
        phi_face = phi_U + 0.5 * flux_lim * (phi_D - phi_U);
        face_values.SetValues(Grid::ToCartesianNeighbor(d * 2 + 2), phi_face);
    }

    return face_values;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N>
dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter>::Interpolate(
    const dare::Data::GridVector<GridType, SC, N>& field) const {
    dare::Data::FaceValueStencil<GridType, SC, N> face_values;
    dare::utils::Vector<N, SC> phi_UU, phi_U, phi_D, r_f, flux_lim;
    Index ind_UU, ind_U, ind_D;
    for (std::size_t d{0}; d < Dim; d++) {
        /*
         * at lower face
         * 1) Get indices depending on upwind
         * 2) Query field values
         * 3) Compute face gradient and insert into flux limiter
         * 4) Determine face values
         */
        // 1) get indices
        ind_UU = ind_U = ind_D = ind;
        LO vel_pos{upwind[d * 2]};
        LO vel_neg{!upwind[d * 2]};
        // LO upwind_dir{!upwind[d * 2] - upwind[2 * d]};
        ind_UU[d] += vel_neg - 2 * vel_pos;
        ind_U[d] -= vel_pos;
        ind_D[d] -= vel_neg;

        // 2) Query field values
        phi_UU = field.GetValues(ind_UU);
        phi_U = field.GetValues(ind_U);
        phi_D = field.GetValues(ind_D);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Determine face values
        auto phi_face = phi_U + 0.5 * flux_lim * (phi_D - phi_U);
        face_values.SetValues(Grid::ToCartesianNeighbor(d * 2 + 1), phi_face);

        /*
         * The same is happening at the upper face, just one cell further
         */
        // 1) reset indices
        ind_UU = ind_U = ind_D = ind;
        vel_pos = upwind[d * 2 + 1];
        vel_neg = !upwind[d * 2 + 1];
        ind_UU[d] += 2 * vel_neg - vel_pos;
        ind_U[d] += vel_neg;
        ind_D[d] += vel_pos;

        // 2) Query field values
        phi_UU = field.GetValues(ind_UU);
        phi_U = field.GetValues(ind_U);
        phi_D = field.GetValues(ind_D);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Determine face values
        phi_face = phi_U + 0.5 * flux_lim * (phi_D - phi_U);
        face_values.SetValues(Grid::ToCartesianNeighbor(d * 2 + 2), phi_face);
    }
    return face_values;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
template <std::size_t N>
dare::Data::FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N>
TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter>::operator*(
    const dare::Data::GridVector<GridType, SC, N>& field) const {
    dare::Data::FaceMatrixStencil<GridType, SC, N> s;
    dare::utils::Vector<N, SC> ONES;
    ONES.SetAllValues(static_cast<SC>(1.));
    dare::utils::Vector<N, SC> phi_UU, phi_U, phi_D, r_f, flux_lim;
    Index ind_UU, ind_U, ind_D;
    for (std::size_t d{0}; d < Dim; d++) {
        /*
         * at lower face
         * 1) Get indices depending on upwind
         * 2) Query field values
         * 3) Compute face gradient and insert into flux limiter
         * 4) Set face values and deferred correction
         */
        // 1) get indices
        ind_UU = ind_U = ind_D = ind;
        LO vel_pos{upwind[d * 2]};
        LO vel_neg{!upwind[d * 2]};
        Grid::CartesianNeighbor face = Grid::ToCartesianNeighbor(d * 2 + 1);
        ind_UU[d] += vel_neg - 2 * vel_pos;
        ind_U[d] -= vel_pos;
        ind_D[d] -= vel_neg;

        // 2) Query field values
        phi_UU = field.GetValues(ind_UU);
        phi_U = field.GetValues(ind_U);
        phi_D = field.GetValues(ind_D);
        SC vel{velocity.GetValue(face, 0)};

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Set face values
        s.SetValueNeighbor(face, vel_pos * vel * ONES);
        s.SetValueCenter(face, vel_neg * vel * ONES);

        // deferred correction, note the negative sign to account for rhs
        auto phi_explicit = -0.5 * flux_lim * (phi_D - phi_U) * vel;
        s.SetRHS(face, phi_explicit);

        /*
         * The same is happening at the upper face, just one cell further
         */
        // 1) reset indices
        ind_UU = ind_U = ind_D = ind;
        vel_pos = upwind[d * 2 + 1];
        vel_neg = !upwind[d * 2 + 1];
        ind_UU[d] += 2 * vel_neg - vel_pos;
        ind_U[d] += vel_neg;
        ind_D[d] += vel_pos;
        face = Grid::ToCartesianNeighbor(d * 2 + 2);

        // 2) Query field values
        phi_UU = field.GetValues(ind_UU);
        phi_U = field.GetValues(ind_U);
        phi_D = field.GetValues(ind_D);
        vel = velocity.GetValue(face, 0);

        // 3) Compute face gradient and flux-limiter
        r_f = (phi_U - phi_UU) / (phi_D - phi_U);
        flux_lim = FluxLimiter::GetValue(r_f);

        // 4) Set face values
        s.SetValueNeighbor(face, vel_neg * vel * ONES);
        s.SetValueCenter(face, vel_pos * vel * ONES);

        // deferred correction, note the negative sign to account for rhs
        phi_explicit = -0.5 * flux_lim * (phi_D - phi_U) * vel;
        s.SetRHS(face, phi_explicit);
    }
    return s;
}

template <std::size_t Dim, typename SC, typename FluxLimiter>
const dare::Data::FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, 1>&
TVD<dare::Grid::Cartesian<Dim>, SC, FluxLimiter>::GetVelocities() const {
    return velocity;
}

}  // end namespace dare::Matrix
