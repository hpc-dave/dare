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
#ifndef GRID_CARTESIAN_TVD_CARTESIAN_H_
#define GRID_CARTESIAN_TVD_CARTESIAN_H_

#include <utility>

#include "Utilities/CompileTimeFunctions.h"
#include "Data/GridVector.h"
#include "Equations/Operators.h"
#include "Grid/Cartesian/Interpolation_Cartesian.h"
#include "Grid/Cartesian/MatrixBlock_Cartesian.h"
#include "Grid/Cartesian/Stencils_Cartesian.h"

namespace dare {

template <std::size_t Dim, typename SC, typename FluxLimiter>
class TVD<dare::Cartesian<Dim>, SC, FluxLimiter> {
public:
    static const std::size_t NUM_FACES{2 * Dim};
    using GridType = dare::Cartesian<Dim>;
    using GridRepresentation = typename GridType::Representation;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using Index = typename GridType::Index;
    using Options = typename GridType::Options;
    using Positions = typename GridType::NeighborID;

    /*!
     * @brief initialization with references to velocity fields
     * @param grid representation of the target grid
     * @param ordinal_internal internal ordinal
     * @param v vector with velocity fields
     */
    TVD(const GridRepresentation& grid,
        LO ordinal_internal,
        dare::Vector<Dim, const dare::GridVector<GridType, SC, 1>*> v);

    /*!
     * @brief initialization with constant velocities
     * @param grid
     * @param ordinal_internal
     * @param v
     */
    TVD(const GridRepresentation& grid,
        LO ordinal_internal,
        const dare::Vector<Dim, SC>& v);

    /*!
     * @brief destructor
     */
    ~TVD();

    /*!
     * @brief When used for interpolation
     * @tparam N number of components
     * @param s_close CenterMatrixStencil for close neighbors
     * @param s_far CenterValueStencil for far neighbors
     * Here, the face values are computed according to the TVD scheme,
     * excluding the velocity component
     */
    template <std::size_t N>
    [[nodiscard]] dare::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::ExtendedStencil<dare::CenterValueStencil<GridType, SC, N>>& s) const;

    /*!
     * @brief When used for interpolation from a field
     * @tparam N number of components
     * @param field field
     * Here, the face values are computed according to the TVD scheme,
     * excluding the velocity component
     */
    template <std::size_t N>
    [[nodiscard]] dare::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::GridVector<GridType, SC, N>& field) const;

    /*!
     * @brief When used for interpolation from a field
     * @tparam N number of components
     * @param field field
     * Here, the face values are computed according to the TVD scheme,
     * excluding the velocity component
     */
    template <std::size_t N>
    [[nodiscard]] dare::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::GridVector<GridType, SC, N>* field) const;

    /*!
     * @brief When used for interpolation from a field
     * @tparam N number of components
     * @param value constant values
     * Here, the face values are set to the values depending on the component
     */
    template <std::size_t N>
    [[nodiscard]] dare::FaceValueStencil<GridType, SC, N> Interpolate(
        const dare::Vector<N, SC>& values) const;

    /*!
     * @brief When used for interpolation from a field
     * @tparam N number of components
     * @param value constant value
     * Here, the face values are set to the value
     */
    [[nodiscard]] dare::FaceValueStencil<GridType, SC, 1> Interpolate(
                    SC value) const;

    /*!
     * @brief When used for interpolation from a field
     * @tparam N number of components
     * @param v none type
     * Here, the face values are set to 1
     */
    [[nodiscard]] dare::FaceValueStencil<GridType, SC, 1> Interpolate(
        dare::None v) const;

    /*!
     * \brief when used for matrix assembly
     * @tparam N number of components
     * @param field reference to relevant field
     * Here, the flux is computed by the TVD scheme, including the velocity component.
     * \f[
     *    J = u \cdot \phi
     * \f]
     * An example code can look like this:
     * @code{.cpp}
     *   GridVector<...> phi;
     *   TVD<...> u(...);
     *   FaceValueStencil<...> rho;
     *
     *   FaceMatrixStencil<...> J = rho * u * phi
     *
     * @endcode
     */
    template <std::size_t N>
    [[nodiscard]] dare::FaceMatrixStencil<GridType, SC, N> operator*(
        const dare::GridVector<GridType, SC, N>& field) const;

    /*!
     * @brief a convenience function for accessing the Apply  function
     * @tparam ...Args parameter types going in
     * @param ...args the parameter pack
     * @return a FaceMatrixStencil
     */
    template <typename... Args>
    auto operator()(const Args&... args) const;

    /*!
     * @brief Applies the TVD scheme to a set of parameters and computes the Fluxes at the faces
     * @tparam ...Args parameter types to go in
     * @param ...args parameters to evaluation
     * @return A FaceMatrixStencils representing the fluxes at the faces
     * \note This function already applies the velocity at the faces, do not add those to the arguments!
     * 
     * @code{.cpp}
     *     // assuming that epsilon, rho and cp are parameters of type
     *     // GridVector, double, None or dare::Vector<N, SC>
     *     // we can compute the fluxes of epsilon * rho * cp * T * u by:
     *     J = tvd.Apply(epsilon, rho, cp, T->GetField()->GetDataVector(1))
     * @endcode
     */
    template <typename... Args>
    auto Apply(const Args&... args) const;

    /*!
     * @brief returns the velocity values
     */
    [[nodiscard]] const dare::FaceValueStencil<GridType, SC, 1>& GetVelocityStencil() const;

    /*!
     * @brief returns local index
     */
    [[nodiscard]] const Index& GetIndex() const;

private:
    /*!
     * @brief computes an extended stencil (close and far) by multiplying the relevant values
     * @tparam ...Args parameter pack types
     * @tparam N number of components
     * @param for_num_components just provided for giving the number of components, not used!
     * @param ...values the values which will be multiplied at the faces
     * @return an extended stencil (std::pair, where first -> close FaceValueStencil and second -> far FaceValueStencil)
     */
    template <std::size_t N, typename... Args>
    dare::ExtendedStencil<dare::CenterValueStencil<GridType, SC, N>>
    ComputeExtendedValueStencil(bool ignore_last,
                                const dare::GridVector<GridType, SC, N>& for_num_components,
                                const Args&... values) const;

    Index ind;                                         //!< triplet of indices
    dare::FaceValueStencil<GridType, SC, 1> velocity;  //!< stencil with velocity
    dare::Vector<NUM_FACES, bool> upwind;              //!< identifier for upwind at each face
    const GridRepresentation* grep;                    //!< grid information
};

namespace detail {

/*!
 * @brief end point for parameter unpacking of the extended stencils
 * @tparam SC scalar type
 * @tparam ...Args parameter types
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * @param grep representation of the target grid
 * @param ind triplet of indices of the center value
 * @param s stencil to multiply with
 * @param ...args value for multiplying with the stencil
 * 
 * This is the endpoint for unpacking of the parameters. It also serves
 * as an elegant way for checking, if we actually can use a provided value
 * for determining the extended stencil. If no appropriate overload is provided,
 * the compiler will end up here and throw an error in the static_assert!
 */
template <std::size_t Dim, typename SC, std::size_t N, typename... Args>
void free_tvd_cartesian_get_extended_stencil(
    const typename dare::Cartesian<Dim>::Representation& grep,
    const typename dare::Cartesian<Dim>::Index& ind,
    bool ignore_last,
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>* s,
    const Args&... args) {
    static_assert(sizeof...(args) == 0, "Cannot interpret the provided arguments!");
}

/*!
 * @brief overload for a multiplying the stencil with the appropriate GridVector values
 * @tparam SC scalar type
 * @tparam ...Args parameter types
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * @param grep representation of the target grid
 * @param ind triplet of indices of the center value
 * @param s stencil to multiply with
 * @param data reference to the grid vector
 * @param ...args remaining values for multiplying with the stencil
 */
template <std::size_t Dim, typename SC, std::size_t N, typename... Args>
void free_tvd_cartesian_get_extended_stencil(
    const typename dare::Cartesian<Dim>::Representation& grep,
    const typename dare::Cartesian<Dim>::Index& ind,
    bool ignore_last,
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>* s,
    const typename dare::GridVector<dare::Cartesian<Dim>, SC, N>& data,
    const Args&... args) {
    if constexpr (sizeof...(args) == 0) {
        if (ignore_last)
            return;
    }
    free_tvd_cartesian_get_extended_stencil(grep, ind, ignore_last, s, args...);
    s->first *= dare::InterpolateToCenterStencil(grep, ind, data, 0);
    s->second *= dare::InterpolateToCenterStencil(grep, ind, data, 1);
}

/*!
 * @brief overload for a multiplying the stencil with the appropriate GridVector values
 * @tparam SC scalar type
 * @tparam ...Args parameter types
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * @param grep representation of the target grid
 * @param ind triplet of indices of the center value
 * @param s stencil to multiply with
 * @param data pointer to the grid vector
 * @param ...args remaining values for multiplying with the stencil
 */
template <std::size_t Dim, typename SC, std::size_t N, typename... Args>
void free_tvd_cartesian_get_extended_stencil(
    const typename dare::Cartesian<Dim>::Representation& grep,
    const typename dare::Cartesian<Dim>::Index& ind,
    bool ignore_last,
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>* s,
    const typename dare::GridVector<dare::Cartesian<Dim>, SC, N>* data,
    const Args&... args) {
    free_tvd_cartesian_get_extended_stencil(grep, ind, ignore_last, s, *data, args...);
}

/*!
 * @brief overload for a multiplying the stencil with the appropriate Field values
 * @tparam SC scalar type
 * @tparam ...Args parameter types
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * @param grep representation of the target grid
 * @param ind triplet of indices of the center value
 * @param s stencil to multiply with
 * @param data reference to the field
 * @param ...args remaining values for multiplying with the stencil
 * 
 * \note The newest timestep of the Field will be used to compute the stencil!
 */
template <std::size_t Dim, typename SC, std::size_t N, typename... Args>
void free_tvd_cartesian_get_extended_stencil(
    const typename dare::Cartesian<Dim>::Representation& grep,
    const typename dare::Cartesian<Dim>::Index& ind,
    bool ignore_last,
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>* s,
    const typename dare::Field<dare::Cartesian<Dim>, SC, N>& data,
    const Args&... args) {
    if constexpr (sizeof...(args) == 0) {
        if (ignore_last)
            return;
    }
    free_tvd_cartesian_get_extended_stencil(grep, ind, ignore_last, s, data.GetDataVector(), args...);
}

/*!
 * @brief overload for a multiplying the stencil with the appropriate Field values
 * @tparam SC scalar type
 * @tparam ...Args parameter types
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * @param grep representation of the target grid
 * @param ind triplet of indices of the center value
 * @param s stencil to multiply with
 * @param data pointer to the field
 * @param ...args remaining values for multiplying with the stencil
 *
 * \note The newest timestep of the Field will be used to compute the stencil!
 */
template <std::size_t Dim, typename SC, std::size_t N, typename... Args>
void free_tvd_cartesian_get_extended_stencil(
    const typename dare::Cartesian<Dim>::Representation& grep,
    const typename dare::Cartesian<Dim>::Index& ind,
    bool ignore_last,
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>* s,
    const typename dare::Field<dare::Cartesian<Dim>, SC, N>* data,
    const Args&... args) {
    free_tvd_cartesian_get_extended_stencil(grep, ind, ignore_last, s, *data, args...);
}

/*!
 * @brief overload for a multiplying the stencil with a certain value
 * @tparam SC scalar type
 * @tparam ...Args parameter types
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * @param grep representation of the target grid
 * @param ind triplet of indices of the center value
 * @param s stencil to multiply with
 * @param value value to multiply the stencil with
 * @param ...args remaining values for multiplying with the stencil
 */
template <std::size_t Dim, typename SC, std::size_t N, typename... Args>
void free_tvd_cartesian_get_extended_stencil(
    const typename dare::Cartesian<Dim>::Representation& grep,
    const typename dare::Cartesian<Dim>::Index& ind,
    bool ignore_last,
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>* s,
    SC value,
    const Args&... args) {
    if constexpr (sizeof...(args) == 0) {
        if (ignore_last)
            return;
    }
    s->first *= value;
    s->second *= value;
    free_tvd_cartesian_get_extended_stencil(grep, ind, ignore_last, s, args...);
}

/*!
 * @brief overload for ignoring a None value
 * @tparam SC scalar type
 * @tparam ...Args parameter types
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * @param grep representation of the target grid
 * @param ind triplet of indices of the center value
 * @param s stencil to multiply with
 * @param value value to multiply the stencil with
 * @param ...args remaining values for multiplying with the stencil
 */
template <std::size_t Dim, typename SC, std::size_t N, typename... Args>
void free_tvd_cartesian_get_extended_stencil(
    const typename dare::Cartesian<Dim>::Representation& grep,
    const typename dare::Cartesian<Dim>::Index& ind,
    bool ignore_last,
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>* s,
    dare::None value,
    const Args&... args) {
    free_tvd_cartesian_get_extended_stencil(grep, ind, ignore_last, s, args...);
}

/*!
 * @brief a helper for extracting the number of component values
 * @tparam T
 * 
 * This is the default and will stop compilation, because no appropriate
 * specialization was found to get the value
 */
template<typename T>
struct free_tvd_cartesian_extract_num_components {
    /*!
     * @brief compile time function for getting the number of components (stops compilation)
     */
    static constexpr std::size_t GetValue() {
        static_assert(dare::always_false<T>, "Cannot extract the number of components from provided type");
        return 0;
    }
};

/*!
 * @brief a helper for extracting the number of component values
 * @tparam SC type of scalar value
 * @tparam Dim dimension fo the Cartesian grid
 * @tparam N number of components
 */
template <std::size_t Dim, typename SC, std::size_t N>
struct free_tvd_cartesian_extract_num_components<dare::GridVector<dare::Cartesian<Dim>, SC, N>> {
    /*!
     * @brief compile time function for getting the number of components
     */
    static constexpr std::size_t GetValue() {
        return N;
    }
};

/*!
 * @brief a helper for extracting the number of component values
 * @tparam SC type of scalar value
 * @tparam Dim dimension fo the Cartesian grid
 * @tparam N number of components
 */
template <std::size_t Dim, typename SC, std::size_t N>
struct free_tvd_cartesian_extract_num_components<
    dare::ExtendedStencil<dare::CenterValueStencil<dare::Cartesian<Dim>, SC, N>>> {
    /*!
     * @brief compile time function for getting the number of components
     */
    static constexpr std::size_t GetValue() {
        return N;
    }
};

}  // namespace detail

}  // end namespace dare

#include "TVD_Cartesian.inl"

#endif  // GRID_CARTESIAN_TVD_CARTESIAN_H_
