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

#ifndef GRID_CARTESIAN_CARTESIANMESHUTILS_H_
#define GRID_CARTESIAN_CARTESIANMESHUTILS_H_

#include <iostream>
#include <list>
#include <string>
#include <array>

#include "Utilities/CompileTimeFunctions.h"

namespace dare::Grid {

namespace details::Cartesian {
/*!
 * @brief A small object to avoid double definition of grids of same names
 */
class AllocationManager {
    static std::list<std::string> reg;  //!< register for grid names
public:
    /*!
     * @brief registering a name
     * @param gname name of the grid to register
     * @return true, if not registered before, false if another grid exists
     */
    static bool RegisterGrid(const std::string& gname);

    /*!
     * @brief removing grid from register
     * @param gname name of the grid to deregister
     * @return true, if successful, false if gridname unknown
     */
    static bool DeregisterGrid(const std::string& gname);
};

}  // end namespace details::Cartesian

/*!
 * @brief identifiers for neighbor IDs of a Cartesian cell
 */
enum class CartesianNeighbor : char {
    CENTER = 0,
    WEST = 1,
    EAST = 2,
    SOUTH = 3,
    NORTH = 4,
    BOTTOM = 5,
    TOP = 6,
    FOURD_LOW = 7,
    FOURD_UP = 8,
};

/*!
 * @brief converts enum to char
 * @param pos Neighbor enum
 */
[[nodiscard]] inline char ToNum(CartesianNeighbor pos) {
    return static_cast<char>(pos);
}

/*!
 * @brief converts enum to char excluding center
 * @param face enum of face
 * This is mainly meant for array access of faces!
 */
[[nodiscard]] inline char ToFace(CartesianNeighbor face) {
#ifndef DARE_NDEBUG
    if (face == CartesianNeighbor::CENTER) {
        std::cerr << "In " << __func__ << ": Center is not a face!\n";
    }
#endif
    return ToNum(face) - 1;
}

/*!
 * @brief Converts enum to normal
 * @param nb neighbor id
 */
[[nodiscard]] inline char ToNormal(CartesianNeighbor nb) {
    char n = (ToNum(nb) % 2 == 0) - (ToNum(nb) % 2 > 0) - (nb == CartesianNeighbor::CENTER);
    return n;
}

template <std::size_t N>
using CartesianRangeType = std::array<CartesianNeighbor, N>;

/*!
 * @brief provides a range for looping through positions
 * @tparam Dim dimension of the grid
 * @return array with cartesian neighbors
 */
template <std::size_t Dim>
constexpr CartesianRangeType<Dim * 2 + 1> GetCartesianPositionRange() {
    if constexpr (Dim == 0) {
        return {CartesianNeighbor::CENTER};
    } else {
        return dare::utils::convert_tuple_to_array(
            std::tuple_cat(GetCartesianPositionRange<Dim - 1>(),
                           std::make_tuple(static_cast<CartesianNeighbor>(Dim * 2 - 1),
                                           static_cast<CartesianNeighbor>(Dim * 2))));
    }
}

/*!
 * @brief provides a range for looping through faces
 * @tparam Dim dimension of the grid
 * @return array with cartesian neighbors
 */
template <std::size_t Dim>
constexpr CartesianRangeType<Dim * 2> GetCartesianFaceRange() {
    static_assert(Dim > 0, "Cannot provide any faces for a 0D Cartesian grid");
    if constexpr (Dim == 1) {
        return {CartesianNeighbor::WEST, CartesianNeighbor::EAST};
    } else {
        return dare::utils::convert_tuple_to_array(
            std::tuple_cat(GetCartesianFaceRange<Dim - 1>(),
                           std::make_tuple(static_cast<CartesianNeighbor>(Dim * 2 - 1),
                                           static_cast<CartesianNeighbor>(Dim * 2))));
    }
}

/*!
 * \brief converts numerical id to CartesianNeighbor enum
 * @param id numerical identifier
 */
[[nodiscard]] inline CartesianNeighbor ToCartesianNeighbor(char id) {
    using CNB = CartesianNeighbor;
    switch (id) {
    case 0:
        return CNB::CENTER;
        break;
    case 1:
        return CNB::WEST;
        break;
    case 2:
        return CNB::EAST;
        break;
    case 3:
        return CNB::SOUTH;
        break;
    case 4:
        return CNB::NORTH;
        break;
    case 5:
        return CNB::BOTTOM;
        break;
    case 6:
        return CNB::TOP;
        break;
    case 7:
        return CNB::FOURD_LOW;
        break;
    case 8:
        return CNB::FOURD_UP;
        break;
    }
    ERROR << "invalid ID provided (" << id << ")" << ERROR_CLOSE;
    return CNB::FOURD_UP;  // most unlikely to be ever used
}

/*!
 * @brief converts an integer to enum
 * @tparam ID signed 8 bit integer
 */
template <char ID>
[[nodiscard]] constexpr CartesianNeighbor ToCartesianNeighbor() {
    static_assert(ID >= 0 && ID < 9);
    using CNB = CartesianNeighbor;
    if constexpr (ID == 0)
        return CNB::CENTER;
    else if constexpr (1)
        return CNB::WEST;
    else if constexpr (2)
        return CNB::EAST;
    else if constexpr (3)
        return CNB::SOUTH;
    else if constexpr (4)
        return CNB::NORTH;
    else if constexpr (5)
        return CNB::BOTTOM;
    else if constexpr (6)
        return CNB::TOP;
    else if constexpr (7)
        return CNB::FOURD_LOW;
    else if constexpr (8)
        return CNB::FOURD_UP;
}

/*!
 * \brief maps a Cartesian face identifier to a dimension
 * @param face face to refer to
 * 
 * this maps the faces of a Cartesian grid to a dimension, e.g.
 * EAST/WEST -> 0
 * SOUTH/NORTH -> 1
 * BOTTOM/TOP -> 2
 */
[[nodiscard]] inline std::size_t MapCartesianFaceToDim(CartesianNeighbor face) {
    using size_t = std::size_t;
#ifndef DARE_NDEBUG
    if (face == CartesianNeighbor::CENTER) {
        ERROR << "The input has to be a Cartesian face and may not be CENTER!" << ERROR_CLOSE;
    }
#endif
    return (static_cast<size_t>(face) - static_cast <size_t>(1)) / static_cast<size_t>(2);
}

/*!
 * \brief maps a Cartesian face identifier to a dimension at compile time
 * @param face face to refer to
 *
 * this maps the faces of a Cartesian grid to a dimension, e.g.
 * EAST/WEST -> 0
 * SOUTH/NORTH -> 1
 * BOTTOM/TOP -> 2
 */
template<CartesianNeighbor CNB>
constexpr std::size_t MapCartesianFaceToDim() {
    using size_t = std::size_t;
    static_assert(CNB != CartesianNeighbor::CENTER, "The input has to be a Cartesian face and may not be CENTER!");
    return (static_cast<size_t>(CNB) - static_cast<size_t>(1)) / static_cast<size_t>(2);
}

}  // namespace dare::Grid

#endif  // GRID_CARTESIAN_CARTESIANMESHUTILS_H_
