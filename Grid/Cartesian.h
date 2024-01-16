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

#ifndef GRID_CARTESIAN_H_
#define GRID_CARTESIAN_H_

#include <unordered_map>

#include "MPI/ExecutionManager.h"
#include "Utilities/Vector.h"
#include "Utilities/InitializationTracker.h"
#include "Grid/DefaultTypes.h"
#include "Grid/CartesianDistribution.h"
#include "Grid/CartesianRepresentation.h"

namespace dare::Grid {

enum class CartesianNeighbor : char {
    CENTER = 0,
    WEST = 1,
    EAST = 2,
    SOUTH = 3,
    NORTH = 4,
    BOTTOM = 5,
    TOP = 6,
    FOURD_LOW = 7,
    FOURD_UP = 8
};

[[nodiscard]] inline char ToNum(CartesianNeighbor pos) {
    return static_cast<char>(pos);
}
[[nodiscard]] inline char ToFace(CartesianNeighbor face) {
#ifndef DARE_NDEBUG
    if (face == CartesianNeighbor::CENTER) {
        std::cerr << "In " << __func__ << ": Center is not a face!\n";
    }
#endif
    return ToNum(face) - 1;
}

[[nodiscard]] inline char ToNormal(CartesianNeighbor nb) {
    char n = (ToNum(nb) % 2 == 0) - (ToNum(nb) % 2 > 0) - (nb == CartesianNeighbor::CENTER);
    return n;
}

template <std::size_t Dim>
class Cartesian : public dare::utils::InitializationTracker {
public:
    enum {
        BOUNDARIES_NONE = 0b00000000,
        BOUNDARIES_WEST = 0b00000001,
        BOUNDARIES_EAST = 0b00000010,
        BOUNDARIES_SOUTH = 0b00000100,
        BOUNDARIES_NORTH = 0b00001000,
        BOUNDARIES_BOTTOM = 0b00010000,
        BOUNDARIES_TOP = 0b00100000
    };

    enum {
        PERIODIC_NONE = 0b00000000,
        PERIODIC_X = 0b00000001,
        PERIODIC_Y = 0b00000010,
        PERIODIC_Z = 0b00000100
    };

    using LO = dare::defaults::LocalOrdinalType;
    using GO = dare::defaults::GlobalOrdinalType;
    using SC = dare::defaults::ScalarType;
    using GlobalOrdinalType = GO;
    using LocalOrdinalType = LO;
    static const std::size_t Dimension = Dim;
    static const std::size_t STENCIL_SIZE = Dim * 2 + 1;
    static const std::size_t NUM_FACES = Dim * 2;
    using ScalarType = SC;
    using VecGO = utils::Vector<Dim, GO>;
    using VecLO = utils::Vector<Dim, LO>;
    using VecSC = utils::Vector<Dim, SC>;
    using Representation = CartesianRepresentation<Dim>;
    using Options = VecLO;
    using Index = VecLO;
    using IndexGlobal = VecGO;
    using NeighborID = CartesianNeighbor;

    template < typename O>
    struct GetIndexType {
        using type = utils::Vector<Dim, O>;
    };

    /*!
     * \brief default constructor
     */
    Cartesian();

    template <typename Distributor>
    Cartesian(
        mpi::ExecutionManager* exec_man,
        const VecGO& resolution,
        const VecSC& size,
        const LO num_ghost,
        const VecLO& periodic,
        Distributor distributor);

    Cartesian(
        mpi::ExecutionManager* exec_man,
        const VecGO& resolution,
        const VecSC& size,
        const LO num_ghost,
        const VecLO& periodic);

    Cartesian(
        mpi::ExecutionManager* exec_man,
        const VecGO& resolution,
        const VecSC& size,
        const LO num_ghost);

    Cartesian(const Cartesian<Dim>& other) = delete;
    void operator=(const Cartesian<Dim>& other) = delete;

    /*!
     * \brief destructor
     */
    virtual ~Cartesian();

    /*!
     * \brief returns Representation object
     * \brief opt Options for the target grid
     */
    [[nodiscard]] Representation GetRepresentation(Options opt);

    /*!
     * \brief returns cell width
     */
    [[nodiscard]] const VecSC& GetCellWidth() const;

    /*!
     * \brief returns face area
     */
    [[nodiscard]] const VecSC& GetFaceArea() const;

    /*!
     * \brief returns cell volume
     */
    [[nodiscard]] SC GetCellVolume() const;

    /*!
     * \brief return number of ghost cells
     */
    [[nodiscard]] LO GetNumGhost() const;

    /*!
     * \brief return offset of subgrid in cells
     */
    [[nodiscard]] const VecGO& GetOffsetCells() const;

    /*!
     * \brief return offset of subgrid in distance
     */
    [[nodiscard]] const VecSC& GetOffsetSize() const;

    /*!
     * @brief returns the local resolution
     */
    [[nodiscard]] const VecLO& GetLocalResolution() const;
    /*!
     * @brief returns global resolution
     */
    [[nodiscard]] const VecGO& GetGlobalResolution() const;

    /*!
     * @brief returns local size
     */
    [[nodiscard]] const VecSC& GetLocalSize() const;
    /*!
     * @brief returns global size
     */
    [[nodiscard]] const VecSC& GetGlobalSize() const;

    /*!
     * @brief provides boundary identifier
     */
    [[nodiscard]] uint8_t GetBoundaryID() const;

    /*!
     * @brief returns periodicity
     */
    [[nodiscard]] const VecLO& GetPeriodicity() const;
    /*!
     * @brief quick inquiriy, if the grid is periodic
     */
    [[nodiscard]] bool IsPeriodic() const;

    /*!
     * @brief returns pointer to execution manager
     */
    [[nodiscard]] mpi::ExecutionManager* GetExecutionManager() const;

private:
    VecLO resolution_local;     //!< resolution of local internal grid
    VecGO resolution_global;    //!< resolution of global internal grid
    VecGO offset_cells;         //!< offset of local grid to origin

    VecSC size_local;           //!< dimensions of the local grid
    VecSC size_global;          //!< dimensions of the global grid
    VecSC offset_size;          //!< offset to origin

    VecSC cell_width;           //!< width of a cell in each dimension
    VecSC face_area;            //!< face area of a cell in each dimension
    SC cell_volume;             //!< cell volume

    LO num_ghost;               //!< number of ghost cells

    uint8_t id_boundaries;      //!< identifier of boundaries
    VecLO periodicity;          //!< identifier of periodicity

    std::unordered_map<Options, Representation> map_representations;  //!< representation storage

    mpi::ExecutionManager* exec_man;    //!< pointer to execution manager
};

}  // namespace dare::Grid

#include "Cartesian.inl"

#endif  // GRID_CARTESIAN_H_
