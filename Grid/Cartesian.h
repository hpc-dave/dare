/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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

#include "../MPI/ExecutionManager.h"
#include "../Utilities/Vector.h"
#include "../Utilities/InitializationTracker.h"
#include "DefaultTypes.h"
#include "CartesianDistribution.h"
#include "CartesianRepresentation.h"

namespace dare::Grid {

template <std::size_t Dim,
          class LO = dare::Grid::details::LocalOrdinalType,
          class GO = dare::Grid::details::GlobalOrdinalType,
          class SC = double>
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

    using GlobalOrdinalType = GO;
    using LocalOrdinalType = LO;
    const std::size_t Dimension = Dim;
    using ScalarType = SC;
    using VecGO = utils::Vector<Dim, GO>;
    using VecLO = utils::Vector<Dim, LO>;
    using VecSC = utils::Vector<Dim, SC>;
    using Representation = CartesianRepresentation<Dim, LO, GO, SC>;
    using Options = VecLO;
    using Index = VecLO;
    using IndexGlobal = VecGO;

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

    Cartesian(const Cartesian<Dim, LO, GO, SC>& other) = delete;
    void operator=(const Cartesian<Dim, LO, GO, SC>& other) = delete;

    /*!
     * \brief destructor
     */
    virtual ~Cartesian();

    /*!
     * \brief returns Representation object
     * \brief opt Options for the target grid
     */
    Representation GetRepresentation(Options opt);

    /*!
     * \brief returns cell width
     */
    const VecSC& GetCellWidth() const;

    /*!
     * \brief returns face area
     */
    const VecSC& GetFaceArea() const;

    /*!
     * \brief returns cell volume
     */
    SC GetCellVolume() const;

    /*!
     * \brief return number of ghost cells
     */
    LO GetNumGhost() const;

    /*!
     * \brief return offset of subgrid in cells
     */
    const VecGO& GetOffsetCells() const;

    /*!
     * \brief return offset of subgrid in distance
     */
    const VecSC& GetOffsetSize() const;

    /*!
     * @brief returns the local resolution
     */
    const VecLO& GetLocalResolution() const;
    /*!
     * @brief returns global resolution
     */
    const VecGO& GetGlobalResolution() const;

    /*!
     * @brief returns local size
     */
    const VecSC& GetLocalSize() const;
    /*!
     * @brief returns global size
     */
    const VecSC& GetGlobalSize() const;

    /*!
     * @brief provides boundary identifier
     */
    uint8_t GetBoundaryID() const;
    /*!
     * @brief returns periodicity
     */
    const VecLO& GetPeriodicity() const;
    /*!
     * @brief quick inquiriy, if the grid is periodic
     */
    bool IsPeriodic() const;

    /*!
     * @brief returns pointer to execution manager
     */
    mpi::ExecutionManager* GetExecutionManager() const;

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
