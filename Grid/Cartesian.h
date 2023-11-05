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

#include "../MPI/ExecutionManager.h"
#include "../Utilities/Vector.h"
#include "CartesianDistribution.h"
#include "CartesianRepresentation.h"

namespace dare::Grid {

template <std::size_t Dim, class LO = int32_t, class GO = int64_t, class SC = double>
class Cartesian {
public:
    // enum class STAGGERED : uint8_t {
    //     SCALAR = 0,
    //     XSTAGGERED = 1,
    //     YSTAGGERED = 2,
    //     ZSTAGGERED = 3
    // };

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
    Representation GetRepresentation(Options opt) const;

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

    VecSC GetPosition(const VecLO& ind) const;

    const VecLO& GetLocalResolution() const;
    const VecGO& GetGlobalResolution() const;

    const VecSC& GetLocalSize() const;
    const VecSC& GetGlobalSize() const;

    uint8_t GetBoundaryID() const;
    const VecLO& GetPeriodicity() const;
    bool IsPeriodic() const;

    // mpi::ExecutionManager* GetExecutionManager();
    mpi::ExecutionManager* GetExecutionManager() const;

private:
    VecLO resolution_local;
    VecGO resolution_global;
    VecGO offset_cells;

    VecSC size_local;
    VecSC size_global;
    VecSC offset_size;

    VecSC cell_width;
    VecSC face_area;
    SC cell_volume;

    LO num_ghost;

    uint8_t id_boundaries;
    VecLO periodicity;

    mpi::ExecutionManager* exec_man;
};

}  // namespace dare::Grid

#include "Cartesian.inl"

#endif  // GRID_CARTESIAN_H_
