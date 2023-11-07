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

#ifndef GRID_CARTESIANREPRESENTATION_H_
#define GRID_CARTESIANREPRESENTATION_H_

#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>

#include "../MPI/HaloBuffer.h"
#include "../Utilities/InitializationTracker.h"
namespace dare::Grid {

template <std::size_t Dim, class LO, class GO, class SC>
class Cartesian;

template <std::size_t Dim, class LO, class GO, class SC>
class CartesianRepresentation : public dare::utils::InitializationTracker {
public:
    using GridType = Cartesian<Dim, LO, GO, SC>;
    using VecLO = typename GridType::VecLO;
    using VecGO = typename GridType::VecGO;
    using VecSC = typename GridType::VecSC;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;

    /*!
     * @brief default construction without initialization
     */
    CartesianRepresentation();

    /*!
     * @brief constructor with initialization
     * @param grid pointer to the grid
     * @param opt options related to the grid
     */
    explicit CartesianRepresentation(const GridType* grid, typename GridType::Options opt);

    /*!
     * @brief default copy constructor
     * @param other instance to copy from
     */
    CartesianRepresentation(const CartesianRepresentation<Dim, LO, GO, SC>& other) = default;

    /*!
     * @brief default copy assignment operator
     * @param other instance to copy from
     */
    CartesianRepresentation<Dim, LO, GO, SC>&
    operator=(const CartesianRepresentation<Dim, LO, GO, SC>& other) = default;

    VecSC GetPositionCenter(const Index& ind) const;

    VecSC GetPositionFace(const Index& ind, std::size_t dir) const;

    LO GetNumberLocalCells() const;
    LO GetNumberLocalCellsInternal() const;

    LO MapInternalToLocal(LO n_internal) const;

    LO MapIndexToOrdinalLocal(const Index& ind) const;
    LO MapIndexToOrdinalLocalInternal(const Index& ind) const;
    GO MapIndexToOrdinalGlobal(const IndexGlobal& ind) const;
    GO MapIndexToOrdinalGlobalInternal(const IndexGlobal& ind) const;

    LO MapGlobalToLocal(GO id_glob) const;
    Index MapGlobalToLocal(const IndexGlobal& ind_glob) const;
    LO MapGlobalToLocalInternal(GO id_glob) const;

    bool IsLocal(GO id_glob) const;
    bool IsLocal(const IndexGlobal& ind_glob) const;
    bool IsLocalInternal(GO id_glob) const;
    bool IsLocalInternal(IndexGlobal ind_glob) const;

    bool IsInternal(const Index& ind_loc) const;
    bool IsInternal(const IndexGlobal& ind_glob) const;
    bool IsInternal(LO n) const;

    Index MapOrdinalToIndexLocal(const LO n_loc) const;
    Index MapOrdinalToIndexLocalInternal(const LO n_loc) const;
    IndexGlobal MapOrdinalToIndexGlobal(GO n_loc) const;
    IndexGlobal MapOrdinalToIndexGlobalInternal(GO n_loc) const;

    IndexGlobal MapLocalToGlobal(const Index& ind) const;

    const Index& GetLocalResolution() const;

    const IndexGlobal& GetGlobalResolution() const;

    /*!
     * @brief prints distribution of domains without ghost cells to file as CSV
     * @param fname filename (should be CSV)
     */
    void PrintDistribution(std::string fname) const;

private:
    void TestIfInitialized(std::string) const;

    const Cartesian<Dim, LO, GO, SC>* grid;  //!< pointer to grid
    typename GridType::Options options;      //!< grid specific options
    VecLO resolution_local;                  //!< resolution of the local grid
    VecLO resolution_local_internal;         //!< resolution of the local grid without halo/ghost cells
    VecGO resolution_global;                 //!< resolution of the global grid
    VecGO resolution_global_internal;        //!< resolution of the global grid without halo/ghost cells
    VecGO offset_cells;                      //!< offset of the cells from local to global
    VecSC offset_size;                       //!< offset in size from local to global

    VecGO hierarchic_sum_glob;                //!< precomputed values to account for ordering
    VecGO hierarchic_sum_glob_internal;       //!< same as above, just for the internal cells
    VecLO hierarchic_sum_loc;                 //!< precomputed values to account for ordering
    VecLO hierarchic_sum_loc_internal;        //!< same as above, just for the internal cells

    mpi::HaloBuffer<LO, GO, SC> halo_buffer;  //!< maps with Halo Buffers
};

}  // namespace dare::Grid

#include "CartesianRepresentation.inl"

#endif  // GRID_CARTESIANREPRESENTATION_H_
