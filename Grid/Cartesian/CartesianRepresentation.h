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

#ifndef GRID_CARTESIAN_CARTESIANREPRESENTATION_H_
#define GRID_CARTESIAN_CARTESIANREPRESENTATION_H_

#include <map>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>

#include "MPI/HaloBuffer.h"
#include "Utilities/InitializationTracker.h"
namespace dare::Grid {

// forward declaration for the grid
template <std::size_t Dim>
class Cartesian;

// forward declarations for CartesianNeighbor
enum class CartesianNeighbor : char;
char ToNum(CartesianNeighbor pos);
char ToFace(CartesianNeighbor face);
char ToNormal(CartesianNeighbor nb);

/*!
 * @brief Representation of Cartesian grid
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 * @tparam SC scalar type
 * @tparam Dim dimension of grid
 */
template <std::size_t Dim>
class CartesianRepresentation : public dare::utils::InitializationTracker {
public:
    using GridType = Cartesian<Dim>;
    using LocalOrdinalType = typename GridType::LocalOrdinalType;
    using GlobalOrdinalType = typename GridType::GlobalOrdinalType;
    using ScalarType = typename GridType::ScalarType;
    using LO = LocalOrdinalType;
    using GO = GlobalOrdinalType;
    using SC = ScalarType;
    using VecLO = typename GridType::VecLO;
    using VecGO = typename GridType::VecGO;
    using VecSC = typename GridType::VecSC;
    using Index = typename GridType::Index;
    using IndexGlobal = typename GridType::IndexGlobal;
    using Options = typename GridType::Options;

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
    CartesianRepresentation(const CartesianRepresentation<Dim>& other) = default;

    /*!
     * @brief default copy assignment operator
     * @param other instance to copy from
     */
    CartesianRepresentation<Dim>&
    operator=(const CartesianRepresentation<Dim>& other) = default;

    /*!
     * @brief provides spatial position of the specified cell
     * @param ind local index
     * @return vector of size Dim with coordinates
     */
    VecSC GetCoordinatesCenter(const Index& ind) const;

    /*!
     * @brief provides spatial position of specified face
     * @param ind local index
     * @param cnb identifier
     * @return
     */
    VecSC GetCoordinatesFace(const Index& ind, CartesianNeighbor cnb) const;

    /*!
     * @brief returns the distances between two cell centers
     * @return
     */
    const VecSC& GetDistances() const;

    /*!
     * @brief returns the face areas in each dimension
     */
    const VecSC& GetFaceArea() const;

    /*!
     * @brief returns volume of cell
     * @param local_ordinal local index of cell
     */
    SC GetCellVolume(LO local_ordinal = 0) const;

    /*!
     * @brief number of cells in local subgrid including ghost/halo cells
     */
    LO GetNumberLocalCells() const;

    /*!
     * @brief number of cells in local subgrid excluding ghost/halo cells
     */
    LO GetNumberLocalCellsInternal() const;

    /*!
     * @brief number of cells in global grid including ghost/halo cells
     */
    GO GetNumberGlobalCells() const;

    /*!
     * @brief number of cells in global grid excluding ghost/halo cells
     */
    GO GetNumberGlobalCellsInternal() const;

    /*!
     * @brief returns the struct with provided options
     */
    const Options& GetOptions() const;

    /*!
     * @brief adds ghost/halo cells to ordinal
     * @param n_internal local ordinal without ghost/halo cells
     */
    LO MapInternalToLocal(LO n_internal) const;

    /*!
     * @brief removes ghost/halo cells from local index
     * @param ind_local local index
     */
    Index MapLocalToInternal(Index ind_local) const;

    /*!
     * @brief adds ghost/halo cells to local index
     * @param ind_local local index
     */
    Index MapInternalToLocal(Index ind_internal) const;

    /*!
     * @brief removes ghost cells from global index
     * @param ind_global global index
     */
    IndexGlobal MapGlobalToInternal(IndexGlobal ind_global) const;

    /*!
     * @brief adds ghost cells to global index
     * @param ind_global global index
     */
    IndexGlobal MapInternalToGlobal(IndexGlobal ind_internal) const;

    LO MapIndexToOrdinalLocal(const Index& ind) const;
    LO MapIndexToOrdinalLocalInternal(const Index& ind) const;
    GO MapIndexToOrdinalGlobal(const IndexGlobal& ind) const;
    GO MapIndexToOrdinalGlobalInternal(const IndexGlobal& ind) const;

    LO MapGlobalToLocal(GO id_glob) const;
    Index MapGlobalToLocal(const IndexGlobal& ind_glob) const;
    LO MapGlobalToLocalInternal(GO id_glob) const;
    GO MapLocalToGlobalInternal(LO id_glob) const;

    /*!
     * @brief tests, if a global ordinals belongs to subgrid
     * @param id_glob global ordinal
     * \note assumes the grid including ghost cells
     */
    bool IsLocal(GO id_glob) const;

    /*!
     * @brief tests, if a global index belongs to subgrid
     * @param ind_glob global index
     * \note assumes the grid including ghost cells
     */
    bool IsLocal(const IndexGlobal& ind_glob) const;

    /*!
     * @brief tests, if the ordinal belongs to internal subgrid
     * @param id_glob global internal ordinal
     */
    bool IsLocalInternal(GO id_glob) const;

    /*!
     * @brief tests, if index belongs to internal subgrid
     * @param ind_glob global internal index
     */
    bool IsLocalInternal(IndexGlobal ind_glob) const;

    bool IsInternal(const Index& ind_loc) const;
    bool IsInternal(const IndexGlobal& ind_glob) const;
    bool IsInternal(LO n) const;

    Index MapOrdinalToIndexLocal(const LO n_loc) const;
    Index MapOrdinalToIndexLocalInternal(const LO n_loc) const;
    IndexGlobal MapOrdinalToIndexGlobal(GO n_loc) const;
    IndexGlobal MapOrdinalToIndexGlobalInternal(GO n_loc) const;

    IndexGlobal MapLocalToGlobal(const Index& ind) const;
    GO MapLocalToGlobal(LO n_loc) const;

    /*!
     * @brief returns local Index of cell containing the point
     * @param point point in space
     * Make sure that the point is on the grid!
     */
    Index GetCell(VecSC point) const;

    /*!
     * @brief provides resolution of local grid
     * @return local index with number of cells in each direction
     */
    const Index& GetLocalResolution() const;

    /*!
     * @brief provides resolution of local grid
     * @return local index with number of cells in each direction
     */
    const Index& GetLocalResolutionInternal() const;

    /*!
     * @brief provides resolution of global grid
     * @return global index with number of cells in each direction
     */
    const IndexGlobal& GetGlobalResolution() const;

    /*!
     * @brief provides resolution of global grid
     * @return global index with number of cells in each direction
     */
    const IndexGlobal& GetGlobalResolutionInternal() const;

    /*!
     * @brief prints distribution of domains without ghost cells to file as CSV
     * @param fname filename (should be CSV)
     */
    void PrintDistribution(std::string fname) const;

    /*!
     * @brief Provides Halobuffer for exchange of data across processes
     */
    mpi::HaloBuffer<SC>& GetHaloBuffer();

    /*!
     * @brief provides name of the grid including options
     */
    std::string GetName() const;

    /*!
     * @brief returns dimensional offset cell to origin
     */
    const VecSC& GetOffsetSize() const;

    /*!
     * @brief returns offset to origin in cells
     * @return 
     */
    const VecGO& GetOffsetCells() const;

    /*!
     * @brief returns the number of ghost cells
     */
    [[nodiscard]] LO GetNumberGhostCells() const;

    /*!
     * @brief provides identifier for the boundaries of the domain
     */
    [[nodiscard]] uint8_t GetBoundaryID() const;

    /*!
     * @brief identifier of periodicity
     */
    [[nodiscard]] const VecLO& GetPeriodicity() const;

private:
    /*!
     * @brief Tests if grid is initialized
     * @param function name of function to show in error message
     */
    void TestIfInitialized(std::string function) const;

    std::string name;                  //!< Identifier of the grid instance
    const Cartesian<Dim>* grid;        //!< pointer to grid
    Options options;                   //!< grid specific options
    VecLO resolution_local;            //!< resolution of the local grid
    VecLO resolution_local_internal;   //!< resolution of the local grid without halo/ghost cells
    VecGO resolution_global;           //!< resolution of the global grid
    VecGO resolution_global_internal;  //!< resolution of the global grid without halo/ghost cells
    VecGO offset_cells;                //!< offset of the cells from local to global
    VecSC offset_size;                 //!< offset in size from local to global

    VecGO hierarchic_sum_glob;           //!< precomputed values to account for ordering
    VecGO hierarchic_sum_glob_internal;  //!< same as above, just for the internal cells
    VecLO hierarchic_sum_loc;            //!< precomputed values to account for ordering
    VecLO hierarchic_sum_loc_internal;   //!< same as above, just for the internal cells

    mpi::HaloBuffer<SC> halo_buffer;  //!< maps with Halo Buffers
};

}  // namespace dare::Grid

#include "CartesianRepresentation.inl"

#endif  // GRID_CARTESIAN_CARTESIANREPRESENTATION_H_
