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

namespace dare::Grid {

template <std::size_t Dim>
CartesianRepresentation<Dim>::CartesianRepresentation()
    : grid(nullptr) {}

template <std::size_t Dim>
CartesianRepresentation<Dim>::CartesianRepresentation(const GridType* grid,
                                                      typename GridType::Options opt)
            : grid(grid),
              options(opt) {
    /*
     * Here we set all the determine all the required information for this instance
     * of the grid, including information containing staggered grids. Following steps are
     * conducted:
     * - get information from general grid
     * - determine additional cells in the case of the staggered grid
     * - add ghost/halo cells
     * - precompute helpers for more efficient conversion of cell ids to indices
     */

    // for now, the boundary id cannot handle more than 8 bits, which seems fine for any intended application
    // either way, the user should be warned before doing something stupid
    static_assert(Dim <= 4,
        "Currently, grids cannot have more than 4 dimensions. Not sure what you're trying to do here, weirdo?!");
    char staggered_IDs[] = {'X', 'Y', 'Z'};
    name = grid->GetName();
    for (std::size_t d{0}; d < Dim; d++) {
        if (opt[d] != 0) {
            name += "_";
            name += staggered_IDs[d];
        }
    }

    // Get information from main grid
    resolution_local = grid->GetLocalResolution();
    resolution_global = grid->GetGlobalResolution();
    offset_cells = grid->GetOffsetCells();
    offset_size = grid->GetOffsetSize();

    // Determine additional cell for staggered grid at the upper ends of the grid
    // Note, that this is only required for non-periodical grids
    for (std::size_t dim{0}; dim < Dim; dim++) {
        if (options[dim] == 0 || grid->GetPeriodicity()[dim] != 0)
            continue;

        if (grid->GetBoundaryID() & (1 << (dim * 2 + 1))) {
            resolution_local[dim] += 1;
        }
        resolution_global[dim] += 1;
        offset_size[dim] -= grid->GetCellWidth()[dim] * 0.5;
    }

    // Adding Ghost/Halo cells to the grids
    resolution_local_internal = resolution_local;
    resolution_global_internal = resolution_global;

    for (std::size_t dim{0}; dim < Dim; dim++) {
        resolution_local[dim] += 2 * grid->GetNumGhost();
        resolution_global[dim] += 2 * grid->GetNumGhost();
    }

    // Precomputing helpers for conversion of cell ID to indices and vice-versa
    for (std::size_t dim{0}; dim < (Dim - 1); dim++) {
        hierarchic_sum_loc_internal[dim] = resolution_local_internal[dim + 1];
        for (std::size_t n{dim + 2}; n < Dim; n++) {
            hierarchic_sum_loc_internal[dim] *= resolution_local_internal[n];
        }
        hierarchic_sum_glob_internal[dim] = resolution_global_internal[dim + 1];
        for (std::size_t n{dim + 2}; n < Dim; n++) {
            hierarchic_sum_glob_internal[dim] *= resolution_global_internal[n];
        }
    }
    hierarchic_sum_loc_internal[Dim - 1] = 1;
    hierarchic_sum_glob_internal[Dim - 1] = 1;

    for (std::size_t dim{0}; dim < (Dim - 1); dim++) {
        hierarchic_sum_loc[dim] = resolution_local[dim + 1];
        for (std::size_t n{dim + 2}; n < Dim; n++) {
            hierarchic_sum_loc[dim] *= resolution_local[n];
        }
        hierarchic_sum_glob[dim] = resolution_global[dim + 1];
        for (std::size_t n{dim + 2}; n < Dim; n++) {
            hierarchic_sum_glob[dim] *= resolution_global[n];
        }
    }
    hierarchic_sum_loc[Dim - 1] = 1;
    hierarchic_sum_glob[Dim - 1] = 1;

    this->dare::utils::InitializationTracker::Initialize();

    // Prepare halo-buffers
    std::vector<GO> required_halo_IDs;

    // calculating approximate number of halo cells (ignoring double count of edges)
    LO num_cells{0};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        // each dimension has two faces to account for
        LO num_face{1};
        for (std::size_t n{0}; n < Dim; n++) {
            if (n != dim)
                num_face *= resolution_local[n];
        }
        num_cells += 2 * grid->GetNumGhost() * num_face;
    }
    required_halo_IDs.reserve(num_cells);

    std::unordered_map<GO, GO> map_periodic;
    for (LO n{0}; n < GetNumberLocalCells(); n++) {
        Index ind = MapOrdinalToIndexLocal(n);
        IndexGlobal ind_glob = MapLocalToGlobal(ind);
        GO id_glob = MapIndexToOrdinalGlobal(ind_glob);
        if (!IsInternal(ind) && IsInternal(ind_glob)) {
                // The cell is locally not internal but globally it is an internal cell
                // -> classic halo cell
                required_halo_IDs.push_back(id_glob);
        }
        if (grid->IsPeriodic()) {
            // In the periodic case, we need to check if the cells is potentially
            // one that in absolute counting is outside of the global grid, but
            // with the periodic mapping actually reflects an internal cell.
            // Additionally, we have to check that the cells are also located at
            // a location, that can be mapped to an internal cell, to do so
            // we count the number of dimensions where the cell exceeds internal region:
            // if more than 1 condition is met, the cell is on an edge or face
            // or whatever that would be in 4D
            // Those do not have any counterpart which it could be mapped to and
            // need to be ignored
            //
            // Especially in the periodic case, we need to take into account, that the
            // periodic cells at the opposing end (potentially on another rank) need to be mapped
            // to internal cells e.g.
            //      rank 0                                      rank 1
            // ---------------------               ----------------------------------
            //| 0 | 1 || 2 | 3 |                       | n-4 | n-3 || n-2 | n-1 |
            // Here, cells 2 & 3 have to be mapped with n-2 & n-1, as well as
            // n-4 & n-3 to 0 & 1. However, due to the conditional above, the cells
            // 2, 3, n-4 and n-3 won't be mapped in if those are located on different ranks.
            // So we need to add those additionally to the map

            std::size_t count{0};
            Index periodicity = grid->GetPeriodicity();
            Index offset_periodic = periodicity * grid->GetNumGhost();
            for (std::size_t dim{0}; dim < Dim; dim++) {
                count += ind_glob[dim] < grid->GetNumGhost();
                count += ind_glob[dim] >= (resolution_global[dim] - grid->GetNumGhost());
            }

            if (count > 1)
                continue;

            bool low_period[Dim], high_period[Dim];
            bool is_internal = IsInternal(ind);
            bool is_not_receiving_cell{is_internal};
            for (std::size_t dim{0}; dim < Dim; dim++) {
                low_period[dim] = false;
                high_period[dim] = false;
                if (periodicity[dim] != 0) {
                    // Determine, if the periodic boundary is at lower or upper boundary
                    low_period[dim] = ind_glob[dim] < (grid->GetNumGhost() + offset_periodic[dim]);
                    high_period[dim] = ind_glob[dim] >= (resolution_global[dim]
                                                         - grid->GetNumGhost()
                                                         - offset_periodic[dim]);
                    if (low_period[dim] || high_period[dim]) {
                        IndexGlobal ind_period = ind_glob;
                        if (low_period[dim]) {
                            ind_period[dim] += resolution_global_internal[dim];
                        } else {
                            ind_period[dim] -= resolution_global_internal[dim];
                        }
                        GO id_glob_period = MapIndexToOrdinalGlobal(ind_period);
                        map_periodic[id_glob_period] = id_glob;
                    }
                }
                // is_not_receiving_cell &= !low_period[dim] || !high_period[dim];
                // is_not_receiving_cell &= !high_period[dim];
            }
            if (!is_not_receiving_cell && !IsInternal(ind_glob))
                required_halo_IDs.push_back(id_glob);
        }
    }

    auto is_internal_global = [&](GO id) {
        bool is_local_internal{IsLocal(id)};
        if (is_local_internal) {
            IndexGlobal ind_glob = MapOrdinalToIndexGlobal(id);
            Index ind_loc = MapGlobalToLocal(ind_glob);
            is_local_internal &= IsInternal(ind_loc);
        }
        return is_local_internal;
    };

    auto is_internal_local = [&](LO id) {
        auto ind = MapOrdinalToIndexLocal(id);
        return IsInternal(ind);
    };

    auto map_global_to_local = [&](GO id) {
        return MapGlobalToLocal(id);
    };
    halo_buffer.Initialize(grid->GetExecutionManager(), required_halo_IDs, map_periodic,
                           is_internal_global, is_internal_local, map_global_to_local);
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::VecSC
CartesianRepresentation<Dim>::GetCoordinatesCenter(const Index& ind) const {
    TestIfInitialized(__func__);
    VecSC pos;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        pos[dim] = offset_size[dim] + (ind[dim] - grid->GetNumGhost() + 0.5) * grid->GetCellWidth()[dim];
    }
    return pos;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::VecSC
CartesianRepresentation<Dim>::GetCoordinatesFace(const Index& ind, CartesianNeighbor cnb) const {
    TestIfInitialized(__func__);
    VecSC pos = GetCoordinatesCenter(ind);
    std::size_t dim = ToFace(cnb) / 2;
    pos[dim] += ToNormal(cnb) * 0.5 * this->GetDistances()[dim];
    return pos;
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::VecSC&
CartesianRepresentation<Dim>::GetDistances() const {
    return grid->GetCellWidth();
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::VecSC&
CartesianRepresentation<Dim>::GetFaceArea() const {
    return grid->GetFaceArea();
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::SC
CartesianRepresentation<Dim>::GetCellVolume(LO local_ordinal) const {
    return grid->GetCellVolume();
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::GetNumberLocalCellsInternal() const {
    TestIfInitialized(__func__);
    LO num_cells{1};
    for (LO dim : resolution_local)
        num_cells *= (dim - 2 * grid->GetNumGhost());
    return num_cells;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::GetNumberLocalCells() const {
    TestIfInitialized(__func__);
    LO num_cells{1};
    for (LO dim : resolution_local)
        num_cells *= dim;
    return num_cells;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::GO
CartesianRepresentation<Dim>::GetNumberGlobalCells() const {
    TestIfInitialized(__func__);
    GO num_cells{1};
    for (GO dim : resolution_global)
        num_cells *= dim;
    return num_cells;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::GO
CartesianRepresentation<Dim>::GetNumberGlobalCellsInternal() const {
    TestIfInitialized(__func__);
    GO num_cells{1};
    for (GO dim : resolution_global)
        num_cells *= (dim - 2 * grid->GetNumGhost());
    return num_cells;
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::Options&
CartesianRepresentation<Dim>::GetOptions() const {
    return options;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::MapInternalToLocal(LO n_internal) const {
    TestIfInitialized(__func__);
    if constexpr (Dim == 1) {
        return n_internal + grid->GetNumGhost();
    } else {
        Index ind;

        for (std::size_t dim{0}; dim < Dim; dim++) {
            ind[dim] = n_internal / hierarchic_sum_loc_internal[dim];
            n_internal -= ind[dim] * hierarchic_sum_loc_internal[dim];
        }

        for (LO& e : ind)
            e += grid->GetNumGhost();

        return MapIndexToOrdinalLocal(ind);
    }
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::Index
CartesianRepresentation<Dim>::MapLocalToInternal(Index ind_local) const {
    for (auto& e : ind_local) {
        e -= grid->GetNumGhost();
    }
    return ind_local;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::Index
CartesianRepresentation<Dim>::MapInternalToLocal(Index ind_local) const {
    for (auto& e : ind_local) {
        e += grid->GetNumGhost();
    }
    return ind_local;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::IndexGlobal
CartesianRepresentation<Dim>::MapGlobalToInternal(IndexGlobal ind_global) const {
    for (auto& e : ind_global) {
        e -= grid->GetNumGhost();
    }
    return ind_global;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::IndexGlobal
CartesianRepresentation<Dim>::MapInternalToGlobal(IndexGlobal ind_global) const {
    for (auto& e : ind_global) {
        e += grid->GetNumGhost();
    }
    return ind_global;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::MapIndexToOrdinalLocal(const Index& ind) const {
    TestIfInitialized(__func__);
    if constexpr (Dim == 1) {
        return ind[0];
    } else {
        LO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_loc[dim] * ind[dim];

        return index;
    }
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::MapIndexToOrdinalLocalInternal(const Index& ind) const {
    TestIfInitialized(__func__);
    if constexpr (Dim == 1) {
        return ind[0];
    } else {
        LO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_loc_internal[dim] * ind[dim];

        return index;
    }
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::GO
CartesianRepresentation<Dim>::MapIndexToOrdinalGlobal(const IndexGlobal& ind) const {
    TestIfInitialized(__func__);
    if constexpr (Dim == 1) {
        return ind[0];
    } else {
        LO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_glob[dim] * ind[dim];

        return index;
    }
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::GO
CartesianRepresentation<Dim>::MapIndexToOrdinalGlobalInternal(const IndexGlobal& ind, bool map_periodic) const {
    TestIfInitialized(__func__);
    if constexpr (Dim == 1) {
        GO o = ind[0];
        if (map_periodic && this->grid->IsPeriodic()) {
            GO o_max = this->GetGlobalResolutionInternal()[0];
            if (o < 0) {
                o = o_max + o;
            } else if (o >= o_max) {
                o = o - o_max;
            }
        }
        return o;
    } else {
        VecGO ind_manip{ind};
        if (map_periodic && this->grid->IsPeriodic()) {
            VecGO res{this->GetGlobalResolutionInternal()};
            for (std::size_t d{0}; d < Dim; d++) {
                GO o_max = res[d];
                if (ind[d] < 0) {
                    ind_manip[d] = o_max + ind[d];
                } else if (ind[d] >= o_max) {
                    ind_manip[d] = ind[d] - o_max;
                }
            }
        }
        GO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_glob_internal[dim] * ind_manip[dim];

        return index;
    }
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::IndexGlobal
CartesianRepresentation<Dim>::MapOrdinalToIndexGlobal(GO n_glob) const {
    TestIfInitialized(__func__);
    IndexGlobal ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_glob / hierarchic_sum_glob[dim];
        n_glob -= ind[dim] * hierarchic_sum_glob[dim];
    }
    return ind;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::IndexGlobal
CartesianRepresentation<Dim>::MapOrdinalToIndexGlobalInternal(GO n_glob) const {
    TestIfInitialized(__func__);
    IndexGlobal ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_glob / hierarchic_sum_glob_internal[dim];
        n_glob -= ind[dim] * hierarchic_sum_glob_internal[dim];
    }
    return ind;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::Index
CartesianRepresentation<Dim>::MapOrdinalToIndexLocal(LO n_loc) const {
    TestIfInitialized(__func__);
    Index ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_loc / hierarchic_sum_loc[dim];
        n_loc -= ind[dim] * hierarchic_sum_loc[dim];
    }
    return ind;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::Index
CartesianRepresentation<Dim>::MapOrdinalToIndexLocalInternal(LO n_loc) const {
    TestIfInitialized(__func__);
    Index ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_loc / hierarchic_sum_loc_internal[dim];
        n_loc -= ind[dim] * hierarchic_sum_loc_internal[dim];
    }
    return ind;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::MapGlobalToLocal(GO id_glob) const {
    TestIfInitialized(__func__);
    IndexGlobal ind_glob = MapOrdinalToIndexGlobal(id_glob);
    Index ind_loc = MapGlobalToLocal(ind_glob);
    return MapIndexToOrdinalLocal(ind_loc);
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::Index
CartesianRepresentation<Dim>::MapGlobalToLocal(const IndexGlobal& ind_glob) const {
    TestIfInitialized(__func__);
    Index ind_loc;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind_loc[dim] = ind_glob[dim] - offset_cells[dim];
    }
    return ind_loc;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::MapGlobalToLocalInternal(GO id_glob) const {
    IndexGlobal ind_glob = MapOrdinalToIndexGlobalInternal(id_glob);
    Index ind_loc = MapGlobalToLocal(ind_glob);
    return MapIndexToOrdinalLocalInternal(ind_loc);
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::GO
CartesianRepresentation<Dim>::MapLocalToGlobalInternal(LO id_loc) const {
    Index ind_loc = MapOrdinalToIndexLocalInternal(id_loc);
    IndexGlobal ind_glob = MapLocalToGlobal(ind_loc);
    return MapIndexToOrdinalGlobalInternal(ind_glob);
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::Index
CartesianRepresentation<Dim>::GetCell(typename CartesianRepresentation<Dim>::VecSC point) const {
    Index ind;

    // computation of the cell indices as floating point values
    // the offset size is referring to the origin of the INTERNAL grid
    // and therefore we need to correct for the ghost cell layer
    point -= offset_size;
    point += this->GetDistances() * static_cast<SC>(grid->GetNumGhost());
    point /= this->GetDistances();

    for (std::size_t d{0}; d < Dim; d++) {
        ind[d] = static_cast<LO>(point[d]);
    }

    return ind;
}

template <std::size_t Dim>
bool CartesianRepresentation<Dim>::IsLocal(GO id_glob) const {
    IndexGlobal ind_glob = MapOrdinalToIndexGlobal(id_glob);
    return IsLocal(ind_glob);
}

template <std::size_t Dim>
bool CartesianRepresentation<Dim>::IsLocalInternal(GO id_glob) const {
    IndexGlobal ind_glob = MapOrdinalToIndexGlobalInternal(id_glob);
    return IsLocalInternal(ind_glob);
}

template <std::size_t Dim>
bool CartesianRepresentation<Dim>::IsLocal(const IndexGlobal& ind_glob) const {
    TestIfInitialized(__func__);
    bool is_local{true};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        is_local &= ind_glob[dim] >= offset_cells[dim];
        is_local &= ind_glob[dim] < offset_cells[dim] + resolution_local[dim];
    }
    return is_local;
}

template <std::size_t Dim>
bool CartesianRepresentation<Dim>::IsLocalInternal(IndexGlobal ind_glob) const {
    TestIfInitialized(__func__);
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind_glob[dim] += 2 * grid->GetNumGhost();
    }
    return IsLocal(ind_glob);
}

template <std::size_t Dim>
bool CartesianRepresentation<Dim>::IsInternal(const Index& ind_loc) const {
    TestIfInitialized(__func__);
    bool is_internal{true};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        is_internal &= ind_loc[dim] >= (grid->GetNumGhost());
        is_internal &= ind_loc[dim] < (resolution_local[dim] - grid->GetNumGhost());
    }
    return is_internal;
}

template <std::size_t Dim>
bool CartesianRepresentation<Dim>::IsInternal(LO id) const {
    Index ind = MapOrdinalToIndexLocal(id);
    return IsInternal(ind);
}

template <std::size_t Dim>
bool CartesianRepresentation<Dim>::IsInternal(const IndexGlobal& ind_glob) const {
    TestIfInitialized(__func__);
    bool is_internal{true};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        is_internal &= ind_glob[dim] >= (grid->GetNumGhost());
        is_internal &= ind_glob[dim] < (resolution_global[dim] - grid->GetNumGhost());
    }
    return is_internal;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::IndexGlobal
CartesianRepresentation<Dim>::MapLocalToGlobal(const Index& ind) const {
    TestIfInitialized(__func__);
    IndexGlobal ind_g;
    for (std::size_t dim{0}; dim < Dim; dim++)
        ind_g[dim] = ind[dim] + offset_cells[dim];
    return ind_g;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::GO
CartesianRepresentation<Dim>::MapLocalToGlobal(LO n_loc) const {
    TestIfInitialized(__func__);
    Index ind_l = MapOrdinalToIndexLocal(n_loc);
    IndexGlobal ind_g = MapLocalToGlobal(ind_l);
    GO n_glob = MapIndexToOrdinalGlobal(ind_g);
    return n_glob;
}
template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::Index&
CartesianRepresentation<Dim>::GetLocalResolution() const {
    return resolution_local;
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::Index&
CartesianRepresentation<Dim>::GetLocalResolutionInternal() const {
    return resolution_local_internal;
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::IndexGlobal&
CartesianRepresentation<Dim>::GetGlobalResolution() const {
    return resolution_global;
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::IndexGlobal&
CartesianRepresentation<Dim>::GetGlobalResolutionInternal() const {
    return resolution_global_internal;
}

template <std::size_t Dim>
void CartesianRepresentation<Dim>::PrintDistribution(std::string fname) const {
    for (int n{0}; n < grid->GetExecutionManager()->GetNumberProcesses(); n++) {
        if (n == grid->GetExecutionManager()->GetRank()) {
            std::ofstream ofs;
            if (n == 0) {
                ofs.open(fname);
                ofs << "Proc,";
                for (std::size_t dim{0}; dim < Dim; dim++)
                    ofs << "Delta_" << dim << ',';
                for (std::size_t dim{0}; dim < Dim - 1; dim++)
                    ofs << "offset_" << dim << ',';
                ofs << "offset_" << Dim - 1 << '\n';
            } else {
                ofs.open(fname, std::ios_base::app);
            }
            if (!ofs)
                std::cout << "Could not open file " << fname << std::endl;
            ofs << n << ',';
            for (std::size_t dim{0}; dim < Dim; dim++)
                ofs << resolution_local[dim] - 2 * grid->GetNumGhost() << ',';
            for (std::size_t dim{0}; dim < Dim - 1; dim++)
                ofs << offset_cells[dim] << ',';
            ofs << offset_cells[Dim - 1] << '\n';
        }
        grid->GetExecutionManager()->Barrier();
    }
}

template <std::size_t Dim>
void CartesianRepresentation<Dim>::TestIfInitialized(std::string func) const {
#ifndef DARE_NDEBUG
    if (!dare::utils::InitializationTracker::IsInitialized()) {
        grid->GetExecutionManager()->Terminate(func, "Cannot continue without initialization");
    }
#endif
}

template <std::size_t Dim>
mpi::HaloBuffer<typename CartesianRepresentation<Dim>::SC>&
CartesianRepresentation<Dim>::GetHaloBuffer() {
    return halo_buffer;
}

template <std::size_t Dim>
std::string CartesianRepresentation<Dim>::GetName() const {
    return name;
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::VecSC&
CartesianRepresentation<Dim>::GetOffsetSize() const {
    return offset_size;
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::VecGO&
CartesianRepresentation<Dim>::GetOffsetCells() const {
    return offset_cells;
}

template <std::size_t Dim>
typename CartesianRepresentation<Dim>::LO
CartesianRepresentation<Dim>::GetNumberGhostCells() const {
    return grid->GetNumGhost();
}

template <std::size_t Dim>
uint8_t CartesianRepresentation<Dim>::GetBoundaryID() const {
    return grid->GetBoundaryID();
}

template <std::size_t Dim>
const typename CartesianRepresentation<Dim>::VecLO&
CartesianRepresentation<Dim>::GetPeriodicity() const {
    return grid->GetPeriodicity();
}

}  // namespace dare::Grid
