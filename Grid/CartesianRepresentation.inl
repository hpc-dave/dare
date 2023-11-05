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

namespace dare::Grid {

template <std::size_t Dim, class LO, class GO, class SC>
CartesianRepresentation<Dim, LO, GO, SC>::CartesianRepresentation(const GridType* grid,
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
        VecLO ind = MapCellToIndexLocal(n);
        if (!IsInternal(ind)) {
            VecGO ind_glob = MapLocalToGlobal(ind);
            GO id_glob = MapIndexToCellGlobal(ind_glob);
            if (IsInternal(ind_glob)) {
                required_halo_IDs.push_back(id_glob);
            } else if (grid->IsPeriodic()) {
                // counts the number of met conditions
                // if more than 1 condition is met, the cell is on an edge or face
                // or whatever that would be in 4D
                // Those do not have any counterpart which it could be mapped to and
                // need to be ignored
                std::size_t count{0};
                VecLO periodicity = grid->GetPeriodicity();
                for (std::size_t dim{0}; dim < Dim; dim++) {
                        count += ind_glob[dim] < grid->GetNumGhost();
                        count += ind_glob[dim] >= (resolution_global[dim] - grid->GetNumGhost());
                }
                if (count != 1)
                    continue;
                for (std::size_t dim{0}; dim < Dim; dim++) {
                    if (periodicity[dim] != 0) {
                        bool low_period{ind_glob[dim] < grid->GetNumGhost()};
                        bool high_period{ind_glob[dim] >= (resolution_global[dim] - grid->GetNumGhost())};
                        if (low_period || high_period) {
                            VecGO ind_period = ind_glob;
                            if (low_period) {
                                ind_period[dim] += resolution_global_internal[dim];
                            } else {
                                ind_period[dim] -= resolution_global_internal[dim];
                            }
                            GO id_glob_period = MapIndexToCellGlobal(ind_period);
                            required_halo_IDs.push_back(id_glob_period);
                            map_periodic[id_glob_period] = id_glob;
                            break;
                        }
                    }
                }
            }
        }
    }

    auto is_local = [&](GO id) {
        return IsLocal(id);
    };
    auto map_global_to_local = [&](GO id) {
        return MapGlobalToLocal(id);
    };
    halo_buffer.Initialize(grid->GetExecutionManager(), required_halo_IDs, map_periodic,
                           is_local, map_global_to_local);
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecSC
CartesianRepresentation<Dim, LO, GO, SC>::GetPositionCenter(const VecLO& ind) const {
    VecSC pos;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        pos[dim] = offset_size[dim] + (ind[dim] - grid->GetNumGhost() + 0.5) * grid->GetCellWidth()[dim];
    }
    return pos;
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecSC
CartesianRepresentation<Dim, LO, GO, SC>::GetPositionFace(const VecLO& ind, std::size_t dir) const {
    VecSC pos;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        pos[dim] = offset_size[dim] + (ind[dim] - grid->GetNumGhost()) * grid->GetCellWidth()[dim];
    }
    pos[dir] += grid->GetCellWidth()[dir] * ((ind[dir] + 1) % 2);
    return pos;
}

template <std::size_t Dim, class LO, class GO, class SC>
LO CartesianRepresentation<Dim, LO, GO, SC>::GetNumberLocalCellsInternal() const {
    LO num_cells{1};
    for (LO dim : resolution_local)
        num_cells *= (dim - 2 * grid->GetNumGhost());
    return num_cells;
}

template <std::size_t Dim, class LO, class GO, class SC>
LO CartesianRepresentation<Dim, LO, GO, SC>::GetNumberLocalCells() const {
    LO num_cells{1};
    for (LO dim : resolution_local)
        num_cells *= dim;
    return num_cells;
}

template <std::size_t Dim, class LO, class GO, class SC>
LO CartesianRepresentation<Dim, LO, GO, SC>::MapInternalToLocal(LO n_internal) const {
    if constexpr (Dim == 1) {
        return n_internal + grid->GetNumGhost();
    } else {
        VecLO ind;

        for (std::size_t dim{0}; dim < Dim; dim++) {
            ind[dim] = n_internal / hierarchic_sum_loc_internal[dim];
            n_internal -= ind[dim] * hierarchic_sum_loc_internal[dim];
        }

        for (LO& e : ind)
            e += grid->GetNumGhost();

        return MapIndexToCellLocal(ind);
    }
}

template <std::size_t Dim, class LO, class GO, class SC>
LO CartesianRepresentation<Dim, LO, GO, SC>::MapIndexToCellLocal(const VecLO& ind) const {
    if constexpr (Dim == 1) {
        return ind[0];
    } else {
        LO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_loc[dim] * ind[dim];

        return index;
    }
}

template <std::size_t Dim, class LO, class GO, class SC>
LO CartesianRepresentation<Dim, LO, GO, SC>::MapIndexToCellLocalInternal(const VecLO& ind) const {
    if constexpr (Dim == 1) {
        return ind[0];
    } else {
        LO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_loc_internal[dim] * ind[dim];

        return index;
    }
}

template <std::size_t Dim, class LO, class GO, class SC>
GO CartesianRepresentation<Dim, LO, GO, SC>::MapIndexToCellGlobal(const VecGO& ind) const {
    if constexpr (Dim == 1) {
        return ind[0];
    } else {
        LO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_glob[dim] * ind[dim];

        return index;
    }
}

template <std::size_t Dim, class LO, class GO, class SC>
GO CartesianRepresentation<Dim, LO, GO, SC>::MapIndexToCellGlobalInternal(const VecGO& ind) const {
    if constexpr (Dim == 1) {
        return ind[0];
    } else {
        LO index{0};
        for (std::size_t dim{0}; dim < Dim; dim++)
            index += hierarchic_sum_glob_internal[dim] * ind[dim];

        return index;
    }
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecGO
CartesianRepresentation<Dim, LO, GO, SC>::MapCellToIndexGlobal(GO n_glob) const {
    VecGO ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_glob / hierarchic_sum_glob[dim];
        n_glob -= ind[dim] * hierarchic_sum_glob[dim];
    }
    return ind;
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecGO
CartesianRepresentation<Dim, LO, GO, SC>::MapCellToIndexGlobalInternal(GO n_glob) const {
    VecGO ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_glob / hierarchic_sum_glob_internal[dim];
        n_glob -= ind[dim] * hierarchic_sum_glob_internal[dim];
    }
    return ind;
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecLO
CartesianRepresentation<Dim, LO, GO, SC>::MapCellToIndexLocal(LO n_loc) const {
    VecLO ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_loc / hierarchic_sum_loc[dim];
        n_loc -= ind[dim] * hierarchic_sum_loc[dim];
    }
    return ind;
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecLO
CartesianRepresentation<Dim, LO, GO, SC>::MapCellToIndexLocalInternal(LO n_loc) const {
    VecLO ind;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind[dim] = n_loc / hierarchic_sum_loc_internal[dim];
        n_loc -= ind[dim] * hierarchic_sum_loc_internal[dim];
    }
    return ind;
}

template <std::size_t Dim, class LO, class GO, class SC>
LO CartesianRepresentation<Dim, LO, GO, SC>::MapGlobalToLocal(GO id_glob) const {
    VecGO ind_glob = MapCellToIndexGlobal(id_glob);
    VecLO ind_loc = MapGlobalToLocal(ind_glob);
    return MapIndexToCellLocal(ind_loc);
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecLO
CartesianRepresentation<Dim, LO, GO, SC>::MapGlobalToLocal(const VecGO& ind_glob) const {
    VecLO ind_loc;
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind_loc[dim] = ind_glob[dim] - offset_cells[dim];
    }
    return ind_loc;
}

template <std::size_t Dim, class LO, class GO, class SC>
LO CartesianRepresentation<Dim, LO, GO, SC>::MapGlobalToLocalInternal(GO id_glob) const {
    VecGO ind_glob = MapCellToIndexGlobalInternal(id_glob);
    VecLO ind_loc = MapGlobalToLocal(ind_glob);
    return MapIndexToCellLocalInternal(ind_loc);
}

template <std::size_t Dim, class LO, class GO, class SC>
bool CartesianRepresentation<Dim, LO, GO, SC>::IsLocal(GO id_glob) const {
    VecGO ind_glob = MapCellToIndexGlobal(id_glob);
    return IsLocal(ind_glob);
}

template <std::size_t Dim, class LO, class GO, class SC>
bool CartesianRepresentation<Dim, LO, GO, SC>::IsLocalInternal(GO id_glob) const {
    VecGO ind_glob = MapCellToIndexGlobalInternal(id_glob);
    return IsLocalInternal(ind_glob);
}

template <std::size_t Dim, class LO, class GO, class SC>
bool CartesianRepresentation<Dim, LO, GO, SC>::IsLocal(const VecGO& ind_glob) const {
    bool is_local{true};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        is_local &= ind_glob[dim] >= offset_cells[dim];
        is_local &= ind_glob[dim] < offset_cells[dim] + resolution_local[dim];
    }
    return is_local;
}

template <std::size_t Dim, class LO, class GO, class SC>
bool CartesianRepresentation<Dim, LO, GO, SC>::IsLocalInternal(VecGO ind_glob) const {
    for (std::size_t dim{0}; dim < Dim; dim++) {
        ind_glob[dim] += 2 * grid->GetNumGhost();
    }
    return IsLocal(ind_glob);
}

template <std::size_t Dim, class LO, class GO, class SC>
bool CartesianRepresentation<Dim, LO, GO, SC>::IsInternal(const VecLO& ind_loc) const {
    bool is_internal{true};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        is_internal &= ind_loc[dim] >= (grid->GetNumGhost());
        is_internal &= ind_loc[dim] < (resolution_local[dim] - grid->GetNumGhost());
    }
    return is_internal;
}

template <std::size_t Dim, class LO, class GO, class SC>
bool CartesianRepresentation<Dim, LO, GO, SC>::IsInternal(LO id) const {
    VecLO ind = MapCellToIndexLocal(id);
    return IsInternal(ind);
}

template <std::size_t Dim, class LO, class GO, class SC>
bool CartesianRepresentation<Dim, LO, GO, SC>::IsInternal(const VecGO& ind_glob) const {
    bool is_internal{true};
    for (std::size_t dim{0}; dim < Dim; dim++) {
        is_internal &= ind_glob[dim] >= (grid->GetNumGhost());
        is_internal &= ind_glob[dim] < (resolution_global[dim] - grid->GetNumGhost());
    }
    return is_internal;
}

template <std::size_t Dim, class LO, class GO, class SC>
typename CartesianRepresentation<Dim, LO, GO, SC>::VecGO
CartesianRepresentation<Dim, LO, GO, SC>::MapLocalToGlobal(const VecLO& ind) const {
    VecGO ind_g;
    for (std::size_t dim{0}; dim < Dim; dim++)
        ind_g[dim] = ind[dim] + offset_cells[dim];
    return ind_g;
}

template <std::size_t Dim, class LO, class GO, class SC>
const typename CartesianRepresentation<Dim, LO, GO, SC>::VecLO&
CartesianRepresentation<Dim, LO, GO, SC>::GetLocalResolution() const {
    return resolution_local;
}

template <std::size_t Dim, class LO, class GO, class SC>
const typename CartesianRepresentation<Dim, LO, GO, SC>::VecGO&
CartesianRepresentation<Dim, LO, GO, SC>::GetGlobalResolution() const {
    return resolution_global;
}

template <std::size_t Dim, class LO, class GO, class SC>
void CartesianRepresentation<Dim, LO, GO, SC>::PrintDistribution(std::string fname) const {
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

}  // namespace dare::Grid
