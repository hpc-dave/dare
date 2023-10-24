/*
 * MIT License
 *
 * Copyright (c) 2022 David Rieder

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
#include <iostream>
#include <fstream>
#include <string>
namespace dare::Grid {

template <std::size_t Dim, class LO, class GO, class SC>
class Cartesian;

template <std::size_t Dim, class LO, class GO, class SC>
class CartesianRepresentation {
public:
    using GridType = Cartesian<Dim, LO, GO, SC>;
    using VecLO = typename GridType::VecLO;
    using VecGO = typename GridType::VecGO;
    using VecSC = typename GridType::VecSC;
    // make a dedicated object to deal with Halo information?

    explicit CartesianRepresentation(const GridType* grid, typename GridType::Options opt)
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
        }
        hierarchic_sum_loc_internal[Dim - 1] = 1;

        for (std::size_t dim{0}; dim < (Dim - 1); dim++) {
            hierarchic_sum_loc[dim] = resolution_local[dim + 1];
            for (std::size_t n{dim + 2}; n < Dim; n++) {
                hierarchic_sum_loc[dim] *= resolution_local[n];
            }
        }
        hierarchic_sum_loc[Dim - 1] = 1;
    }

    CartesianRepresentation(const CartesianRepresentation<Dim, LO, GO, SC>& other) = default;

    CartesianRepresentation<Dim, LO, GO, SC>&
    operator=(const CartesianRepresentation<Dim, LO, GO, SC>& other) = default;

    VecSC GetPositionCenter(const VecLO& ind) const {
        VecSC pos;
        for (std::size_t dim{0}; dim < Dim; dim++) {
            pos[dim] = offset_size[dim] + (ind[dim] - grid->GetNumGhost() + 0.5) * grid->GetCellWidth()[dim];
        }
        return pos;
    }

    VecSC GetPositionFace(const VecLO& ind, std::size_t dir) const {
        VecSC pos;
        for (std::size_t dim{0}; dim < Dim; dim++) {
            pos[dim] = offset_size[dim] + (ind[dim] - grid->GetNumGhost()) * grid->GetCellWidth()[dim];
        }
        pos[dir] += grid->GetCellWidth()[dir] * ((ind[dir] + 1) % 2);
        return pos;
    }

    LO GetNumberLocalCells() const {
        LO num_cells{1};
        for (LO dim : resolution_local)
            num_cells *= (dim - 2 * grid->GetNumGhost());
        return num_cells;
    }

    LO MapInternalToLocal(LO n_internal) {
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

    LO MapIndexToCellLocal(const VecLO& ind) const {
        if constexpr (Dim == 1) {
            return ind[0];
        } else {
            LO index{0};
            for (std::size_t dim{0}; dim < Dim; dim++)
                index += hierarchic_sum_loc[dim] * ind[dim];

            return index;
        }
    }

    VecLO MapCellToIndexLocal(const LO n_loc) const {
        VecLO ind;
        for (std::size_t dim{0}; dim < Dim; dim++) {
            ind[dim] = n_loc / hierarchic_sum_loc[dim];
            n_loc -= ind[dim] * hierarchic_sum_loc[dim];
        }
        return ind;
    }

    VecGO MapLocalToGlobal(const VecLO& ind) const {
        VecGO ind_g;
        for (std::size_t dim{0}; dim < Dim; dim++)
            ind_g[dim] = ind[dim] + offset_cells[dim];
        return ind_g;
    }

    const VecLO& GetLocalResolution() const {
        return resolution_local;
    }
    const VecGO& GetGlobalResolution() const {
        return resolution_global;
    }

    /*!
     * @brief prints distribution of domains without ghost cells to file as CSV
     * @param fname filename (should be CSV)
     */
    void PrintDistribution(std::string fname) const {
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
                    ofs << "offset_" << Dim-1 << '\n';
                } else {
                    ofs.open(fname, std::ios_base::app);
                }
                if (!ofs)
                    std::cout << "Could not open file " << fname << std::endl;
                ofs << n << ',';
                for (std::size_t dim{0}; dim < Dim; dim++)
                    ofs << resolution_local[dim] - 2* grid->GetNumGhost() << ',';
                for (std::size_t dim{0}; dim < Dim - 1; dim++)
                    ofs << offset_cells[dim] << ',';
                ofs << offset_cells[Dim - 1] << '\n';
            }
            grid->GetExecutionManager()->Barrier();
        }
    }

private:
    const Cartesian<Dim, LO, GO, SC>* grid;  //!< pointer to grid
    typename GridType::Options options;      //!< grid specific options
    VecLO resolution_local;                  //!< resolution of the local grid
    VecLO resolution_local_internal;         //!< resolution of the local grid without halo/ghost cells
    VecGO resolution_global;                 //!< resolution of the global grid
    VecGO resolution_global_internal;        //!< resolution of the global grid without halo/ghost cells
    VecGO offset_cells;                      //!< offset of the cells from local to global
    VecSC offset_size;                       //!< offset in size from local to global

    VecLO hierarchic_sum_loc;           //!< precomputed values to account for ordering
    VecLO hierarchic_sum_loc_internal;  //!< same as above, just for the internal cells
};

}  // namespace dare::Grid

#endif  // GRID_CARTESIANREPRESENTATION_H_
