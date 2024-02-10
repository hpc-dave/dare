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
Cartesian<Dim>::Cartesian() : exec_man(nullptr) {
    static_assert(std::is_signed_v<LO> && std::is_integral_v<LO>, "The local ordinal needs to be a signed integer!");
    static_assert(std::is_signed_v<GO> && std::is_integral_v<GO>, "The global ordinal needs to be a signed integer!");
}

template <std::size_t Dim>
template <typename Distributor>
Cartesian<Dim>::Cartesian(std::string _name,
                          mpi::ExecutionManager* _exec_man,
                          const typename Cartesian<Dim>::VecGO& res,
                          const typename Cartesian<Dim>::VecSC& size,
                          const LO _num_ghost,
                          const VecLO& periodic,
                          Distributor dist)
    : name(_name),
      resolution_global(res),
      size_global(size),
      cell_volume(1.),
      num_ghost(_num_ghost),
      id_boundaries(BOUNDARIES_NONE),
      periodicity(periodic),
      exec_man(_exec_man) {
    static_assert(std::is_signed_v<LO> && std::is_integral_v<LO>, "The local ordinal needs to be a signed integer!");
    static_assert(std::is_signed_v<GO> && std::is_integral_v<GO>, "The global ordinal needs to be a signed integer!");

    // test if grid with same name was already allocated and register this one
    if (!details::Cartesian::RegisterGrid(name)) {
        exec_man->Terminate(__func__, "double allocation");
    }

    // Test for some potentially problematic resolutions
    for (std::size_t dim{0}; dim < Dim; dim++) {
        if (res[dim] <= num_ghost) {
            std::ostringstream os;
            os << "The grid needs to have a higher resolution in each direction than the number of ghost/halo cells! "
               << "The resolution in Dimension " << dim << " was found to be " << res[dim]
               << " but should be > " << num_ghost << "!\n"
               << " Otherwise the halo buffers might become problematic\n";
            ERROR << os.str() << ERROR_CLOSE;
            exec_man->Terminate(__func__, "problematic resolution");
        }
    }
    /*
     * 1) compute constant variables like face area, cell volume and the like
     * 2) distribute the grid over all processors
     */
    for (std::size_t dim{0}; dim < Dim; dim++) {
        cell_width[dim] = size_global[dim] / resolution_global[dim];
    }

    for (const auto dn : cell_width)
        cell_volume *= dn;

    for (std::size_t dim{0}; dim < Dim; dim++)
        face_area[dim] = cell_volume / cell_width[dim];

    dist(exec_man, resolution_global, &resolution_local, &offset_cells);

    if (offset_cells[0] == 0)
        id_boundaries |= BOUNDARIES_WEST;
    if ((offset_cells[0] + resolution_local[0]) == resolution_global[0])
        id_boundaries |= BOUNDARIES_EAST;

    if constexpr (Dim > 1) {
        if (offset_cells[1] == 0)
            id_boundaries |= BOUNDARIES_SOUTH;
        if ((offset_cells[1] + resolution_local[1]) == resolution_global[1])
            id_boundaries |= BOUNDARIES_NORTH;
    }

    if constexpr (Dim > 2) {
        if (offset_cells[2] == 0)
            id_boundaries |= BOUNDARIES_BOTTOM;
        if ((offset_cells[2] + resolution_local[2]) == resolution_global[2])
            id_boundaries |= BOUNDARIES_TOP;
    }

    for (std::size_t dim{0}; dim < Dim; dim++) {
        offset_size[dim] = offset_cells[dim] * cell_width[dim];
        size_local[dim] = resolution_local[dim] * cell_width[dim];
    }

    this->dare::utils::InitializationTracker::Initialize();
}

template <std::size_t Dim>
template <typename Distributor>
Cartesian<Dim>::Cartesian(mpi::ExecutionManager* exec_man,
                          const typename Cartesian<Dim>::VecGO& res,
                          const typename Cartesian<Dim>::VecSC& size,
                          const LO num_ghost,
                          const VecLO& periodic,
                          Distributor dist)
                          : Cartesian("Cartesian", exec_man, res, size, num_ghost, periodic, dist) {
}

template <std::size_t Dim>
Cartesian<Dim>::Cartesian(mpi::ExecutionManager* _exec_man,
                          const typename Cartesian<Dim>::VecGO& res,
                          const typename Cartesian<Dim>::VecSC& size,
                          const LO _num_ghost,
                          const VecLO& periodic)
    : Cartesian(_exec_man,
                res,
                size,
                _num_ghost,
                periodic,
                [](mpi::ExecutionManager* a,
                   const typename Cartesian<Dim>::VecGO& b,
                   typename Cartesian<Dim>::VecLO* c,
                   typename Cartesian<Dim>::VecGO* d) {
                    dare::Grid::CartesianDistribution_MPI_Dims_create(a, b, c, d);
                }) {
}

template <std::size_t Dim>
Cartesian<Dim>::Cartesian(const std::string& _name,
                          mpi::ExecutionManager* _exec_man,
                          const typename Cartesian<Dim>::VecGO& res,
                          const typename Cartesian<Dim>::VecSC& size,
                          const LO _num_ghost,
                          const VecLO& periodic)
    : Cartesian(_name,
                _exec_man,
                res,
                size,
                _num_ghost,
                periodic,
                [](mpi::ExecutionManager* a,
                   const typename Cartesian<Dim>::VecGO& b,
                   typename Cartesian<Dim>::VecLO* c,
                   typename Cartesian<Dim>::VecGO* d) {
                    dare::Grid::CartesianDistribution_MPI_Dims_create(a, b, c, d);
                }) {
}

template <std::size_t Dim>
Cartesian<Dim>::Cartesian(mpi::ExecutionManager* _exec_man,
                                      const typename Cartesian<Dim>::VecGO& res,
                                      const typename Cartesian<Dim>::VecSC& size,
                                      const LO _num_ghost)
    : Cartesian(_exec_man,
                res,
                size,
                _num_ghost,
                VecLO()) {
}

template <std::size_t Dim>
Cartesian<Dim>::Cartesian(const std::string& _name,
                          mpi::ExecutionManager* _exec_man,
                          const typename Cartesian<Dim>::VecGO& res,
                          const typename Cartesian<Dim>::VecSC& size,
                          const LO _num_ghost)
    : Cartesian(_name,
                _exec_man,
                res,
                size,
                _num_ghost,
                VecLO()) {
}

template <std::size_t Dim>
Cartesian<Dim>::~Cartesian() {
    if (!details::Cartesian::DeregisterGrid(name)) {
        exec_man->Terminate(__func__, "register issue");
    }
}

template <std::size_t Dim>
typename Cartesian<Dim>::Representation Cartesian<Dim>::GetRepresentation(Options opt) {
    if (!this->dare::utils::InitializationTracker::IsInitialized()) {
        exec_man->Terminate(__func__, "Grid needs to be initialized before providing a Representation");
    }
    if (map_representations.find(opt) == map_representations.end()) {
        map_representations.insert({opt, Representation(this, opt)});
    }
    return map_representations[opt];
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecSC& Cartesian<Dim>::GetCellWidth() const {
    return cell_width;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecSC& Cartesian<Dim>::GetFaceArea() const {
    return face_area;
}

template <std::size_t Dim>
typename Cartesian<Dim>::SC
Cartesian<Dim>::GetCellVolume() const {
    return cell_volume;
}

template <std::size_t Dim>
typename Cartesian<Dim>::LO
Cartesian<Dim>::GetNumGhost() const {
    return num_ghost;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecGO& Cartesian<Dim>::GetOffsetCells() const {
    return offset_cells;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecSC& Cartesian<Dim>::GetOffsetSize() const {
    return offset_size;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecLO& Cartesian<Dim>::GetLocalResolution() const {
    return resolution_local;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecGO& Cartesian<Dim>::GetGlobalResolution() const {
    return resolution_global;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecSC& Cartesian<Dim>::GetLocalSize() const {
    return size_local;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecSC& Cartesian<Dim>::GetGlobalSize() const {
    return size_global;
}

template <std::size_t Dim>
uint8_t Cartesian<Dim>::GetBoundaryID() const {
    return id_boundaries;
}

template <std::size_t Dim>
const typename Cartesian<Dim>::VecLO& Cartesian<Dim>::GetPeriodicity() const {
    return periodicity;
}

template <std::size_t Dim>
bool Cartesian<Dim>::IsPeriodic() const {
    bool is_periodic{false};
    for (auto e : periodicity)
        is_periodic |= (e != 0);
    return is_periodic;
}

template <std::size_t Dim>
mpi::ExecutionManager* Cartesian<Dim>::GetExecutionManager() const {
    return exec_man;
}

template <std::size_t Dim>
std::string Cartesian<Dim>::GetName() const {
    return name;
}

}  // namespace dare::Grid
