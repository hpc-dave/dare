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

#ifndef IO_VTKOPTIONS_CARTESIAN_H_
#define IO_VTKOPTIONS_CARTESIAN_H_
#include <vtkInformation.h>
#include <vtkMPIController.h>
#include <vtkProgrammableFilter.h>
#include <vtkStructuredGrid.h>

#include <utility>
#include <string>
#include <memory>

#include "Grid/Cartesian.h"
#include "VTKOptions.h"

namespace dare::io {

namespace details::Cartesian {

/*!
 * @brief arguments for the programmable filter
 */
struct PFArgs {
    explicit PFArgs(const VTKExtent& lextent) : local_extent(lextent) {}
    vtkNew<vtkProgrammableFilter> pf;   //!< programmable filter
    VTKExtent local_extent;             //!< local extent
};

/*!
 * @brief function to set up the vtkFilter
 * @param arg pointer to the arguments
 */
void inline execute(void* arg) {
    PFArgs* args = reinterpret_cast<PFArgs*>(arg);
    int extent[6] = {0};
    for (std::size_t n{0}; n < 6; n++)
        extent[n] = args->local_extent[n];
    auto output_tmp = args->pf->GetOutput();
    auto input_tmp = args->pf->GetInput();
    vtkStructuredGrid* output = dynamic_cast<vtkStructuredGrid*>(output_tmp);
    vtkStructuredGrid* input = dynamic_cast<vtkStructuredGrid*>(input_tmp);
    output->ShallowCopy(input);
    output->SetExtent(extent);
}

/*!
 * @brief small object for mapping grid ordinals to vtk ordinals
 * @tparam Dim dimension of the grid
 * This is required, since the Cartesian grid orders x-y-z, whereas
 * vtk orders z-y-x
 */
template <std::size_t Dim>
class Mapper {
public:
    using CartesianGrid = Grid::Cartesian<Dim>;
    using Representation = typename CartesianGrid::Representation;
    using LO = typename CartesianGrid::LocalOrdinalType;
    using Index = typename Representation::Index;

    /*!
     * @brief constructor
     * @param representation representation of the Cartesian grid
     */
    explicit Mapper(const Representation* representation)
        : grep(representation), hsum(1, 1, 1) {
        Index res = grep->GetLocalResolution();
        dare::utils::Vector<3, vtkOrdinal> res_3d(1, 1, 1);
        for (std::size_t d{0}; d < Dim; d++) {
            res_3d[d] = res[d];
        }
        hsum[1] = res_3d[0];              // jmax
        hsum[2] = res_3d[0] * res_3d[1];  // imax * jmax
    }

    /*!
     * @brief maps a local ordinal to the vtk ordinal
     * @param local_ordinal local ordinal of the grid
     * @return mapped ordinal
     */
    [[nodiscard]] vtkOrdinal operator()(LO local_ordinal) const {
        Index ind = grep->MapOrdinalToIndexLocal(local_ordinal);
        dare::utils::Vector<3, vtkOrdinal> ind_3d;
        for (std::size_t d{0}; d < Dim; d++) {
            ind_3d[d] = ind[d];
        }
        vtkOrdinal mapped_ordinal = hsum[0] * ind_3d[0] + hsum[1] * ind_3d[1] + hsum[2] * ind_3d[2];
        return mapped_ordinal;
    }

private:
    const Representation* grep;               //!< representation of the grid
    dare::utils::Vector<3, vtkOrdinal> hsum;  //!< hierarchical sum for mapping
};

}  // end namespace details::Cartesian

/*!
 * @brief specialization for the Cartesian grid type
 * @tparam Dim dimension of the Cartesian grid
 */
template <std::size_t Dim>
struct VTKOptions<Grid::Cartesian<Dim>> {
    using CartesianGrid = Grid::Cartesian<Dim>;
    using Representation = typename CartesianGrid::Representation;
    using GridType = vtkStructuredGrid;
    using LO = typename CartesianGrid::LocalOrdinalType;
    using SC = typename CartesianGrid::ScalarType;
    using vtkOrdinal = vtkIdType;
    using Mapper = typename details::Cartesian::Mapper<Dim>;
    using Writer = typename VTKWriterMapper<GridType>::type;

    /*!
     * @brief provides the grid for IO
     * @param grep representation of the grid
     * @param vtkgrid pointer to the grid which requires the information
     * VTK always requires a 3d version of the grid (using 0 for missing dimensions) and
     * sorts strongly with i, then j and finally z. Additionally, VTK requires
     * the node points, which form the cells and not the cells itself, although the extent
     * is given in number of cells and not points. Somebody understand that....
     */
    static bool AllocateGrid(const Representation& grep, GridType* vtkgrid) {
        using Index = typename Representation::Index;
        Index res = grep.GetLocalResolution();
        dare::utils::Vector<3, vtkOrdinal> res_3d(1, 1, 1);
        for (std::size_t d{0}; d < Dim; d++) {
            res_3d[d] = res[d] + 1;  // +1 for the number of points, required by VTK
        }

        // Here we compute some relevant data, including the cell dimensions and offset
        // Note, that we need to move the offset by half a cell size in each dimension
        // since the offset refers to the cell center!
        auto dn = grep.GetDistances();
        auto offset = grep.GetOffsetSize();
        offset -= static_cast<SC>(grep.GetNumberGhostCells()) * dn;
        dare::utils::Vector<3, double> offset_3d(0., 0., 0.), dn_3d(0., 0., 0.);
        dare::utils::Vector<3, int> offset_cells(0, 0, 0);
        for (std::size_t d{0}; d < Dim; d++) {
            dn_3d[d] = dn[d];
            offset_3d[d] = offset[d];
            offset_cells[d] = grep.GetOffsetCells()[d];
        }
        VTKExtent extent = GetPExtentGlobal(grep);

        // for 1D grids, we allocate a pseudo 2D grid
        if constexpr (Dim == 1) {
            res_3d.j() = 2;
            dn_3d.y() = dn_3d.x();
        }

        vtkNew<vtkPoints> points;
        points->Allocate(res_3d.i() * res_3d.j() * res_3d.k());
        vtkgrid->SetExtent(extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);

        int rank{0};
        int nranks{0};
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nranks);

        int count = 0;
        for (int k = 0; k < res_3d.k(); k++) {
            for (int j = 0; j < res_3d.j(); j++) {
                for (int i = 0; i < res_3d.i(); i++) {
                    vtkOrdinal n = count;
                    // double tuple[3];
                    double x = dn_3d.x() * i + offset_3d.x();
                    double y = dn_3d.y() * j + offset_3d.y();
                    double z = dn_3d.z() * k + offset_3d.z();
                    points->InsertPoint(n, x, y, z);
                    count++;
                }
            }
        }

        vtkgrid->SetPoints(points);

        vtkOrdinal num_cells_glob = vtkgrid->GetNumberOfCells();
        bool success = num_cells_glob == grep.GetNumberGlobalCells();
        if (!success) {
            ERROR << "Number of cells are incompatible! VTK computed "
                  << num_cells_glob << " vs the number of global cells: "
                  << grep.GetNumberGlobalCells() << ERROR_CLOSE;
        }
        return success;
    }

    /*!
     * @brief adds specialized data to the parallel writer
     * @param grep grid representation
     * @param exec_man execution manager
     * @param data vtkGridtype
     * @param writer the parallel writer
     * @return data which needs to persist until the files are written
     * Here, the halo cells are explicitly removed from the displayed extents int the root file.
     * This does not mean, that they don't exist in the separate pieces!
     */
    static std::unique_ptr<details::Cartesian::PFArgs> AddDataToWriter(const Representation& grep,
                                dare::mpi::ExecutionManager* exec_man,
                                vtkStructuredGrid* data,
                                Writer* writer) {
        using G = CartesianGrid;
        auto args = std::make_unique<details::Cartesian::PFArgs>(GetPExtentLocal(grep));
        VTKExtent extent_plocal = GetPExtentLocal(grep);
        // for improved output we remove the halo cell layer from
        // the extent in the parallel file
        // Warning! This only applies to the halo cells, NOT the ghost cells!
        VTKExtent halo_correction(0, 0, 0, 0, 0, 0);
        LO num_ghost = grep.GetNumberGhostCells();
        uint8_t b_id = grep.GetBoundaryID();
        halo_correction[0] = !(b_id & G::BOUNDARIES_WEST) * num_ghost;
        halo_correction[1] = !(b_id & G::BOUNDARIES_EAST) * (-num_ghost);
        if constexpr (Dim > 1) {
            halo_correction[2] = !(b_id & G::BOUNDARIES_SOUTH) * num_ghost;
            halo_correction[3] = !(b_id & G::BOUNDARIES_NORTH) * (-num_ghost);
        }
        if constexpr(Dim > 2) {
            halo_correction[4] = !(b_id & G::BOUNDARIES_BOTTOM) * num_ghost;
            halo_correction[5] = !(b_id & G::BOUNDARIES_TOP) * (-num_ghost);
        }
        extent_plocal += halo_correction;

        args->pf->SetExecuteMethod(details::Cartesian::execute, args.get());
        args->pf->SetInputData(data);
        writer->SetInputConnection(args->pf->GetOutputPort());
        writer->SetNumberOfPieces(exec_man->GetNumberProcesses());
        writer->SetStartPiece(exec_man->GetRank());
        writer->SetEndPiece(exec_man->GetRank());
        writer->SetGhostLevel(grep.GetNumberGhostCells());
        writer->SetPPieceExtent(extent_plocal, exec_man);
        return args;
    }

    /*!
     * @brief computes the global extent of the grid
     * @param grep representation of the grid
     */
    static VTKExtent GetPExtentGlobal(const Representation& grep) {
        using Index = typename Representation::IndexGlobal;
        Index res = grep.GetGlobalResolution();
        dare::utils::Vector<3, int> res_3d(0, 0, 0);
        // VTK always expects 3D data, so we need accomodate that
        for (std::size_t d{0}; d < Dim; d++) {
            res_3d[d] = static_cast<int>(res[d]);
        }
        VTKExtent extent(0, 0, 0, 0, 0, 0);

        // in x-direction
        extent[1] += res_3d.i();

        // in y-direction
        if constexpr (Dim > 1)
            extent[3] += res_3d.j();
        else
            extent[3] = 1;  // pseudo 2D grid for 1D

        // in z-direction
        if constexpr (Dim > 2)
            // extent[1] += res_3d.k();
            extent[5] += res_3d.k();

        return extent;
    }

    /*!
     * @brief returns the local extent of the grid including halo cells
     * @param grep grid representation
     */
    static VTKExtent GetPExtentLocal(const Representation& grep) {
        using Index = typename Representation::Index;

        Index res = grep.GetLocalResolution();
        dare::utils::Vector<3, int> res_3d(0, 0, 0);
        dare::utils::Vector<3, int> offset_cells(0, 0, 0);
        // VTK always expects 3D data, so we need accomodate that
        for (std::size_t d{0}; d < Dim; d++) {
            res_3d[d] = static_cast<int>(res[d]);
            offset_cells[d] = grep.GetOffsetCells()[d];
        }
        VTKExtent extent;
        extent[0] = extent[1] = offset_cells[0];
        extent[2] = extent[3] = offset_cells[1];
        extent[4] = extent[5] = offset_cells[2];

        // in x-direction
        extent[1] += res_3d.i();

        // in y-direction
        if constexpr (Dim > 1)
            extent[3] += res_3d.j();
        else
            extent[3] = 1;

        // in z-direction
        if constexpr (Dim > 2)
            extent[5] += res_3d.k();

        return extent;
    }
};

}  // end namespace dare::io

#endif  // IO_VTKOPTIONS_CARTESIAN_H_
