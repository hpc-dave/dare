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
#include <utility>
#include <vtkStructuredGrid.h>

#include "Grid/Cartesian.h"
#include "VTKOptions.h"

namespace dare::io {

template <std::size_t Dim>
struct VTKOptions<Grid::Cartesian<Dim>> {
    using CartesianGrid = Grid::Cartesian<Dim>;
    using Representation = typename CartesianGrid::Representation;
    using GridType = vtkStructuredGrid;
    using LO = typename CartesianGrid::LocalOrdinalType;
    using SC = typename CartesianGrid::ScalarType;
    using vtkOrdinal = vtkIdType;

    /*!
     * @brief just some default which serves as example
     * @param grep instance of the grid
     * @param local_ordinal input ordinal
     * @return output ordinal according to the requirements of VTK
     * In the VTK structured grid, the x-component advances fastest,
     * while the z-component advances slowest, inverse to how the
     * Cartesian grid deals with that
     */
    static LO Map(const Representation& grep, LO local_ordinal) {
        using Index = typename Representation::Index;
        Index res = grep.GetLocalResolution();
        Index ind = grep.MapOrdinalToIndexLocal(local_ordinal);
        dare::utils::Vector<3, LO> res_3d(1, 1, 1);
        dare::utils::Vector<3, LO> ind_3d(0, 0, 0);
        for (std::size_t d{0}; d < Dim; d++) {
            res_3d[d] = res[d];
            ind_3d[d] = ind[d];
        }

        // std::swap(res_3d[0], res_3d[2]);
        // std::swap(ind_3d[0], ind_3d[2]);

        dare::utils::Vector<3, LO> hsum(1, 1, 1);
        hsum[0] = res_3d[1] * res_3d[2];  // jmax * imax
        hsum[1] = res_3d[2];             // imax

        LO mapped_ordinal = hsum[0] * ind_3d[0] + hsum[1] * ind_3d[1] + hsum[2] * ind_3d[2];
        return mapped_ordinal;
        // using Index = typename Representation::IndexGlobal;
        // auto ind_local = grep.MapOrdinalToIndexLocal(local_ordinal);
        // Index ind = grep.MapLocalToGlobal(ind_local);
        // Index res = grep.GetGlobalResolution();
        // dare::utils::Vector<3, vtkOrdinal> res_3d(1, 1, 1);
        // dare::utils::Vector<3, vtkOrdinal> ind_3d(0, 0, 0);
        // for (std::size_t d{0}; d < Dim; d++) {
        //     res_3d[d] = static_cast<vtkOrdinal>(res[d]);
        //     ind_3d[d] = static_cast<vtkOrdinal>(ind[d]);
        // }

        // std::swap(res_3d[0], res_3d[2]);
        // std::swap(ind_3d[0], ind_3d[2]);

        // dare::utils::Vector<3, vtkOrdinal> hsum(1, 1, 1);
        // hsum[0] = res_3d[1] * res_3d[2];  // jmax * imax
        // hsum[1] = res_3d[2];              // imax

        // LO mapped_ordinal = hsum[0] * ind_3d[0] + hsum[1] * ind_3d[1] + hsum[2] * ind_3d[2];
        // return mapped_ordinal;
    }

    /*!
     * @brief provides the grid for IO
     * @param grep representation of the grid
     */
    static vtkNew<GridType> GetGrid(const Representation& grep) {
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
        // offset -= 0.5 * dn;
        // offset -= static_cast<SC>(grep.GetNumberGhostCells()) * dn;
        dare::utils::Vector<3, double> offset_3d(0., 0., 0.), dn_3d(0., 0., 0.);
        dare::utils::Vector<3, int> offset_cells(0, 0, 0);
        for (std::size_t d{0}; d < Dim; d++) {
            dn_3d[d] = dn[d];
            offset_3d[d] = offset[d];
            offset_cells[d] = grep.GetOffsetCells()[d];
        }
        dare::utils::Vector<6, int> extent;
        // extent[0] = extent[1] = offset_cells[2];
        // extent[2] = extent[3] = offset_cells[1];
        // extent[4] = extent[5] = offset_cells[0];
        // // in x-direction
        // extent[5] += res_3d.i() - 1;
        // // in y-direction
        // if constexpr (Dim > 1)
        //     extent[3] += res_3d.j() - 1;
        // // in z-direction
        // if constexpr (Dim > 2)
        //     extent[1] += res_3d.k() - 1;
        extent = GetPExtentGlobal(grep);

        vtkNew<vtkPoints> points;
        vtkNew<GridType> grid;
        // points->SetDataTypeToDouble();
        points->Allocate(res_3d.i() * res_3d.j() * res_3d.k());
        grid->SetExtent(extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);
        // grid->SetExtent(extent[0], extent[5], extent[2], extent[3], extent[4], extent[1]);

        int rank{0};
        int nranks{0};
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nranks);

        // if(rank == 0)
        //     std::cout << "size of vtkIdType: " << std::to_string(sizeof(vtkOrdinal) * 8) << " bit\n";
        int count = 0;
        for (int nr = 0; nr < nranks; nr++){
            if (nr == rank) {
                for (int k = 0; k < res_3d.k(); k++) {
                    for (int j = 0; j < res_3d.j(); j++) {
                        // if (rank == 1)
                        std::cout << "new layer" << std::endl;
                        // auto n_glob = grep.MapLocalToGlobal(count);
                        // count++;
                        for (int i = 0; i < res_3d.i(); i++) {
                            vtkOrdinal n = count;
                            // double tuple[3];
                            double x = dn_3d.x() * i + offset_3d.x();
                            double y = dn_3d.y() * j + offset_3d.y();
                            double z = dn_3d.z() * k + offset_3d.z();
                            // if (rank == 1)
                            std::cout << "index: " << n << " x: " << x << " y: " << y << " z: " << z << std::endl;
                            points->InsertPoint(n, x, y, z);
                            count++;
                        }
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        // int count = 0;
        // count = 0;
        // dare::utils::Vector<3, vtkOrdinal> hsum(1, 1, 1);
        // hsum[0] = res_3d[0] * res_3d[1];  // jmax * imax
        // hsum[1] = res_3d[0];              // imax

        // for (int i = 0; i < res_3d.i(); i++) {
        //     for (int j = 0; j < res_3d.j(); j++) {
        //         for (int k = 0; k < res_3d.k(); k++) {
        //             // if (rank == 1)
        //             //     std::cout << "new layer" << std::endl;
        //             // auto n_glob = grep.MapLocalToGlobal(count);
        //             // count++;
        //             // vtkOrdinal n = Map(grep, count);
        //             vtkOrdinal n = hsum[0] * k + hsum[1] * j + hsum[2] * i;
        //             // double tuple[3];
        //             double x = dn_3d.x() * i + offset_3d.x();
        //             double y = dn_3d.y() * j + offset_3d.y();
        //             double z = dn_3d.z() * k + offset_3d.z();
        //             // if (rank == 1)
        //             //     std::cout << "index: " << n << " x: " << x << " y: " << y << " z: " << z << std::endl;
        //             points->SetPoint(n, x, y, z);
        //             i_failing[n].SetValues(x,y,z);
        //             count++;
        //         }
        //     }
        // }
        // if(rank == 1)
        // for (std::size_t n = 0; n < i_working.size(); n++){
        //     std::cout << i_working[n] << " : " << i_failing[n] << std::endl;
        // }

        grid->SetPoints(points);
        std::ostringstream os;
        os << "process " << rank << ": local extent -> ";
        for (int i = 0; i < 6; i++)
            os << extent[i] << " ";
        os << '\n';
        // os << "process " << rank << ": global extent -> ";
        // for (int i = 0; i < 6; i++)
        //     os << global_extent[i] << " ";
        // os << '\n';
        std::cout << os.str();
        return grid;
    }

    static VTKExtent GetPExtentGlobal(const Representation& grep) {
        using Index = typename Representation::IndexGlobal;
        Index res = grep.GetGlobalResolution();
        dare::utils::Vector<3, int> res_3d(0, 0, 0);
        // VTK always expects 3D data, so we need accomodate that
        for (std::size_t d{0}; d < Dim; d++) {
            res_3d[d] = static_cast<int>(res[d]);
        }
        VTKExtent extent;
        // just to be consistent with the other formulations
        extent[0] = extent[1] = 0;
        extent[2] = extent[3] = 0;
        extent[4] = extent[5] = 0;
        // in x-direction
        extent[1] += res_3d.i();
        // extent[5] += res_3d.i();
        // in y-direction
        if constexpr (Dim > 1)
            extent[3] += res_3d.j();
        // in z-direction
        if constexpr (Dim > 2)
            // extent[1] += res_3d.k();
            extent[5] += res_3d.k();

        return extent;
    }

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
        // extent[0] = extent[1] = offset_cells[2];
        extent[0] = extent[1] = offset_cells[0];
        extent[2] = extent[3] = offset_cells[1];
        // extent[4] = extent[5] = offset_cells[0];
        extent[4] = extent[5] = offset_cells[2];
        // in x-direction
        // extent[5] += res_3d.i();
        extent[1] += res_3d.i();
        // in y-direction
        if constexpr (Dim > 1)
            extent[3] += res_3d.j();
        // in z-direction
        if constexpr (Dim > 2)
            // extent[1] += res_3d.k();
            extent[5] += res_3d.k();

        // for the parallel output, we need to remove the halo-layers
        // int num_ghost = static_cast<int>(grep.GetNumberGhostCells());
        // uint8_t bc = grep.GetBoundaryID();
        // if (!(bc & CartesianGrid::BOUNDARIES_WEST)) {
        //     extent[4] += num_ghost;
        // }
        // if (!(bc & CartesianGrid::BOUNDARIES_EAST)) {
        //     extent[5] -= num_ghost;
        // }
        // if constexpr (Dim > 1) {
        //     if (!(bc & CartesianGrid::BOUNDARIES_SOUTH)) {
        //         extent[2] += num_ghost;
        //     }
        //     if (!(bc & CartesianGrid::BOUNDARIES_NORTH)) {
        //         extent[3] -= num_ghost;
        //     }
        // }
        // if constexpr (Dim > 2) {
        //     if (!(bc & CartesianGrid::BOUNDARIES_BOTTOM)) {
        //         extent[0] += num_ghost;
        //     }
        //     if (!(bc & CartesianGrid::BOUNDARIES_TOP)) {
        //         extent[1] -= num_ghost;
        //     }
        // }
        return extent;
    }
};

}  // end namespace dare::io

#endif  // IO_VTKOPTIONS_CARTESIAN_H_
