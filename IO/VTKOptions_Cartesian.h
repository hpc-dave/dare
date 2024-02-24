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
    using Representation = typename Grid::Cartesian<Dim>::Representation;
    using GridType = vtkStructuredGrid;
    using LO = typename Grid::Cartesian<Dim>::LocalOrdinalType;

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

        std::swap(res_3d[0], res_3d[2]);
        std::swap(ind_3d[0], ind_3d[2]);

        dare::utils::Vector<3, LO> hsum(1, 1, 1);
        hsum[0] = res_3d[1] * res_3d[2];  // jmax * imax
        hsum[1] = res_3d[2];             // imax

        LO mapped_ordinal = hsum[0] * ind_3d[0] + hsum[1] * ind_3d[1] + hsum[2] * ind_3d[2];
        return mapped_ordinal;
    }

    /*!
     * @brief provides the grid for IO
     * @param grep representation of the grid
     */
    static vtkNew<GridType> GetGrid(const Representation& grep) {
        using Index = typename Representation::Index;
        Index res = grep.GetLocalResolution();
        dare::utils::Vector<3, int> res_3d(1, 1, 1);
        for (std::size_t d{0}; d < Dim; d++) {
            res_3d[d] = res[d] + 1;  // +1 for the number of points, required by VTK
        }
        // Here we compute some relevant data, including the cell dimensions and offset
        // Note, that we need to move the offset by half a cell size in each dimension
        // since the offset refers to the cell center!
        auto dn = grep.GetDistances();
        auto offset = grep.GetOffsetSize();
        offset -= 0.5 * dn;
        dare::utils::Vector<3, double> offset_3d(0., 0., 0.), dn_3d(0., 0., 0.);
        for (std::size_t d{0}; d < Dim; d++) {
            dn_3d[d] = dn[d];
            offset_3d[d] = offset[d];
        }

        vtkNew<vtkPoints> points;
        vtkNew<GridType> grid;
        points->SetNumberOfPoints(res_3d.i() * res_3d.j() * res_3d.k());
        grid->SetDimensions(res_3d.i(), res_3d.j(), res_3d.k());

        int count = 0;
        for (int k = 0; k < res_3d.k(); k++) {
            for (int j = 0; j < res_3d.j(); j++) {
                for (int i = 0; i < res_3d.i(); i++) {
                    double tuple[3];
                    tuple[0] = dn_3d.x() * i + offset_3d.x();
                    tuple[1] = dn_3d.y() * j + offset_3d.y();
                    tuple[2] = dn_3d.z() * k + offset_3d.z();
                    points->SetPoint(count, tuple);
                    count++;
                }
            }
        }
        grid->SetPoints(points);
        return grid;
    }
};

}  // end namespace dare::io

#endif  // IO_VTKOPTIONS_CARTESIAN_H_