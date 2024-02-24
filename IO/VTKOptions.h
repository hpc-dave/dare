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

#ifndef IO_VTKOPTIONS_H_
#define IO_VTKOPTIONS_H_

#include <vtkNew.h>

#include "Grid/DefaultTypes.h"
#include "Utilities/Vector.h"
namespace dare::io {
using VTKExtent = dare::utils::Vector<6, int>;

/*!
 * @brief type of vtk output that is dealt with here
 */
enum class VTKDataAgglomerateType {
    VECTORS,
    SCALARS
};

enum class VTKOutputType {
    CELL_DATA,
    POINT_DATA
};

/*!
 * @brief VTKOptions for a grid, which allows mapping between data and vtk formats
 * @tparam Grid grid type
 */
template <typename Grid>
struct VTKOptions {
    // using GridType = void;  // will lead to compilation error if not set by SFINAE
    // using LO = typename Grid::LocalOrdinalType;

    /*!
     * @brief just some default which serves as example
     * @param grep instance of the grid
     * @param local_ordinal input ordinal
     * @return output ordinal according to the requirements of VTK
     * This may be required, since VTK orders structured grids differently than
     * this implementation
     */
    // static LO Map(const typename Grid::Representation& grep, LO local_ordinal) {
    //     return local_ordinal;
    // }

    /*!
     * @brief provides the grid for IO
     * @param grep representation of the grid
     */
    // static vtkNew<GridType> GetGrid(const typename Grid::Representation& grep) {
    //     return vtkNew<GridType>();
    // }
};

}  // end namespace dare::io

#endif  // IO_VTKOPTIONS_H_
