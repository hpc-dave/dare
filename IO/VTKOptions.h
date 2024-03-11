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
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPStructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "Data/DefaultTypes.h"
#include "Utilities/Vector.h"
#include "VTKTypes.h"
#include "VTKPXMLStructuredGridWriter.h"

namespace dare::io {

/*!
 * @brief provides writer type according to vtkGrid via SFINAE
 * @tparam VTKGridType type of vtk grid
 */
template <typename VTKGridType>
struct VTKWriterMapper {
};

/*!
 * @brief partial specialization for vtkStructuredGrid
 */
template <>
struct VTKWriterMapper<vtkStructuredGrid> {
    using type = VTKPXMLStructuredGridWriter;
};

/*!
 * @brief partial specialization for vtkUnstructured grid
 */
template <>
struct VTKWriterMapper<vtkUnstructuredGrid> {
    using type = vtkXMLPUnstructuredGridWriter;
};

/*!
 * @brief user input to define if data is located at cell centers or nodes
 */
enum class VTKOutputType {
    CELL_DATA,      //!< data at the cell centers
    POINT_DATA      //!< data at the node points
};

/*!
 * @brief VTKOptions for a grid, which allows mapping between data and vtk formats
 * @tparam Grid grid type
 * \note if you end up here during compilation, then check your linking or write the specialization
 * of this struct
 */
template <typename Grid>
struct VTKOptions {
    static_assert(std::is_same_v<Grid, void>, "No specialization for your grid type found!");
};

}  // end namespace dare::io

#endif  // IO_VTKOPTIONS_H_
