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

#ifndef IO_VTKPXMLSTRUCTUREDGRIDWRITER_H_
#define IO_VTKPXMLSTRUCTUREDGRIDWRITER_H_

#include <vtkStructuredGrid.h>
#include <vtkXMLPStructuredGridWriter.h>
#include <vtkInformation.h>
#include <vtkInformationIntegerKey.h>

#include <iostream>
#include <string>
#include <vector>

#include "MPI/ExecutionManager.h"
#include "Utilities/Vector.h"
#include "VTKTypes.h"

namespace dare::io {

/*!
 * \brief custom writer for the parallel output of the structured grid
 */
class VTKPXMLStructuredGridWriter : public vtkXMLPStructuredGridWriter {
public:
    /*!
     * @brief default construction by vtkMacro
     */
    static VTKPXMLStructuredGridWriter* New();

    /*!
     * @brief provide the local extent that you want to see in the root file
     * @param local_extent local extent to print to the file
     * @param exman execution manager
     * \note to developers: invoces communication!
     */
    void SetPPieceExtent(const VTKExtent& local_extent, dare::mpi::ExecutionManager* exman);

    /*!
     * @brief override of the WritePPieceAttributes function
     * @param index index of the pieces (effectively the rank)
     */
    void WritePPieceAttributes(int index);

private:
    std::vector<vtkOrdinal> extent_array;   //!< list of the extents of each piece
};

}  // end namespace dare::io

#endif  // IO_VTKPXMLSTRUCTUREDGRIDWRITER_H_
