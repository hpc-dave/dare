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
#ifndef MATRIXSYSTEM_STENCILS_CARTESIAN_H_
#define MATRIXSYSTEM_STENCILS_CARTESIAN_H_

#include <string>
#include "../Data/Stencil.h"
#include "../Utilities/Vector.h"
#include "../Grid/Cartesian.h"
#include "MatrixBlock_Cartesian.h"

namespace dare::Data {

/*!
 * @brief stores stencil based on cell centers
 * @tparam Dim dimension of grid
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 * @tparam SC type of scalar
 */
template <std::size_t Dim, typename LO, typename GO, typename SC>
class CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>> {
public:
    static const std::size_t NUM_ENTRIES{Dim * 2 + 1};          //!< stencil size
    // static const std::size_t NUM_COMPONENTS{N};                 //!< number of components in stencil
    using GridType = dare::Grid::Cartesian<Dim, LO, GO, SC>;    //!< type of grid
    using Positions = dare::Matrix::CartesianNeighbor;          //!< convenient position definion
    using DataArray = dare::utils::Vector<NUM_ENTRIES, SC>;       //!< data storage
    // using ComponentArray = dare::utils::Vector<NUM_ENTRIES, SC>;  //!< stencil of each component

    /*!
     * @brief default constructor
     */
    CenterMatrixStencil();

    /*!
     * @brief default destructor
     */
    ~CenterMatrixStencil();

    /*!
     * @brief copy assignment constructor
     * @param other 
     */
    CenterMatrixStencil(const CenterMatrixStencil<GridType>& other);

    /*!
     * @brief copy assignment operator
     * @param other instance to copy from
     */
    CenterMatrixStencil<GridType>&
    operator=(const CenterMatrixStencil<GridType>& other);

    /*!
     * @brief multiplication assignment operator for scalar
     * @param v scalar value
     */
    CenterMatrixStencil<GridType>& operator*=(SC v);

    /*!
     * @brief mulitplication operator for scalars
     * @param v scalar to multiply with
     */
    CenterMatrixStencil<GridType> operator*(SC v) const;

    /*!
     * @brief division assignment operator for scalars
     * @param v scalar to divide by
     */
    CenterMatrixStencil<GridType>& operator/=(SC v);

    /*!
     * @brief division operator for scalars
     * @param v scalar to divide by
     */
    CenterMatrixStencil<GridType> operator/(SC v) const;

    /*!
     * @brief addition assignment of other stencil
     * @param other stencil to add from
     */
    CenterMatrixStencil<GridType>& operator+=(const CenterMatrixStencil<GridType>& other);

    /*!
     * @brief addition operator for another stencil
     * @param other stencil to add from
     */
    CenterMatrixStencil<GridType> operator+(const CenterMatrixStencil<GridType>& other) const;

    /*!
     * @brief subtraction assignment of other stencil
     * @param other stencil to subtract from this instance
     */
    CenterMatrixStencil<GridType>& operator-=(const CenterMatrixStencil<GridType>& other);

    /*!
     * @brief subtraction operator for another stencil
     * @param other stencil to subtract from this instance
     */
    CenterMatrixStencil<GridType> operator-(const CenterMatrixStencil<GridType>& other) const;

    CenterMatrixStencil<GridType>& operator*=(const CenterMatrixStencil<GridType>& other);

    CenterMatrixStencil<GridType> operator*(const CenterMatrixStencil<GridType>& other) const;

    CenterMatrixStencil<GridType>& operator/=(const CenterMatrixStencil<GridType>& other);

    CenterMatrixStencil<GridType> operator/(const CenterMatrixStencil<GridType>& other) const;

    void SetValue(Positions pos, SC v);
    void SetAll(SC v);
    SC& GetValue(Positions pos);
    SC GetValue(Positions pos) const;
    DataArray& GetData();
    const DataArray& GetData() const;

private:
    void RangeCheck(std::string func, Positions pos) const;
    DataArray coefficients;
};

/*!
 * @brief For the Cartesian grid, no difference needs to be made for the matrix and value stencil
 * @tparam LO local ordinal
 * @tparam GO global orindal
 * @tparam SC scalar
 * @tparam Dim dimension
 */
template <std::size_t Dim, typename LO, typename GO, typename SC>
class CenterValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>>
    : public CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>> {
};

// /*!
//  * @brief dummy class for SFINAE
//  * @tparam Grid type of grid
//  * @tparam O type of ordinal
//  * @tparam SC type of scalar
//  */
// template <std::size_t Dim, typename LO, typename GO, typename SC>
// class FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>> {
// public:
//     static const std::size_t NUM_FACES{Dim * 2};
//     using GridType = dare::Grid::Cartesian<Dim, LO, GO, SC>;

//     CenterMatrixStencil();
//     ~CenterMatrixStencil();

//     CenterMatrixStencil(const CenterMatrixStencil<GridType>& other);

//     CenterMatrixStencil<GridType, Dim, LO, GO, SC>
//     operator=(const CenterMatrixStencil<GridType>& other);

// private:
//     dare::utils::Vector<NUM_FACES, SC> entries_f;
//     dare::utils::Vector<NUM_FACES, SC> entries_c;
// };

// /*!
//  * @brief dummy class for SFINAE
//  * @tparam Grid type of grid
//  * @tparam O type of ordinal
//  * @tparam SC type of scalar
//  */
// template <std::size_t Dim, typename LO, typename GO, typename SC>
// class FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>> {
// public:
//     static const std::size_t NUM_FACES{Dim * 2};
//     using GridType = dare::Grid::Cartesian<Dim, LO, GO, SC>;

//     CenterMatrixStencil();
//     ~CenterMatrixStencil();

//     CenterMatrixStencil(const CenterMatrixStencil<GridType>& other);

//     CenterMatrixStencil<GridType, Dim, LO, GO, SC>
//     operator=(const CenterMatrixStencil<GridType>& other);

// private:
//     dare::utils::Vector<NUM_FACES, SC> entries;
// };

}  // end namespace dare::Data

#include "Stencils_Cartesian.inl"

#endif  // MATRIXSYSTEM_STENCILS_CARTESIAN_H_
