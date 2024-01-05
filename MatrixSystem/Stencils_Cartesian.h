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
template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
class CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N> {
public:
    static const std::size_t NUM_ENTRIES{Dim * 2 + 1};                      //!< stencil size
    static const std::size_t NUM_COMPONENTS{N};                             //!< number of components in stencil
    using GridType = dare::Grid::Cartesian<Dim, LO, GO, SC>;                //!< type of grid
    using Positions = dare::Matrix::CartesianNeighbor;                      //!< convenient position definion
    using ComponentArray = dare::utils::Vector<NUM_ENTRIES, SC>;            //!< stencil of each component
    using DataArray = dare::utils::Vector<NUM_COMPONENTS, ComponentArray>;  //!< data storage for all entries

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
    CenterMatrixStencil(const CenterMatrixStencil<GridType, N>& other);

    /*!
     * @brief copy assignment operator
     * @param other instance to copy from
     */
    CenterMatrixStencil<GridType, N>&
    operator=(const CenterMatrixStencil<GridType, N>& other);

    /*!
     * @brief multiplication assignment operator for scalar
     * @param v scalar value
     */
    CenterMatrixStencil<GridType, N>& operator*=(SC v);

    /*!
     * @brief mulitplication operator for scalars
     * @param v scalar to multiply with
     */
    CenterMatrixStencil<GridType, N> operator*(SC v) const;

    /*!
     * @brief division assignment operator for scalars
     * @param v scalar to divide by
     */
    CenterMatrixStencil<GridType, N>& operator/=(SC v);

    /*!
     * @brief division operator for scalars
     * @param v scalar to divide by
     */
    CenterMatrixStencil<GridType, N> operator/(SC v) const;

    /*!
     * @brief addition assignment of other stencil
     * @param other stencil to add from
     */
    CenterMatrixStencil<GridType, N>& operator+=(const CenterMatrixStencil<GridType, N>& other);

    /*!
     * @brief addition operator for another stencil
     * @param other stencil to add from
     */
    CenterMatrixStencil<GridType, N> operator+(const CenterMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief subtraction assignment of other stencil
     * @param other stencil to subtract from this instance
     */
    CenterMatrixStencil<GridType, N>& operator-=(const CenterMatrixStencil<GridType, N>& other);

    /*!
     * @brief subtraction operator for another stencil
     * @param other stencil to subtract from this instance
     */
    CenterMatrixStencil<GridType, N> operator-(const CenterMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief multiplication of this instance with another stencil
     * @param other stencil to mulitply with
     */
    CenterMatrixStencil<GridType, N>& operator*=(const CenterMatrixStencil<GridType, N>& other);

    /*!
     * @brief multiplication with another stencil
     * @param other factors to mulitply with
     */
    CenterMatrixStencil<GridType, N> operator*(const CenterMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief division of this instance by another stencil
     * @param other stencil to divide by
     */
    CenterMatrixStencil<GridType, N>& operator/=(const CenterMatrixStencil<GridType, N>& other);

    /*!
     * @brief division of stencil by another stencil
     * @param other stencil to divide by
     */
    CenterMatrixStencil<GridType, N> operator/(const CenterMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief sets a specific value in the stencil
     * @param pos position of the value (e.g. CENTER)
     * @param n component ID
     * @param v value
     */
    void SetValue(Positions pos, std::size_t n, SC v);

    /*!
     * @brief sets the whole stencil to a certain value
     * @param v value
     */
    void SetAll(SC v);

    /*!
     * @brief returns reference to value
     * @param pos position of value (e.g. CENTER)
     * @param n component ID
     */
    SC& GetValue(Positions pos, std::size_t n);

    /*!
     * @brief returns value
     * @param pos position of value (e.g. CENTER)
     * @param n component ID
     */
    SC GetValue(Positions pos, std::size_t n) const;

    /*!
     * @brief returns raw data array
     */
    DataArray& GetData();

    /*!
     * @brief return copy of raw data
     * @return 
     */
    const DataArray& GetData() const;

private:
    /*!
     * @brief internal range check, disabled with DARE_NDEBUG
     * @param func function name
     * @param pos position (e.g. center)
     * @param n component ID
     */
    void RangeCheck(std::string func, Positions pos, std::size_t n) const;
    DataArray coefficients;  //!< raw data array with stencil data
};

/*!
 * @brief For the Cartesian grid, no difference needs to be made for the matrix and value stencil
 * @tparam LO local ordinal
 * @tparam GO global orindal
 * @tparam SC scalar
 * @tparam Dim dimension
 * @tparam N number of components
 */
template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
class CenterValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
    : public CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N> {
public:
    using MatrixStencil = CenterMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>;

    /*!
     * @brief default constructor
     */
    CenterValueStencil() : MatrixStencil() {}

    /*!
     * @brief conversion constructor to avoid compilation issues
     * @param source instance of parent class
     */
    CenterValueStencil(const MatrixStencil& source) : MatrixStencil(source) {} // NOLINT

    /*!
     * @brief default destructor
     */
    virtual ~CenterValueStencil() {}

    /*!
     * @brief assignment operator in the case of a provided Matrix stencil
     * @param source matrix stencil
     */
    CenterValueStencil& operator=(const MatrixStencil& source) {
        if (this == &source)
            return *this;
        MatrixStencil::operator=(source);
        return *this;
    }
};


template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
class FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N> {
public:
    static const std::size_t NUM_FACES{Dim * 2};                            //!< faces in stencil
    static const std::size_t NUM_COMPONENTS{N};                             //!< number of components in stencil
    using GridType = dare::Grid::Cartesian<Dim, LO, GO, SC>;                //!< type of grid
    using Positions = dare::Matrix::CartesianNeighbor;                      //!< convenient position definion
    using ComponentArray = dare::utils::Vector<NUM_FACES, SC>;            //!< stencil of each component
    using DataArray = dare::utils::Vector<NUM_COMPONENTS, ComponentArray>;  //!< data storage for all entries

    FaceMatrixStencil();
    ~FaceMatrixStencil();

    FaceMatrixStencil(const FaceMatrixStencil<GridType, N>& other);

    FaceMatrixStencil<GridType, N>&
    operator=(const FaceMatrixStencil<GridType, N>& other);

    /*!
     * @brief multiplication assignment operator for scalar
     * @param v scalar value
     */
    FaceMatrixStencil<GridType, N>& operator*=(SC v);

    /*!
     * @brief mulitplication operator for scalars
     * @param v scalar to multiply with
     */
    FaceMatrixStencil<GridType, N> operator*(SC v) const;

    /*!
     * @brief division assignment operator for scalars
     * @param v scalar to divide by
     */
    FaceMatrixStencil<GridType, N>& operator/=(SC v);

    /*!
     * @brief division operator for scalars
     * @param v scalar to divide by
     */
    FaceMatrixStencil<GridType, N> operator/(SC v) const;

    /*!
     * @brief addition assignment of other stencil
     * @param other stencil to add from
     */
    FaceMatrixStencil<GridType, N>& operator+=(const FaceMatrixStencil<GridType, N>& other);

    /*!
     * @brief addition operator for another stencil
     * @param other stencil to add from
     */
    FaceMatrixStencil<GridType, N> operator+(const FaceMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief subtraction assignment of other stencil
     * @param other stencil to subtract from this instance
     */
    FaceMatrixStencil<GridType, N>& operator-=(const FaceMatrixStencil<GridType, N>& other);

    /*!
     * @brief subtraction operator for another stencil
     * @param other stencil to subtract from this instance
     */
    FaceMatrixStencil<GridType, N> operator-(const FaceMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief multiplication of this instance with another stencil
     * @param other stencil to mulitply with
     */
    FaceMatrixStencil<GridType, N>& operator*=(const FaceMatrixStencil<GridType, N>& other);

    /*!
     * @brief multiplication with another stencil
     * @param other factors to mulitply with
     */
    FaceMatrixStencil<GridType, N> operator*(const FaceMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief division of this instance by another stencil
     * @param other stencil to divide by
     */
    FaceMatrixStencil<GridType, N>& operator/=(const FaceMatrixStencil<GridType, N>& other);

    /*!
     * @brief division of stencil by another stencil
     * @param other stencil to divide by
     */
    FaceMatrixStencil<GridType, N> operator/(const FaceMatrixStencil<GridType, N>& other) const;

    /*!
     * @brief sets a specific value in the stencil
     * @param pos position of the value (e.g. CENTER)
     * @param n component ID
     * @param v value
     */
    void SetValueNeighbor(Positions pos, std::size_t n, SC v);
    void SetValueCenter(Positions pos, std::size_t n, SC v);
    void SetValues(Positions pos, std::size_t n, SC v_nb, SC v_center);

    /*!
     * @brief sets the whole stencil to a certain value
     * @param v value
     */
    void SetAll(SC v);

    /*!
     * @brief returns reference to value
     * @param pos position of value (e.g. CENTER)
     * @param n component ID
     */
    SC& GetValueNeighbor(Positions pos, std::size_t n);
    SC& GetValueCenter(Positions pos, std::size_t n);

    /*!
     * @brief returns value
     * @param pos position of value (e.g. CENTER)
     * @param n component ID
     */
    SC GetValueNeighbor(Positions pos, std::size_t n) const;
    SC GetValueCenter(Positions pos, std::size_t n) const;

    /*!
     * @brief returns raw data array
     */
    DataArray& GetDataNeighbor();
    DataArray& GetDataCenter();

    /*!
     * @brief return copy of raw data
     * @return
     */
    const DataArray& GetDataNeighbor() const;
    const DataArray& GetDataCenter() const;

private:
    /*!
     * @brief internal range check, disabled with DARE_NDEBUG
     * @param func function name
     * @param pos position (e.g. center)
     * @param n component ID
     */
    void RangeCheck(std::string func, Positions pos, std::size_t n) const;

    DataArray coefficients_nb;    //!< values associated for each neighbor
    DataArray coefficients_c;     //!< values associated with center
};

/*!
 * @brief For the Cartesian grid, no difference needs to be made for the matrix and value stencil
 * @tparam LO local ordinal
 * @tparam GO global orindal
 * @tparam SC scalar
 * @tparam Dim dimension
 * @tparam N number of components
 */
template <std::size_t Dim, typename LO, typename GO, typename SC, std::size_t N>
class FaceValueStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>
    : public FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N> {
public:
    using MatrixStencil = FaceMatrixStencil<dare::Grid::Cartesian<Dim, LO, GO, SC>, N>;

    /*!
     * @brief default constructor
     */
    FaceValueStencil() : MatrixStencil() {}

    /*!
     * @brief conversion constructor to avoid compilation issues
     * @param source instance of parent class
     */
    FaceValueStencil(const MatrixStencil& source) : MatrixStencil(source) {}  // NOLINT

    /*!
     * @brief default destructor
     */
    virtual ~FaceValueStencil() {}

    /*!
     * @brief assignment operator in the case of a provided Matrix stencil
     * @param source matrix stencil
     */
    FaceValueStencil& operator=(const MatrixStencil& source) {
        if (this == &source)
            return *this;
        MatrixStencil::operator=(source);
        return *this;
    }
};

}  // end namespace dare::Data

#include "Stencils_Cartesian.inl"

#endif  // MATRIXSYSTEM_STENCILS_CARTESIAN_H_
