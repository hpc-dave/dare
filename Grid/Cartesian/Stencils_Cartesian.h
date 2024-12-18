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
#ifndef GRID_CARTESIAN_STENCILS_CARTESIAN_H_
#define GRID_CARTESIAN_STENCILS_CARTESIAN_H_

#include <string>
#include "Data/Stencil.h"
#include "Utilities/Vector.h"
#include "Grid/Cartesian.h"

namespace dare::Data {

/*!
 * @brief stores stencil based on cell centers
 * @tparam Dim dimension of grid
 * @tparam LO local ordinal type
 * @tparam GO global ordinal type
 * @tparam SC type of scalar
 */
template <std::size_t Dim, typename SC, std::size_t N>
class CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N> {
public:
    using GridType = dare::Grid::Cartesian<Dim>;                            //!< type of grid
    using ScalarType = SC;                                                  //!< scalar type
    static const std::size_t NUM_ENTRIES{GridType::STENCIL_SIZE};           //!< stencil size
    static const std::size_t NUM_COMPONENTS{N};                             //!< number of components in stencil
    using Positions = typename GridType::NeighborID;                        //!< convenient position definion
    using ComponentArray = dare::utils::Vector<NUM_ENTRIES, SC>;            //!< stencil of each component
    using DataArray = dare::utils::Vector<NUM_COMPONENTS, ComponentArray>;  //!< data storage for all entries
    using RHSType = dare::utils::Vector<NUM_COMPONENTS, SC>;                //!< storage of explicit rhs values

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
    CenterMatrixStencil(const CenterMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief copy assignment operator
     * @param other instance to copy from
     */
    CenterMatrixStencil<GridType, SC, N>&
    operator=(const CenterMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief multiplication assignment operator for scalar
     * @param v scalar value
     */
    CenterMatrixStencil<GridType, SC, N>& operator*=(SC v);

    /*!
     * @brief mulitplication operator for scalars
     * @param v scalar to multiply with
     */
    CenterMatrixStencil<GridType, SC, N> operator*(SC v) const;

    /*!
     * @brief division assignment operator for scalars
     * @param v scalar to divide by
     */
    CenterMatrixStencil<GridType, SC, N>& operator/=(SC v);

    /*!
     * @brief division operator for scalars
     * @param v scalar to divide by
     */
    CenterMatrixStencil<GridType, SC, N> operator/(SC v) const;

    /*!
     * @brief addition assignment of other stencil
     * @param other stencil to add from
     */
    CenterMatrixStencil<GridType, SC, N>& operator+=(const CenterMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief addition operator for another stencil
     * @param other stencil to add from
     */
    CenterMatrixStencil<GridType, SC, N> operator+(const CenterMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief subtraction assignment of other stencil
     * @param other stencil to subtract from this instance
     */
    CenterMatrixStencil<GridType, SC, N>& operator-=(const CenterMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief subtraction operator for another stencil
     * @param other stencil to subtract from this instance
     */
    CenterMatrixStencil<GridType, SC, N> operator-(const CenterMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief multiplication of this instance with another stencil
     * @param other stencil to mulitply with
     * \note the rhs will be added, not multiplied!
     */
    CenterMatrixStencil<GridType, SC, N>& operator*=(const CenterMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief multiplication with another stencil
     * @param other factors to mulitply with
     * \note the rhs will be added, not multiplied!
     */
    CenterMatrixStencil<GridType, SC, N> operator*(const CenterMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief division of this instance by another stencil
     * @param other stencil to divide by
     * \note the rhs will be added, not multiplied!
     */
    CenterMatrixStencil<GridType, SC, N>& operator/=(const CenterMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief division of stencil by another stencil
     * @param other stencil to divide by
     * \note the rhs will be added, not multiplied!
     */
    CenterMatrixStencil<GridType, SC, N> operator/(const CenterMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief returns center value
     * @param n component ID
     */
    SC& Center(std::size_t n);

    /*!
     * @brief returns center value
     * @param n component ID
     */
    SC Center(std::size_t n) const;

    /*!
     * @brief sets a specific value in the stencil
     * @param pos position of the value (e.g. CENTER)
     * @param n component ID
     * @param v value
     */
    void SetValue(Positions pos, std::size_t n, SC v);

    /*!
     * @brief sets value at the rhs
     * @param n component ID
     * @param v value to set
     */
    void SetRHS(std::size_t n, SC v);

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
     * @brief returns values of all components
     * @param pos position of value (e.g. CENTER)
     */
    dare::utils::Vector<N, SC> GetValues(Positions pos) const;

    /*!
     * @brief returns raw data array
     */
    DataArray& GetData();

    /*!
     * @brief return copy of raw data
     * @return 
     */
    const DataArray& GetData() const;

    /*!
     * @brief getter for right hand side
     */
    RHSType& GetRHS();

    /*!
     * @brief const getter for right hand side
     */
    const RHSType& GetRHS() const;

    /*!
     * @brief getter for right hand side
     * @param n component ID
     */
    SC& GetRHS(std::size_t n);

    /*!
     * @brief const getter for right hand side
     * @param n component ID
     */
    SC GetRHS(std::size_t n) const;

    CenterMatrixStencil<GridType, SC, 1> GetSlice(std::size_t n) const;

private:
    /*!
     * @brief internal range check, disabled with DARE_NDEBUG
     * @param func function name
     * @param pos position (e.g. center)
     * @param n component ID
     */
    void RangeCheck(std::string func, Positions pos, std::size_t n) const;
    DataArray coefficients;  //!< raw data array with stencil data
    RHSType rhs;             //!< potential array for storing RHS values
};

/*!
 * @brief For the Cartesian grid, no difference needs to be made for the matrix and value stencil
 * @tparam LO local ordinal
 * @tparam GO global orindal
 * @tparam SC scalar
 * @tparam Dim dimension
 * @tparam N number of components
 */
template <std::size_t Dim, typename SC, std::size_t N>
class CenterValueStencil<dare::Grid::Cartesian<Dim>, SC, N>
    : public CenterMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N> {
public:
    using GridType = dare::Grid::Cartesian<Dim>;                 //!< grid type
    using ScalarType = SC;                                       //!< scalar type
    using MatrixStencil = CenterMatrixStencil<GridType, SC, N>;  //!< matrix stencil type
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


/*!
 * @brief storage structure for matrix stencils referred to with faces of a stencil
 * @tparam SC type of scalar value
 * @tparam Dim dimension of the Cartesian grid
 * @tparam N number of components
 * The matrix entries associated witha face consist of two values. One for the neighbor
 * and one for the center. Therefore a more complex data structure is used here!
 */
template <std::size_t Dim, typename SC, std::size_t N>
class FaceMatrixStencil<dare::Grid::Cartesian<Dim>, SC, N> {
public:
    using GridType = dare::Grid::Cartesian<Dim>;                            //!< type of grid
    using ScalarType = SC;                                                  //!< scalar type
    static const std::size_t NUM_FACES{GridType::NUM_FACES};                //!< faces in stencil
    static const std::size_t NUM_COMPONENTS{N};                             //!< number of components in stencil
    using Positions = typename GridType::NeighborID;                        //!< convenient position definion
    using ComponentArray = dare::utils::Vector<NUM_FACES, SC>;              //!< stencil of each component
    using DataArray = dare::utils::Vector<NUM_COMPONENTS, ComponentArray>;  //!< data storage for all entries
    using RHSType = dare::Data::FaceValueStencil<GridType, SC, N>;          //!< array for explicit rhs values

    /*!
     * @brief default constructor
     */
    FaceMatrixStencil();

    /*!
     * @brief default destructor
     */
    ~FaceMatrixStencil();

    /*!
     * @brief copy constructor
     * @param other instance to copy from
     */
    FaceMatrixStencil(const FaceMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief copy assignment operator
     * @param other instance to copy from
     */
    FaceMatrixStencil<GridType, SC, N>&
    operator=(const FaceMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief multiplication assignment operator for scalar
     * @param v scalar value
     */
    FaceMatrixStencil<GridType, SC, N>& operator*=(SC v);

    /*!
     * @brief mulitplication operator for scalars
     * @param v scalar to multiply with
     */
    FaceMatrixStencil<GridType, SC, N> operator*(SC v) const;

    /*!
     * @brief division assignment operator for scalars
     * @param v scalar to divide by
     */
    FaceMatrixStencil<GridType, SC, N>& operator/=(SC v);

    /*!
     * @brief division operator for scalars
     * @param v scalar to divide by
     */
    FaceMatrixStencil<GridType, SC, N> operator/(SC v) const;

    /*!
     * @brief addition assignment of other stencil
     * @param other stencil to add from
     */
    FaceMatrixStencil<GridType, SC, N>& operator+=(const FaceMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief addition operator for another stencil
     * @param other stencil to add from
     */
    FaceMatrixStencil<GridType, SC, N> operator+(const FaceMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief subtraction assignment of other stencil
     * @param other stencil to subtract from this instance
     */
    FaceMatrixStencil<GridType, SC, N>& operator-=(const FaceMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief subtraction operator for another stencil
     * @param other stencil to subtract from this instance
     */
    FaceMatrixStencil<GridType, SC, N> operator-(const FaceMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief multiplication of this instance with another stencil
     * @param other stencil to mulitply with
     */
    FaceMatrixStencil<GridType, SC, N>& operator*=(const FaceMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief multiplication with another stencil
     * @param other factors to mulitply with
     */
    FaceMatrixStencil<GridType, SC, N> operator*(const FaceMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief division of this instance by another stencil
     * @param other stencil to divide by
     */
    FaceMatrixStencil<GridType, SC, N>& operator/=(const FaceMatrixStencil<GridType, SC, N>& other);

    /*!
     * @brief division of stencil by another stencil
     * @param other stencil to divide by
     */
    FaceMatrixStencil<GridType, SC, N> operator/(const FaceMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief sets a specific value in the stencil
     * @param pos position of the value (e.g. CENTER)
     * @param n component ID
     * @param v value
     */
    void SetValueNeighbor(Positions pos, std::size_t n, SC v);

    /*!
     * @brief Sets specific value for center coefficient
     * @param pos neighbor this contribution is associated with
     * @param n component ID
     * @param v value to set
     */
    void SetValueCenter(Positions pos, std::size_t n, SC v);

    /*!
     * @brief Sets specific value pair
     * @param pos neighbor this contribution is associated with
     * @param n component ID
     * @param v value to set
     */
    void SetValues(Positions pos, std::size_t n, SC v_nb, SC v_center);

    /*!
     * @brief sets value at the rhs
     * @param n component ID
     * @param v value to set
     */
    void SetRHS(Positions pos, std::size_t n, SC v);

    /*!
     * @brief sets a specific value in the stencil
     * @param pos position of the value (e.g. CENTER)
     * @param v vector with values
     */
    void SetValueNeighbor(Positions pos, const dare::utils::Vector<N, SC>& v);

    /*!
     * @brief Sets specific value for center coefficient
     * @param pos neighbor this contribution is associated with
     * @param v vector with values
     */
    void SetValueCenter(Positions pos, const dare::utils::Vector<N, SC>& v);

    /*!
     * @brief Sets specific value pair
     * @param pos neighbor this contribution is associated with
     * @param v vector with values
     */
    void SetValues(Positions pos,
                   const dare::utils::Vector<N, SC>& v_nb,
                   const dare::utils::Vector<N, SC>& v_center);

    /*!
     * @brief sets value at the rhs
     * @param v vector with values
     */
    void SetRHS(Positions pos, const dare::utils::Vector<N, SC>& v);

    /*!
     * @brief sets the whole stencil to a certain value
     * @param v value
     */
    void SetAll(SC v);

    /*!
     * @brief returns reference to value of face
     * @param pos position of value (e.g. WEST)
     * @param n component ID
     */
    SC& GetValueNeighbor(Positions pos, std::size_t n);

    /*!
     * @brief Getter for center associated coefficient
     * @param pos face of value
     * @param n component ID
     */
    SC& GetValueCenter(Positions pos, std::size_t n);

    /*!
     * @brief returns value
     * @param pos position of value (e.g. CENTER)
     * @param n component ID
     */
    SC GetValueNeighbor(Positions pos, std::size_t n) const;

    /*!
     * @brief returns value
     * @param pos position of value (e.g. CENTER)
     * @param n component ID
     */
    SC GetValueCenter(Positions pos, std::size_t n) const;

    /*!
     * @brief returns explicit contributions at face
     * @param pos position of value
     * @param n component ID
     */
    SC& GetRHS(Positions pos, std::size_t n);

    /*!
     * @brief constant getter for contributions at face
     * @param pos position of value
     * @param n component ID
     */
    SC GetRHS(Positions pos, std::size_t n) const;

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

    /*!
     * @brief getter for right hand side
     */
    RHSType& GetRHS();

    /*!
     * @brief const getter for right hand side
     */
    const RHSType& GetRHS() const;

    FaceMatrixStencil<GridType, SC, 1> GetSlice(std::size_t n) const;

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
    RHSType rhs;                  //!< storage for explicit components
};

/*!
 * @brief For the Cartesian grid, no difference needs to be made for the matrix and value stencil
 * @tparam LO local ordinal
 * @tparam GO global orindal
 * @tparam SC scalar
 * @tparam Dim dimension
 * @tparam N number of components
 */
template <std::size_t Dim, typename SC, std::size_t N>
class FaceValueStencil<dare::Grid::Cartesian<Dim>, SC, N> {
public:
    using GridType = dare::Grid::Cartesian<Dim>;                            //!< type of grid
    using ScalarType = SC;                                                  //!< scalar type
    static const std::size_t NUM_FACES{GridType::NUM_FACES};                //!< stencil size
    static const std::size_t NUM_COMPONENTS{N};                             //!< number of components in stencil
    using Positions = typename GridType::NeighborID;                        //!< convenient position definion
    using ComponentArray = dare::utils::Vector<NUM_FACES, SC>;              //!< stencil of each component
    using DataArray = dare::utils::Vector<NUM_COMPONENTS, ComponentArray>;  //!< data storage for all entries

    /*!
     * @brief default constructor
     */
    FaceValueStencil();

    /*!
     * @brief default destructor
     */
    ~FaceValueStencil();

    /*!
     * @brief copy assignment constructor
     * @param other
     */
    FaceValueStencil(const FaceValueStencil<GridType, SC, N>& other);

    /*!
     * @brief copy assignment operator
     * @param other instance to copy from
     */
    FaceValueStencil<GridType, SC, N>&
    operator=(const FaceValueStencil<GridType, SC, N>& other);

    /*!
     * @brief multiplication assignment operator for scalar
     * @param v scalar value
     */
    FaceValueStencil<GridType, SC, N>& operator*=(SC v);

    /*!
     * @brief mulitplication operator for scalars
     * @param v scalar to multiply with
     */
    FaceValueStencil<GridType, SC, N> operator*(SC v) const;

    /*!
     * @brief division assignment operator for scalars
     * @param v scalar to divide by
     */
    FaceValueStencil<GridType, SC, N>& operator/=(SC v);

    /*!
     * @brief division operator for scalars
     * @param v scalar to divide by
     */
    FaceValueStencil<GridType, SC, N> operator/(SC v) const;

    /*!
     * @brief addition assignment of other stencil
     * @param other stencil to add from
     */
    FaceValueStencil<GridType, SC, N>& operator+=(const FaceValueStencil<GridType, SC, N>& other);

    /*!
     * @brief addition operator for another stencil
     * @param other stencil to add from
     */
    FaceValueStencil<GridType, SC, N> operator+(const FaceValueStencil<GridType, SC, N>& other) const;

    /*!
     * @brief subtraction assignment of other stencil
     * @param other stencil to subtract from this instance
     */
    FaceValueStencil<GridType, SC, N>& operator-=(const FaceValueStencil<GridType, SC, N>& other);

    /*!
     * @brief subtraction operator for another stencil
     * @param other stencil to subtract from this instance
     */
    FaceValueStencil<GridType, SC, N> operator-(const FaceValueStencil<GridType, SC, N>& other) const;

    /*!
     * @brief multiplication of this instance with another stencil
     * @param other stencil to mulitply with
     */
    FaceValueStencil<GridType, SC, N>& operator*=(const FaceValueStencil<GridType, SC, N>& other);

    /*!
     * @brief multiplication with another stencil
     * @param other factors to mulitply with
     */
    FaceValueStencil<GridType, SC, N> operator*(const FaceValueStencil<GridType, SC, N>& other) const;

    /*!
     * @brief division of this instance by another stencil
     * @param other stencil to divide by
     */
    FaceValueStencil<GridType, SC, N>& operator/=(const FaceValueStencil<GridType, SC, N>& other);

    /*!
     * @brief division of stencil by another stencil
     * @param other stencil to divide by
     */
    FaceValueStencil<GridType, SC, N> operator/(const FaceValueStencil<GridType, SC, N>& other) const;

    /*!
     * @brief division operation with matrix stencil
     * @param other stencil to divide by
     */
    FaceMatrixStencil<GridType, SC, N> operator*(const FaceMatrixStencil<GridType, SC, N>& other) const;

    /*!
     * @brief sets a specific value in the stencil
     * @param pos position of the value (e.g. CENTER)
     * @param n component ID
     * @param v value
     */
    void SetValue(Positions pos, std::size_t n, SC v);

    /*!
     * @brief sets values in the stencil at certain position
     * @param pos position of the value (e.g. CENTER)
     * @param v value
     */
    void SetValues(Positions pos, const dare::utils::Vector<N, SC>& v);

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

}  // end namespace dare::Data

#include "Stencils_Cartesian.inl"

#endif  // GRID_CARTESIAN_STENCILS_CARTESIAN_H_
