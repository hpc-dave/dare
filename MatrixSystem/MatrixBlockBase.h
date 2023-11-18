/*
 * MIT License
 *
 * Copyright (c) 2023 David Rieder

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

#ifndef MATRIXSYSTEM_MATRIXBLOCKBASE_H_
#define MATRIXSYSTEM_MATRIXBLOCKBASE_H_

#include <Kokkos_Core.hpp>

#include "../Utilities/Vector.h"

namespace dare::Matrix {

/*!
 * @brief basic data structure for matrix blocks
 * @tparam O ordinal type
 * @tparam SC scalar type
 * @tparam N number of components
 */
template <typename O, typename SC, std::size_t N>
class MatrixBlockBase {
public:
    using OrdinalType = O;
    using ScalarType = SC;
    using OrdinalArray = Kokkos::View<O*>;
    using ScalarArray = Kokkos::View<SC*>;

    /*!
     * @brief default constructor
     */
    MatrixBlockBase();

    /*!
     * @brief initializing constructor
     * @param node grid cell ID
     * @param size_hint number of elements which will be allocated
     */
    MatrixBlockBase(const O& node, const dare::utils::Vector<N, std::size_t>& size_hint);

    /*!
     * @brief copy constructor
     * @param other instance to copy from
     */
    MatrixBlockBase(const MatrixBlockBase<O, SC, N>& other);

    /*!
     * @brief copy assignment operator
     * @param other instance to copy from
     * @return 
     */
    MatrixBlockBase<O, SC, N>& operator=(const MatrixBlockBase<O, SC, N>& other);

    /*!
     * @brief destructor
     */
    virtual ~MatrixBlockBase();

    /*!
     * @brief Initializes matrix block and precomputes required values
     * @tparam Field type of field
     * @param g_rep grid representation
     * @param field field with values to copy from
     * @param node grid cell ID
     */
    void Initialize(O node, const dare::utils::Vector<N, std::size_t>& size_hint);

    /*!
     * @brief Allocates memory according to size hint
     * @param hint size hint for each row
     */
    void ProvideSizeHint(const dare::utils::Vector<N, std::size_t>& hint);

    /*!
     * @brief returns the local row of component n
     * @param n ID of component
     */
    O GetRow(std::size_t n) const;

    /*!
     * @brief provides rhs value of component n
     * @param n ID of component
     */
    SC GetRhs(std::size_t n) const;

    /*!
     * @brief provides rhs value of component n
     * @param n ID of component
     */
    SC& GetRhs(std::size_t n);

    /*!
     * @brief Returns initial guess of component n
     * @param n ID of component
     */
    SC GetInitialGuess(std::size_t n) const;

    /*!
     * @brief Returns initial guess of component n
     * @param n ID of component
     */
    SC& GetInitialGuess(std::size_t n);

    /*!
     * @brief returns number of row entries of component n
     * @param n ID of component
     */
    std::size_t GetNumEntries(std::size_t n) const;

    /*!
     * @brief returns array with columns ordinal of the row
     * @param n ID of associated component
     */
    const OrdinalArray& GetColumnOrdinals(std::size_t n) const;

    /*!
     * @brief returns array with column values of the row
     * @param n ID of associated component
     */
    const ScalarArray& GetColumnValues(std::size_t n) const;

    /*!
     * @brief returns reference ordinal according to specified position in array
     * @param n component ID
     * @param pos position in array
     * \warning There is not range check, neither for n nor for pos!
     */
    O& GetOrdinalByPosition(std::size_t n, std::size_t pos);

    /*!
     * @brief returns ordinal according to specified position in array
     * @param n component ID
     * @param pos position in array
     * \warning There is not range check, neither for n nor for pos!
     */
    O GetOrdinalByPosition(std::size_t n, std::size_t pos) const;

    /*!
     * @brief getter for coefficients depending on position in array
     * @param n component ID
     * @param pos position in array
     * @return reference to value
     */
    SC& GetCoefficientByPosition(std::size_t n, std::size_t pos);

    /*!
     * @brief getter for coefficients depending on position in array
     * @param n component ID
     * @param pos position in array
     * @return value
     */
    SC GetCoefficientByPosition(std::size_t n, std::size_t pos) const;

    /*!
     * @brief returns coefficient depending on column ordinal
     * @param n component ID
     * @param ordinal column ordinal
     * @return reference to value
     * \note This will resize and insert a value, if ordinal was not found!
     * This function is comparatively heavy and potentially a performance
     * bottleneck if overused. Try to avoid as much as possible!
     */
    SC& GetCoefficientByOrdinal(std::size_t n, O ordinal);

    /*!
     * @brief returns coefficient depending on column ordinal
     * @param n component ID
     * @param ordinal column ordinal
     * @return reference to value
     * \note if the ordinal was not found, then the value at position 0 is returned
     */
    SC GetCoefficientByOrdinal(std::size_t n, O ordinal) const;

    /*!
     * @brief set rhs value of component n
     * @param n ID of component
     * @param value value of component
     */
    void SetRhs(std::size_t n, SC value);

    /*!
     * @brief set initial guess of component n
     * @param n ID of component
     * @param value value of component
     */
    void SetInitialGuess(std::size_t n, SC value);

    /*!
     * @brief sets all rhs values
     * @param rhs_array pointer to array with values
     * \note assumes that array is of size N and random access
     */
    void SetRhs(const SC* rhs_array);

    /*!
     * @brief sets all initial values
     * @param x_array pointer to array with values
     * \note assumes that array is of size N and random access
     */
    void SetInitialGuess(const SC* x_array);

    /*!
     * @brief set a coefficient
     * @param n_row row component id (<N)
     * @param id_col column id
     * @param value coefficient value
     */
    void SetCoefficient(std::size_t n_row, O id_col, SC value);

    /*!
     * @brief set multiple coefficients
     * @tparam Array1 array type of ordinals
     * @tparam Array2 array type of values
     * @param n_row component number (<N)
     * @param size number of elements in the arrays
     * @param id_col array of column ordinals
     * @param values array of coefficients
     * 
     * The arrays' sizes need to be greater or equal to the provided size
     * parameter. Additionally, the elements in the arrays need to be accessible
     * by array[i], where i is of type std::size_t. It is assumed, that the
     * ordinals and values are related with respect to their position in the
     * arrays.
     */
    template <typename Array1, typename Array2>
    void SetCoefficients(std::size_t n_row, std::size_t size, const Array1& id_col, const Array2& values);

    /*!
     * @brief removes coefficient from array
     * @param n component ID
     * @param pos position in array
     * \warning This function should be used with care! Potentially it's expensive!
     */
    void RemoveCoefficientByPosition(std::size_t n, std::size_t pos);

    /*!
     * @brief removes multiple coefficients
     * @tparam Array some kind of array which can be accessed by []
     * @param n component ID
     * @param positions array with positions
     * @param num_entries number of entries in the array
     */
    template <typename Array>
    void RemoveCoefficientsByPositions(std::size_t n, const Array& positions, std::size_t num_entries);

    /*!
     * @brief removes coefficient from array
     * @param n component ID
     * @param ordinal column ordinal associated with value
     * \warning This function should be used with care! Potentially it's expensive!
     */
    void RemoveCoefficientByOrdinal(std::size_t n, O ordinal);

    /*!
     * @brief removes multiple coefficients
     * @tparam Array some kind of array which can be accessed by []
     * @param n component ID
     * @param ordinals array with column ordinals
     * @param num_entries number of entries in the array
     */
    template <typename Array>
    void RemoveCoefficientsByOrdinals(std::size_t n, const Array& ordinals, std::size_t num_entries);

private:
    OrdinalArray ordinals[N];       //!< arrays with columns indices
    ScalarArray coefficients[N];    //!< arrays with column coefficients
    ScalarArray initial_guess;      //!< array with initial guess
    ScalarArray rhs;                //!< array with rhs
    OrdinalType node;                   //!< node associated with this matrix block
};
}  // namespace dare::Matrix

#include "MatrixBlockBase.inl"
#endif  // MATRIXSYSTEM_MATRIXBLOCKBASE_H_
