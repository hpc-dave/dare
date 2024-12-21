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

#ifndef GRID_CARTESIAN_MATRIXBLOCK_CARTESIAN_H_
#define GRID_CARTESIAN_MATRIXBLOCK_CARTESIAN_H_

#include <vector>
#include <utility>

#include "MatrixSystem/MatrixBlock_Default.h"
#include "Grid/Cartesian/CartesianMesh.h"
#include "Grid/Cartesian/Stencils_Cartesian.h"
#include "Utilities/Errors.h"
#include "Utilities/Array.h"

namespace dare::Matrix {

using CartesianNeighbor = dare::Grid::CartesianNeighbor;

enum class CartesianNeighborBitSet : int16_t {
    CENTER = 1 << static_cast<char>(CartesianNeighbor::CENTER),
    WEST = 1 << static_cast<char>(CartesianNeighbor::WEST),
    EAST = 1 << static_cast<char>(CartesianNeighbor::EAST),
    SOUTH = 1 << static_cast<char>(CartesianNeighbor::SOUTH),
    NORTH = 1 << static_cast<char>(CartesianNeighbor::NORTH),
    BOTTOM = 1 << static_cast<char>(CartesianNeighbor::BOTTOM),
    TOP = 1 << static_cast<char>(CartesianNeighbor::TOP),
    FOURD_LOW = 1 << static_cast<char>(CartesianNeighbor::FOURD_LOW),
    FOURD_UP = 1 << static_cast<char>(CartesianNeighbor::FOURD_UP)
};

template <CartesianNeighbor A, CartesianNeighbor B>
constexpr bool IsSame() { return A == B; }

template <std::size_t Dim, typename O, typename SC, std::size_t N>
class MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N> : public MatrixBlockBase<O, SC, N> {
public:
    static const std::size_t STENCIL_SIZE = 2 * Dim + 1;
    using GridType = dare::Grid::Cartesian<Dim>;
    using GridRepresentation = typename GridType::Representation;
    using LocalOrdinalType = typename GridType::LocalOrdinalType;
    using GlobalOrdinalType = typename GridType::GlobalOrdinalType;
    using OrdinalType = O;
    using GO = GlobalOrdinalType;
    using LO = LocalOrdinalType;
    using SelfType = MatrixBlock<GridType, O, SC, N>;
    using Index = typename GridType::GetIndexType<O>::type;
    using ScalarArray = dare::utils::Vector<STENCIL_SIZE, SC>;

    /*!
     * @brief default constructor
     */
    MatrixBlock();

    /*!
     * @brief initializing constructor
     * @param g_rep grid representation
     * @param node grid node of the matrix block
     * @param size_hint resizing
     */
    MatrixBlock(const GridRepresentation* g_rep,
                O node,
                const dare::utils::Vector<N, std::size_t>& size_hint);

    /*!
     * @brief initializing without size hint
     * @param g_rep reference to the grid
     * @param node grid node of the matrix block
     */
    MatrixBlock(const GridRepresentation* g_rep, O node);

    /*!
     * @brief copy constructor
     * @param other instance to copy from
     */
    MatrixBlock(const SelfType& other);  // NOLINT

    /*!
     * @brief destructor
     */
    virtual ~MatrixBlock();

    /*!
     * @brief copy operator
     * @param other instance to copy from
     */
    SelfType& operator=(SelfType other);

    // /*!
    //  * @brief swap operator for the copy and swap idiom
    //  * @param m1 first block
    //  * @param m2 second block
    //  */
    // template <std::size_t Dims, typename Os, typename SCs, std::size_t Ns>
    // friend void swap(MatrixBlock<dare::Grid::Cartesian<Dims>, Os, SCs, Ns>& m1,
    //                  MatrixBlock<dare::Grid::Cartesian<Dims>, Os, SCs, Ns>& m2);

    /*!
     * @brief initialize the matrix block
     * @param g_rep grid representation
     * @param Node grid node
     */
    void Initialize(const GridRepresentation* g_rep, O Node);

    /*!
     * @brief initializes
     * @param g_rep reference to the grid
     * @param Node node of the grid
     * @param size_hint allocates memory
     */
    void Initialize(const GridRepresentation* g_rep,
                    O Node,
                    const dare::utils::Vector<N, std::size_t>& size_hint);

    /*!
     * @brief identifier, if the ordinal is a global ordinal
     * @return 
     */
    static constexpr bool IsGlobal();

    /*!
     * @brief return indices
     */
    Index GetIndex() const;

    /*!
     * @brief return indices
     */
    Index GetIndexInternal() const;

    /*!
     * @brief returns representation
     */
    const GridRepresentation* GetRepresentation() const;

    /*!
     * @brief tests if all ordinals are local
     */
    bool IsStencilLocal() const;

    /*!
     * @brief access to value of neighbor
     * @tparam CNB neighbor ID
     * @param n component id
     */
    template <CartesianNeighbor CNB>
    SC& Get(std::size_t nr, std::size_t nc = N);

    /*!
     * @brief non-templated access to value of neighbor
     * @tparam CNB neighbor ID
     * @param n component id
     */
    SC& Get(std::size_t nr, std::size_t nc, CartesianNeighbor cnb);

    /*!
     * @brief const access to value of neighbor
     * @tparam CNB neighbor ID
     * @param n component id
     */
    template <CartesianNeighbor CNB>
    SC Get(std::size_t nr, std::size_t nc = N) const;

    /*!
     * @brief non-templated const access to value of neighbor
     * @tparam CNB neighbor ID
     * @param n component id
     */
    SC Get(std::size_t nr, std::size_t nc, CartesianNeighbor cnb) const;

    /*!
     * @brief removes a neighbor from the stencil
     * @tparam CNB cartesian neighbor
     * @param n component ID
     */
    template <CartesianNeighbor CNB>
    void Remove(std::size_t nc, std::size_t nr);

    /*!
     * @brief removes a neighbor from the stencil
     * @param n component ID
     * @param cnb cartesian neighbor ID
     */
    void Remove(std::size_t nr, std::size_t nc, CartesianNeighbor cnb);

    /*!
     * @brief Inquiry, if a value was set
     * @tparam CNB neighbor ID
     * @param n component ID
     */
    template <CartesianNeighbor CNB>
    bool IsSet(std::size_t nr, std::size_t nc) const;

    /*!
     * @brief moves intermediate values to final ordinal & coefficient array
     */
    void Finalize();

    /*!
     * @brief assigns values from a stencil
     * @param s center stencil with matrix values
     */
    SelfType& Set(const dare::Data::CenterMatrixStencil<GridType, SC, N>& s);

    /*!
     * @brief assigns values from a single component stencil to single component
     * @param n component ID to which the value is assigned
     * @param s center stencil with matrix values
     */
    SelfType& Set(std::size_t n, const dare::Data::CenterMatrixStencil<GridType, SC, 1>& s);

    /*!
     * @brief operator for assigning values from a stencil
     * @param s center stencil with matrix values
     */
    SelfType& operator=(const dare::Data::CenterMatrixStencil<GridType, SC, N>& s);

    /*!
     * @brief adds values from stencil
     * @param s center stencil with matrix values
     */
    SelfType& Add(const dare::Data::CenterMatrixStencil<GridType, SC, N>& s);

    /*!
     * @brief adds values from single component stencil to single component
     * @param component ID to which the values are added
     * @param s center stencil with matrix values
     */
    SelfType& Add(std::size_t n, const dare::Data::CenterMatrixStencil<GridType, SC, 1>& s);

    /*!
     * @brief operator for adding values from stencil
     * @param s center stencil with matrix values
     */
    SelfType& operator+=(const dare::Data::CenterMatrixStencil<GridType, SC, N>& s);

    /*!
     * @brief subtracts values of stencil
     * @param s center stencil with matrix values
     */
    SelfType& Subtract(const dare::Data::CenterMatrixStencil<GridType, SC, N>& s);

    /*!
     * @brief subtracts values of single component stencil from single component
     * @param n component ID from which the values will be subtracted
     * @param s center stencil with matrix values
     */
    SelfType& Subtract(std::size_t n, const dare::Data::CenterMatrixStencil<GridType, SC, 1>& s);

    /*!
     * @brief operator for subtracting values of stencil
     * @param s center stencil with matrix values
     */
    SelfType& operator-=(const dare::Data::CenterMatrixStencil<GridType, SC, N>& s);

    /*!
     * @brief returns global internal ordinal
     */
    GO GetGlobalOrdinal() const;

    /*!
     * @brief returns local internal ordinal
     */
    LO GetLocalOrdinal() const;

    template <typename OS>
    friend OS& operator<<(OS& os, const MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N>& m) {
        os << "GO: " << m.GetGlobalOrdinal() << "\tLO: " << m.GetLocalOrdinal() << '\n';
        for (std::size_t n{0}; n < N; n++) {
            if (N > 1)
                os << "[" << n << "] ";
            if (m.IsSet<CartesianNeighbor::WEST>(n))
                os << "\tW: " << m.Get(n, CartesianNeighbor::WEST);

            if (Dim > 1 && m.IsSet<CartesianNeighbor::SOUTH>(n))
                os << "\tS: " << m.Get(n, CartesianNeighbor::SOUTH);

            if (Dim > 2 && m.IsSet<CartesianNeighbor::BOTTOM>(n))
                os << "\tB: " << m.Get(n, CartesianNeighbor::BOTTOM);

            if (m.IsSet<CartesianNeighbor::CENTER>(n))
                os << "\tC: " << m.Get(n, CartesianNeighbor::CENTER);

            if (Dim > 1 && m.IsSet<CartesianNeighbor::NORTH>(n))
                os << "\tN: " << m.Get(n, CartesianNeighbor::NORTH);

            if (Dim > 2 && m.IsSet<CartesianNeighbor::TOP>(n))
                os << "\tT: " << m.Get(n, CartesianNeighbor::TOP);

            if (m.IsSet<CartesianNeighbor::EAST>(n))
                os << "\tE: " << m.Get(n, CartesianNeighbor::EAST);

            os << "\tRHS: " << m.GetRhs(n) << '\n';
        }
        return os;
    }

private:

    /*!
     * @brief Getter for raw neighbor array
     */
    dare::utils::Array<N, N, ScalarArray>& GetNeighbors();

    /*!
     * @brief const getter for raw neighbor array
     */
    const dare::utils::Array<N, N, ScalarArray>& GetNeighbors() const;

    /*!
     * @brief Getter for partial neighbor array of component nr
     * @param nr component id
     */
    dare::utils::Vector<N, ScalarArray>& GetNeighbors(std::size_t nr);

    /*!
     * @brief const getter for partial neighbor array of component nr
     * @param nr component id
     */
    const dare::utils::Vector<N, ScalarArray>& GetNeighbors(std::size_t nr) const;

    /*!
     * @brief Get the neighbor array for a component nr with respect to the component nc
     * @param nr component id associated with the row
     * @param nc component id associated with the column
     */
    ScalarArray& GetNeighbors(std::size_t nr, std::size_t nc);

    /*!
     * @brief const getter the neighbor array for a component nr with respect to the component nc
     * @param nr component id associated with the row
     * @param nc component id associated with the column
     */
    const ScalarArray& GetNeighbors(std::size_t nr, std::size_t nc) const;

    /*!
     * @brief Getter for the raw bitset array
     */
    dare::utils::Array<N, N, char>& GetNeighborBitSet();

    /*!
     * @brief const getter for the raw bitset array
     */
    const dare::utils::Array<N, N, char>& GetNeighborBitSet() const;

    /*!
     * @brief getter for the bitset array associated with the component id nr
     * @param nr component id associated with the row
     */
    dare::utils::Vector<N, char>& GetNeighborBitSet(std::size_t nr);

    /*!
     * @brief const getter for the bitset array associated with the component id nr
     * @param nr component id associated with the row
     */
    const dare::utils::Vector<N, char>& GetNeighborBitSet(std::size_t nr) const;

    /*!
     * @brief getter for the stencil-bitset associated with the row-relative and column-relative component
     * @param nr component id associated with the row
     * @param nc component id associated with the column
     */
    char& GetNeighborBitSet(std::size_t nr, std::size_t nc);

    /*!
     * @brief const getter for the stencil-bitset associated with the row-relative and column-relative component
     * @param nr component id associated with the row
     * @param nc component id associated with the column
     */
    char GetNeighborBitSet(std::size_t nr, std::size_t nc) const;

    const GridRepresentation* g_rep;                      //!< reference to grid representation
    dare::utils::Array<N, N, ScalarArray> neighbors;      //!< holds values referencing to neighbor coefficients
    dare::utils::Array<N, N, char> neighbor_set;          //!< identifiers, if the neighbors were set
    Index ind_internal;                                   //!< indices of internal grid
    Index ind_full;                                       //!< indices of grid including halo/ghost cells
};
}  // end namespace dare::Matrix

#include "MatrixBlock_Cartesian.inl"

#endif  // GRID_CARTESIAN_MATRIXBLOCK_CARTESIAN_H_
