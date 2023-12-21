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

#ifndef MATRIXSYSTEM_MATRIXBLOCK_CARTESIAN_H_
#define MATRIXSYSTEM_MATRIXBLOCK_CARTESIAN_H_

#include <vector>

#include "MatrixBlock_Default.h"
#include "../Grid/Cartesian.h"

namespace dare::Matrix {

enum class CartesianNeighbor : char {
    CENTER = 0,
    WEST = 1,
    EAST = 2,
    SOUTH = 3,
    NORTH = 4,
    BOTTOM = 5,
    TOP = 6
};
enum class CartesianNeighborBitSet : char {
    CENTER = 1 << static_cast<char>(CartesianNeighbor::CENTER),
    WEST = 1 << static_cast<char>(CartesianNeighbor::WEST),
    EAST = 1 << static_cast<char>(CartesianNeighbor::EAST),
    SOUTH = 1 << static_cast<char>(CartesianNeighbor::SOUTH),
    NORTH = 1 << static_cast<char>(CartesianNeighbor::NORTH),
    BOTTOM = 1 << static_cast<char>(CartesianNeighbor::BOTTOM),
    TOP = 1 << static_cast<char>(CartesianNeighbor::TOP)
};

template <CartesianNeighbor A, CartesianNeighbor B>
constexpr bool IsSame() { return A == B; }

template <std::size_t Dim, typename O, typename SC, std::size_t N>
class MatrixBlock<dare::Grid::Cartesian<Dim>, O, SC, N> : public MatrixBlockBase<O, SC, N> {
public:
    using GridType = dare::Grid::Cartesian<Dim>;
    using GridRepresentation = typename GridType::Representation;
    using LocalOrdinalType = typename GridType::LocalOrdinalType;
    using GlobalOrdinalType = typename GridType::GlobalOrdinalType;
    using SelfType = MatrixBlock<GridType, O, SC, N>;
    using HostSpace = typename MatrixBlockBase<O, SC, N>::HostSpace;
    using ExecutionSpace = typename MatrixBlockBase<O, SC, N>::ExecutionSpace;
    template <typename DView, typename Space>
    using ReturnTypeDV = typename MatrixBlockBase<O, SC, N>::ReturnTypeDV<DView, Space>;
    using Index = typename GridType::GetIndexType<O>::type;

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
    SelfType& operator=(const SelfType& other);

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

    template <CartesianNeighbor CNB>
    SC& Get(std::size_t n);

    SC& Get(std::size_t n, CartesianNeighbor cnb);

    template <CartesianNeighbor CNB>
    SC Get(std::size_t n) const;

    SC Get(std::size_t n, CartesianNeighbor cnb) const;

    template <CartesianNeighbor CNB>
    bool IsSet(std::size_t n) const;

    void Finalize();

private:
    template <typename Space>
    typename Kokkos::DualView<SC**>::t_host& GetNeighbors();

    template <typename Space>
    const typename Kokkos::DualView<SC**>::t_host& GetNeighbors() const;

    template <typename Space>
    std::vector<char>& GetNeighborBitSet();

    template <typename Space>
    const std::vector<char>& GetNeighborBitSet() const;
    // const typename Kokkos::DualView<char*>::t_host& GetNeighborBitSet() const;

    const GridRepresentation* g_rep;   //!< reference to grid representation
    Kokkos::DualView<SC**> neighbors;  //!< holds values referencing to neighbor coefficients
    std::vector<char> neighbor_set;  //!< identifiers, if the neighbors were set
    // Kokkos::DualView<char*> neighbor_set;  //!< identifiers, if the neighbors were set
    Index ind_internal;                //!< indices of internal grid
    Index ind_full;                    //!< indices of grid including halo/ghost cells
};
}  // end namespace dare::Matrix

#include "MatrixBlock_Cartesian.inl"

#endif  // MATRIXSYSTEM_MATRIXBLOCK_CARTESIAN_H_
