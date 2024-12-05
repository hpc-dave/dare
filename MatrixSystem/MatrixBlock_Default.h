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

#ifndef MATRIXSYSTEM_MATRIXBLOCK_DEFAULT_H_
#define MATRIXSYSTEM_MATRIXBLOCK_DEFAULT_H_

#include <type_traits>

#include "MatrixBlockBase.h"
#include "Utilities/Vector.h"

namespace dare::Matrix {

/*!
 * @brief default version of MatrixBlock
 * @tparam Grid type of grid
 * @tparam O type of ordinal
 * @tparam SC type of scalar
 * @tparam N number of components
 */
template<typename Grid, typename O, typename SC, std::size_t N>
class MatrixBlock : public MatrixBlockBase<O, SC, N> {
public:
    using GridRepresentation = typename Grid::Representation;
    using LocalOrdinalType = typename Grid::LocalOrdinalType;
    using GlobalOrdinalType = typename Grid::GlobalOrdinalType;
    using SelfType = MatrixBlock<Grid, O, SC, N>;
    using GO = GlobalOrdinalType;
    using LO = LocalOrdinalType;

    /*!
     * @brief default constructor
     */
    MatrixBlock();

    /*!
     * @brief initializing constructor
     * @param g_rep refernce to the grid
     * @param node grid node of the matrix block
     * @param size_hint number of matrix elements per component
     */
    MatrixBlock(const GridRepresentation* g_rep,
                O node,
                const dare::utils::Vector<N, std::size_t>& size_hint);

    /*!
     * @brief initialzing constructor without memory allocation
     * @param g_rep reference to the grid
     * @param node grid node of the matrix block
     */
    MatrixBlock(const GridRepresentation* g_rep, O node);

    /*!
     * @brief copy constructor
     * @param other instance to copy from
     */
    MatrixBlock(const SelfType& other);     // NOLINT

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
     * @brief Initializes the matrix block
     * @param g_rep reference to grid
     * @param Node grid node of matrix block
     */
    void Initialize(const GridRepresentation* g_rep, O Node);

    /*!
     * @brief initializes matrix block and allocated memory
     * @param g_rep reference to grid
     * @param Node grid node of matrix block
     * @param size_hint number of matrix elements for each component
     */
    void Initialize(const GridRepresentation* g_rep,
                    O Node,
                    const dare::utils::Vector<N, std::size_t>& size_hint);

    /*!
     * @brief checks, if ordinal is of global ordinal type
     */
    bool IsGlobal() const;

    /*!
     * @brief returns reference to grid
     */
    const GridRepresentation* GetRepresentation() const;

    /*!
     * @brief tests if all ordinals are local
     */
    bool IsStencilLocal() const;

    /*!
     * @brief provides local internal ordinal
     */
    LO GetLocalOrdinal() const;

    /*!
     * @brief provides global internal ordinal
      */
    GO GetGlobalOrdinal() const;

    /*!
     * @brief post assembly step
     * \note can be used for sorting and so on
     */
    void Finalize();

private:
    const GridRepresentation* g_rep;  //!< reference to grid representation
};

/*!
 * @brief Converts a local matrix block to a global one and vice-versa
 * @tparam Grid type of grid
 * @tparam Osource ordinal type of source
 * @tparam Otarget ordinal type of target
 * @tparam SC scalar type
 * @tparam N number of components
 * @param source input matrix block
 */
template <typename Otarget, typename Grid, typename Osource, typename SC, std::size_t N>
MatrixBlock<Grid, Otarget, SC, N> Convert(const MatrixBlock<Grid, Osource, SC, N>& source) {
    using LO = typename Grid::LocalOrdinalType;
    using GO = typename Grid::GlobalOrdinalType;
    static_assert(std::is_same_v<Otarget, LO> || std::is_same_v<Otarget, GO>,
                  "the target ordinal is not part of the grid");
    if constexpr (std::is_same_v<Osource, Otarget>) {
        return source;
    }
    const typename Grid::Representation* g_rep = source.GetRepresentation();
    dare::utils::Vector<N, std::size_t> size_hint;
    for (std::size_t n{0}; n < N; n++)
        size_hint[n] = source.GetNumEntries(n);
    if constexpr (std::is_same_v<Otarget, GO>) {
        MatrixBlock<Grid, Otarget, SC, N> target(g_rep, g_rep->MapLocalToGlobalInternal(source.GetNode()), size_hint);
        for (std::size_t n{0}; n < N; n++) {
            Kokkos::View<GO*> ordinals_global("global_ordinals", size_hint[n]);
            for (std::size_t p{0}; p < size_hint[n]; p++) {
                LO col_loc = source.GetOrdinalByPosition(n, p);
                GO col_glob = col_loc / N;
                GO component_id = col_loc - (col_glob * N);
                col_glob = g_rep->MapLocalToGlobalInternal(col_glob);
                col_glob = col_glob * N + component_id;
                ordinals_global[p] = col_glob;
            }
            target.SetCoefficients(n, size_hint[n], ordinals_global, source.GetColumnValues(n));
        }
        return target;
    } else {
        MatrixBlock<Grid, Otarget, SC, N> target(g_rep, g_rep->MapGlobalToLocalInternal(source.GetNode()), size_hint);
        for (std::size_t n{0}; n < N; n++) {
            Kokkos::View<GO*> ordinals_local("local_ordinals", size_hint[n]);
            for (std::size_t p{0}; p < size_hint[n]; p++) {
                GO col_glob = source.GetOrdinalByPosition(n, p);
                GO row_glob = col_glob / N;
                LO component_id = col_glob - (row_glob * N);
                LO row_loc = g_rep->MapGlobalToLocalInternal(row_glob);
                row_loc = row_loc * N + component_id;
                ordinals_local[p] = row_loc;
            }
            target.SetCoefficients(n, size_hint[n], ordinals_local, source.GetColumnValues(n));
        }
        MatrixBlock<Grid, Otarget, SC, N> target2 = target;
        return target;
    }
}
}  // end namespace dare::Matrix

#include "MatrixBlock_Default.inl"

#endif  // MATRIXSYSTEM_MATRIXBLOCK_DEFAULT_H_
