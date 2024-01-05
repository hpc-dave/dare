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

#include <gtest/gtest.h>

#include "../../Utilities/Vector.h"
#include "../GridVector.h"

namespace dare::Data::test {

/*!
 * @brief small mock grid with minimal interface for the GridVector type
 */
class TestGrid {
public:
    using LocalOrdinalType = int32_t;
    using GlobalOrdinalType = int64_t;
    using LO = LocalOrdinalType;
    using Index = dare::utils::Vector<3, LocalOrdinalType>;
    using IndexGlobal = dare::utils::Vector<3, GlobalOrdinalType>;

    /*!
     * @brief represenation of mock grid
     */
    class TestGridRepresentation {
    public:
        using LO = TestGrid::LocalOrdinalType;

        /*!
         * @brief default constructor
         */
        TestGridRepresentation() : TestGridRepresentation(nullptr) {}

        /*!
         * @brief initializing constructor
         * @param _g pointer to grid
         */
        explicit TestGridRepresentation(TestGrid* _g) : grid(_g) {}

        /*!
         * @brief returns grid size
         */
        LO GetNumberLocalCells() const {
            LO n{1};
            for (auto e : grid->size)
                n *= e;
            return n;
        }

        /*!
         * @brief small mapping (Cartesian)
         * @param ind index for input
         */
        LO MapIndexToOrdinalLocal(const Index& ind) {
            LO n{0};
            for (std::size_t i{0}; i < ind.size(); i++) {
                LO n_sub{1};
                for (std::size_t j{i + 1}; j < grid->size.size(); j++)
                    n_sub *= grid->size[j];
                n += n_sub * ind[i];
            }
            return n;
        }

    private:
        TestGrid* grid;  //!< pointer to the test grid
    };
    using Representation = TestGridRepresentation;

    /*!
     * @brief default constructor
     */
    TestGrid() {}

    /*!
     * @brief default destructor
     */
    ~TestGrid() {}

    /*!
     * @brief Get simple representation
     */
    Representation GetRepresentation() {
        return Representation(this);
    }

    /*!
     * @brief set the size of the grid
     * @param _s size
     */
    void SetSize(Index _s) { size = _s; }

    Index size;  //!< size of the grid
};

}  // namespace dare::Data::test

/*!
 * @brief Testing suite for GridVector
 */
class GridVectorTest : public testing::Test {
protected:
    using TestGrid = dare::Data::test::TestGrid;
    using LO = typename TestGrid::LocalOrdinalType;
    using SC = double;
    static const std::size_t num_eq{5};
    using GridVector = dare::Data::GridVector<TestGrid, SC, num_eq>;

    void SetUp() override {
        resolution = TestGrid::Index(10, 12, 13);
        grid.SetSize(resolution);
        test_vector = GridVector("InitializationTest", grid.GetRepresentation());
    }

    typename TestGrid::Index resolution;
    TestGrid grid;
    GridVector test_vector;
};

TEST_F(GridVectorTest, Initialization) {
    std::size_t num_eq_comp = test_vector.GetNumEquations();
    LO num_expected_entries = num_eq_comp * grid.GetRepresentation().GetNumberLocalCells();
    EXPECT_EQ(test_vector.GetSize(), num_expected_entries);
}

TEST_F(GridVectorTest, GetSetValues) {
    for (std::size_t n{0}; n < test_vector.GetSize(); n++) {
        test_vector[n] = static_cast<SC>(n);
    }
    for (std::size_t n{0}; n < test_vector.GetSize(); n++) {
        EXPECT_EQ(test_vector[n], static_cast<SC>(n)) << "Error for access with local ordinal";
    }

    auto grid_rep = grid.GetRepresentation();
    for (LO i{0}; i < grid.size.i(); i++) {
        for (LO j{0}; j < grid.size.j(); j++) {
            for (LO k{0}; k < grid.size.k(); k++) {
                TestGrid::Index ind(i, j, k);
                for (std::size_t e{0}; e < test_vector.GetNumEquations(); e++) {
                    SC value_expect = grid_rep.MapIndexToOrdinalLocal(ind) * test_vector.GetNumEquations() + e;
                    SC value_vector = test_vector.At(ind, e);
                    EXPECT_EQ(value_expect, value_vector) << "Error for access with index";
                }
            }
        }
    }
}
