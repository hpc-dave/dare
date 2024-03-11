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

#include "Data/DefaultTypes.h"
#include "Grid/Cartesian.h"
#include "IO/VTKWriter.h"
#include "IO/FileSystemManager.h"
#include "Math/Interpolation_Cartesian.h"
#include "Math/Random.h"

namespace dare::test {

template <std::size_t Dim>
dare::utils::Vector<Dim, defaults::GlobalOrdinalType> GetResolutionTestVTKWriterCart() {
    dare::utils::Vector<Dim, defaults::GlobalOrdinalType> res;
    for (std::size_t n{0}; n < Dim; n++)
        res[n] = 10 + n;
    return res;
}
template <std::size_t Dim>
dare::utils::Vector<Dim, defaults::ScalarType> GetSizeTestVTKWriterCart() {
    dare::utils::Vector<Dim, defaults::ScalarType> size;
    for (std::size_t n{0}; n < Dim; n++)
        size[n] = 1. + n;
    return size;
}

}  // namespace dare::test

template<std::size_t Dim>
class VTKWriterTestsCartesian : public testing::Test {
public:
    static const std::size_t N{Dim};
    using GridType = dare::Grid::Cartesian<Dim>;
    using LO = typename GridType::LocalOrdinalType;
    using GO = typename GridType::GlobalOrdinalType;
    using SC = typename GridType::ScalarType;
    using Index = typename GridType::Index;
    using GridVector = dare::Data::GridVector<GridType, SC, N>;
    using Options = typename GridType::Options;
    using VTKWriter = dare::io::VTKWriter<GridType>;

    void SetUp() {
        const LO num_ghost{2};
        std::string name_grid = "Cartesian";
        name_grid += std::to_string(Dim) + "D";
        grid = std::make_unique<GridType>(
                                          name_grid,
                                          &exec_man,
                                          dare::test::GetResolutionTestVTKWriterCart<Dim>(),
                                          dare::test::GetSizeTestVTKWriterCart<Dim>(),
                                          num_ghost);
        fsys_man = std::make_unique<dare::io::FileSystemManager>(&exec_man);
        fsys_man->SetOutputPath("cartesian_test");
        fsys_man->CheckWithUser(false);
        fsys_man->OverwriteFiles(true);
        fsys_man->ClearContents(false);
    }

    std::unique_ptr<GridType> grid;        //!< the grid
    dare::mpi::ExecutionManager exec_man;  //!< the execution manager
    std::unique_ptr<dare::io::FileSystemManager> fsys_man;  //!< a file system manager
};

using VTKWriterTestsCartesian1Dim = VTKWriterTestsCartesian<1>;
using VTKWriterTestsCartesian2Dim = VTKWriterTestsCartesian<2>;
using VTKWriterTestsCartesian3Dim = VTKWriterTestsCartesian<3>;
TEST_F(VTKWriterTestsCartesian1Dim, GridVectorOutputTest) {
    Options opt(0);     // not staggered option
    Options opt_s(1);   // staggered options
    double time = 0.2;
    int step = 100;

    dare::math::Randomizer dice(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * dice.Get(); };
    auto grep = grid->GetRepresentation(opt);
    auto grep_s = grid->GetRepresentation(opt_s);
    GridVector vector_data("3D_vector", grep);
    GridVector scalar_data("3D_scalar", grep_s);
    VTKWriter writer(&exec_man, time, step);
    for (LO i{0}; i < grep.GetNumberLocalCells(); i++) {
        for (std::size_t n{0}; n < N; n++) {
            vector_data.At(i, n) = GetRandValue();
        }
    }
    for (LO i{0}; i < grep_s.GetNumberLocalCells(); i++) {
        for (std::size_t n{0}; n < N; n++) {
            scalar_data.At(i, n) = GetRandValue();
        }
    }
    scalar_data.SetComponentName(0, "myscalar");

    bool res = writer.Write(*fsys_man,
                 &vector_data,
                 &scalar_data);
    EXPECT_TRUE(res);
}

TEST_F(VTKWriterTestsCartesian2Dim, GridVectorOutputTest) {
    Options opt(0, 0);    // not staggered option
    Options opt_s(1, 0);  // staggered options
    double time = 0.5;
    int step = 200;

    dare::math::Randomizer dice(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * dice.Get(); };
    auto grep = grid->GetRepresentation(opt);
    auto grep_s = grid->GetRepresentation(opt_s);
    GridVector vector_data("2D_vector", grep);
    GridVector scalar_data("2D_scalar", grep_s);
    VTKWriter writer(&exec_man, time, step);
    for (LO i{0}; i < grep.GetNumberLocalCells(); i++) {
        for (std::size_t n{0}; n < N; n++) {
            vector_data.At(i, n) = GetRandValue();
        }
    }
    for (LO i{0}; i < grep_s.GetNumberLocalCells(); i++) {
        for (std::size_t n{0}; n < N; n++) {
            scalar_data.At(i, n) = GetRandValue();
        }
    }
    scalar_data.SetComponentName(0, "myscalar_0");
    scalar_data.SetComponentName(1, "myscalar_1");

    bool res = writer.Write(*fsys_man,
                            &vector_data,
                            &scalar_data);
    EXPECT_TRUE(res);
}

TEST_F(VTKWriterTestsCartesian3Dim, GridVectorOutputTest) {
    Options opt(0, 0, 0);    // not staggered option
    Options opt_s(1, 0, 0);  // staggered options
    double time = 0.67;
    int step = 300;

    dare::math::Randomizer dice(-10000, 10000);
    auto GetRandValue = [&]() { return 1e-4 * dice.Get(); };
    auto grep = grid->GetRepresentation(opt);
    auto grep_s = grid->GetRepresentation(opt_s);
    GridVector vector_data("1D_vector", grep);
    GridVector scalar_data("1D_scalar", grep_s);
    VTKWriter writer(&exec_man, time, step);
    for (LO i{0}; i < grep.GetLocalResolution().i(); i++) {
        for (LO j{0}; j < grep.GetLocalResolution().j(); j++) {
            for (LO k{0}; k < grep.GetLocalResolution().k(); k++) {
                Index ind(i, j, k);
                vector_data.At(ind, 0) = i;
                vector_data.At(ind, 1) = j * 2;
                vector_data.At(ind, 2) = k * 3;
            }
        }
    }
    for (LO i{0}; i < grep_s.GetNumberLocalCells(); i++) {
        for (std::size_t n{0}; n < N; n++) {
            scalar_data.At(i, n) = GetRandValue();
        }
    }
    scalar_data.SetComponentName(0, "myscalar_0");
    scalar_data.SetComponentName(1, "myscalar_1");
    scalar_data.SetComponentName(2, "myscalar_2");

    bool res = writer.Write(*fsys_man,
                            &vector_data,
                            &scalar_data);
    EXPECT_TRUE(res);
}
