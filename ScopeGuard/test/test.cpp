
#include <mpi.h>
#include <string>
#include <gtest/gtest.h>
#include "../ScopeGuard.h"

TEST(ScopeGuardTest, IsRoot)
{
    int argc = 0;
    std::string args = "";
    char* argval = args.data();
    char** argv = &argval;
    dare::ScopeGuard scope(argc, argv);

    int rank{-1};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_root = rank == 0;
    ASSERT_EQ(is_root, scope.AmIRoot());
}