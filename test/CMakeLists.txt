target_sources(unit_test
    PRIVATE main_test.cpp
    PRIVATE test_integration_Cartesian_Field.cpp
    PRIVATE test_integration_Cartesian_GridVector.cpp
    PRIVATE test_integration_Cartesian_MatrixBlock.cpp
    PRIVATE test_integration_Cartesian_Stencils.cpp
    PRIVATE test_integration_Cartesian_Trilinos.cpp
)
