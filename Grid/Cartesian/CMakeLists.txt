set(DaRe_HEADER_FILES
    "CartesianMesh.h"
    "CartesianDistribution.h"
    "CartesianRepresentation.h"
    "CartesianMesh.inl"
    "CartesianRepresentation.inl"
    "Divergence_Cartesian.h"
    "Divergence_Cartesian.inl"
    "Gradient_Cartesian.h"
    "Gradient_Cartesian.inl"
    "Interpolation_Cartesian.h"
    "MatrixBlock_Cartesian.h"
    "MatrixBlock_Cartesian.inl"
    "Operators_Cartesian.h"
    "Stencils_Cartesian.h"
    "Stencils_Cartesian.inl"
    "TVD_Cartesian.h"
    "TVD_Cartesian.inl"
    "VTKOptions_Cartesian.h"
)

set(DaRe_SOURCE_FILES
    "Cartesian.cpp"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

target_sources(exec
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

if(USE_GTEST)
    target_sources(unit_test
        PRIVATE ${DaRe_HEADER_FILES}
        PRIVATE ${DaRe_SOURCE_FILES}
    )
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Grid/Cartesian")
