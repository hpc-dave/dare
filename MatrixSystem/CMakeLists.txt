set(DaRe_HEADER_FILES
    "BiCGStab2.h"
    "BiCGStab2.inl"
    "Trilinos.h"
    "Trilinos.inl"
    "TrilinosSolver.h"
    "TrilinosSolver.inl"
    "MatrixBlockBase.h"
    "MatrixBlockBase.inl"
    "MatrixBlock.h"
    "MatrixBlock_Default.h"
    "MatrixBlock_Default.inl"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

target_sources(exec
    PRIVATE ${DaRe_HEADER_FILES}
)

if(USE_GTEST)
    target_sources(unit_test
        PRIVATE ${DaRe_HEADER_FILES}
    )
    add_subdirectory(test)
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/MatrixSystem")
