set(DaRe_HEADER_FILES
    "Operators.h"
    "Gradient_Cartesian.h"
    "Gradient_Cartesian.inl"
    "Divergence_Cartesian.h"
    "Divergence_Cartesian.inl"
    "FluxLimiter.h"
    "TVD_Cartesian.h"
    "TVD_Cartesian.inl"
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
