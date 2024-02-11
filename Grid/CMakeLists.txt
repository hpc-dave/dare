set(DaRe_HEADER_FILES
    "Cartesian.h"
    "CartesianDistribution.h"
    "CartesianRepresentation.h"
    "Cartesian.inl"
    "CartesianRepresentation.inl"
    "DefaultTypes.h"
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
    add_subdirectory(test)
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Grid")
