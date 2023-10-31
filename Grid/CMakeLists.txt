set(DaRe_HEADER_FILES
    "Cartesian.h"
    "CartesianDistribution.h"
    "CartesianRepresentation.h"
    "Cartesian.inl"
    "CartesianRepresentation.inl"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

target_sources(exec
    PRIVATE ${DaRe_HEADER_FILES}
)

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Grid")
