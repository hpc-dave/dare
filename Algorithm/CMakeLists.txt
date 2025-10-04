set(DaRe_HEADER_FILES
    "AlgorithmTraits.h"
    "ConstantTimeStep.h"
    "StaticInformation.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

add_subdirectory(NavierStokes)

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Algorithm")
