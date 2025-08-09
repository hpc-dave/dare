set(DaRe_HEADER_FILES
    "Diffusion.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/AnalyticalSolutions")
