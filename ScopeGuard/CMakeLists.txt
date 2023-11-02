set(DaRe_HEADER_FILES
    "ScopeGuard.h"
)

set(DaRe_SOURCE_FILES
    "ScopeGuard.cpp"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

target_sources(exec
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

# Install header files
install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/ScopeGuard")
