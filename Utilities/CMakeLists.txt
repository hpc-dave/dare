set(DaRe_HEADER_FILES
    "Vector.h"
    "Vector.inl"
    "Vector_traits.h"
    "Hashes.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

# Install header files
install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Utilities")
