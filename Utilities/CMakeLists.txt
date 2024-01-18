set(DaRe_HEADER_FILES
    "Vector.h"
    "Vector.inl"
    "Vector_traits.h"
    "Hashes.h"
    "InitializationTracker.h"
    "Errors.h"
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
endif()

# Install header files
install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Utilities")
