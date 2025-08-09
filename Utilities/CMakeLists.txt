set(DaRe_HEADER_FILES
    "Array.h"
    "Array.inl"
    "CompileTimeFunctions.h"
    "Vector.h"
    "Vector.inl"
    "Vector_traits.h"
    "Hashes.h"
    "InitializationTracker.h"
    "Observer.h"
    "RandomNumberGenerator.h"
    "Errors.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

# Install header files
install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Utilities")
