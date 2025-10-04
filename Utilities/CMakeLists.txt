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
    "ParameterList.h"
    "RandomNumberGenerator.h"
    "Errors.h"
    "Timer.h"
)

set(DaRe_SOURCE_FILES
    "Timer.cpp"
)


target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

# Install header files
install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Utilities")
