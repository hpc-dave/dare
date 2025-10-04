set(DaRe_HEADER_FILES
    "DDT.h"
    "DDT.inl"
    "GenericEquation.h"
    "GenericEquation.inl"
    "Operators.h"
    "FluxLimiter.h"
    "TimeDiscretizationSchemes.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

if(USE_GTEST)
    add_subdirectory(test)
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Equations")
