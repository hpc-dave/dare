set(DaRe_HEADER_FILES
    "ProjectionMethod.h"
    "ProjectionMethod.inl"
    "ProjectionMethod_freefunc.h"
    "ProjectionMethod_Cartesian.h"
    "PM_Information.h"
    "PM_details.h"
    "PM_Continuity.h"
    "PM_Continuity.inl"
    "PMTerminalOutputDefault.h"
)

set(DaRe_SOURCE_FILES
    "PMTerminalOutputDefault.cpp"
    "ProjectionMethod.cpp"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

if(USE_GTEST)
   add_subdirectory(test)
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Algorithm/NavierStokes")
