set(DaRe_HEADER_FILES
    "Cartesian.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

if(USE_GTEST)
    add_subdirectory(test)
endif()

add_subdirectory(Cartesian)

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Grid")
