set(DaRe_HEADER_FILES
    "Field.h"
    "Field.inl"
    "GridVector.h"
    "GridVector.inl"
    "Stencil.h"
    "DefaultTypes.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

if(USE_GTEST)
    add_subdirectory(test)
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Data")
