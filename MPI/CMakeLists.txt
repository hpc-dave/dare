set(DaRe_HEADER_FILES
    "OBlackHoleStream.h"
    "MPITypeConverter.h"
    "ExecutionManager.h"
    "ExecutionManager.inl"
    "ExecutionManager_MPI.inl"
    "HaloBuffer.h"
    "HaloBuffer.inl"
    "SingleHaloBuffer.h"
    "SingleHaloBuffer.inl"
)

set(DaRe_SOURCE_FILES
    "ExecutionManager.cpp"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

target_sources(exec
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

if(USE_GTEST)
    target_sources(unit_test
        PRIVATE ${DaRe_HEADER_FILES}
        PRIVATE ${DaRe_SOURCE_FILES}
    )
    add_subdirectory(test)
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/MPI")
