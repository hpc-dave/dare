set(DaRe_HEADER_FILES
    "OBlackHoleStream.h"
    "MPITypeConverter.h"
    "ExecutionManager.h"
    "ExecutionManager.inl"
    "ExecutionManager_MPI.inl"
)

set(DaRe_SOURCE_FILES
    "ExecutionManager.cpp"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
    PRIVATE ${DaRe_SOURCE_FILES}
)

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/MPI")
