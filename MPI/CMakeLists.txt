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
