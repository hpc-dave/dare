set(DaRe_HEADER_FILES
    "OBlackHoleStream.h"
    "MPITypeConverter.h"
)

target_sources(dare
    PRIVATE ${YAPoNeMo_HEADER_FILES}
)
