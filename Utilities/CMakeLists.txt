set(DaRe_HEADER_FILES
    "Vector.h"
    "Vector.inl"
    "Hashes.h"
)

target_sources(dare
    PRIVATE ${YAPoNeMo_HEADER_FILES}
)
