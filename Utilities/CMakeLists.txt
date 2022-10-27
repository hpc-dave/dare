set(DaRe_HEADER_FILES
    "Vector.h"
    "Vector.inl"
    "Vector_traits.h"
    "Hashes.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)
