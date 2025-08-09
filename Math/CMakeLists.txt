set(DaRe_HEADER_FILES
    "Divisors.h"
    "Interpolation.h"
    "Pow.h"
    "Random.h"
)

target_sources(dare
    PRIVATE ${DaRe_HEADER_FILES}
)

if(USE_GTEST)
   add_subdirectory(test)
endif()

install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/Math")
