target_sources(unit_test
    PRIVATE test.cpp
)

# Install header files
install(FILES ${DaRe_HEADER_FILES} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}/test")
