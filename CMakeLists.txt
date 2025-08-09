cmake_minimum_required(VERSION 3.13)

project(DaRe
    VERSION 0.4.0
    DESCRIPTION "Dave's repository"
    LANGUAGES C CXX
)

# setting it here generates the compile_commands.json at first invokation
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(DEBUG
    "Enables debugging statements for the repository"
    OFF
)

option(USE_GTEST
    "Enables debugging statements for the repository"
    OFF
)

option(ENABLE_COVERAGE
    "Enables code coverage for the repository"
    OFF
)

include(GNUInstallDirs)
# a few settings

# define main target
add_executable(exec)
add_library(dare STATIC)
if(USE_GTEST)
    add_executable(unit_test)
endif()

# exports compile_command.json with relevant parths
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
set(CMAKE_FIND_PACKAGE_SORT_ORDER NATURAL)
set(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)

# set the properties of the target (might be moved to another file)

set_target_properties(dare
           PROPERTIES CXX_STANDARD 20
           VERSION 0.1
)

set_target_properties(exec
           PROPERTIES CXX_STANDARD 20
           CXX_STANDARD_REQUIRED TRUE
)

if(USE_GTEST)
    set_target_properties(unit_test
           PROPERTIES CXX_STANDARD 20
)
endif()

target_include_directories(dare PRIVATE .)
set(PACKAGE_NAME_EXPORT     "${CMAKE_PROJECT_NAME}Targets"                  )
# set(PACKAGE_CMAKEDIR_NAME   "${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}")
set(PACKAGE_CMAKEDIR_NAME   "${CMAKE_PROJECT_NAME}")
set(PACKAGE_BUILD_TREE_NAME "${CMAKE_PROJECT_NAME}-Exports"                 )
# set(PACKAGE_VERSION_DIR     "${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}")
set(PACKAGE_VERSION_DIR     "${CMAKE_PROJECT_NAME}")

# set compilation flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(dare
        PRIVATE -march=native
        PRIVATE -Wall
        PRIVATE -faligned-new
        PRIVATE -fmessage-length=0
    )
    target_compile_options(exec
        PRIVATE -march=native
        PRIVATE -Wall
        PRIVATE -faligned-new
        PRIVATE -fmessage-length=0
    )
    if(USE_GTEST)
        target_compile_options(unit_test
        PRIVATE -march=native
        PRIVATE -Wall
        PRIVATE -faligned-new
        PRIVATE -fmessage-length=0
        )
    endif()
    if(DEBUG)
        target_compile_options(dare
        PRIVATE -Werror
        )
        target_compile_options(exec
        PRIVATE -Werror
        )
        if(USE_GTEST)
            target_compile_options(unit_test
            PRIVATE -Werror
            )
        endif()
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
            if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
                target_compile_options(dare PRIVATE -Wno-maybe-uninitialized)
                target_compile_options(exec PRIVATE -Wno-maybe-uninitialized)
                if(USE_GTEST)
                    target_compile_options(unit_test PRIVATE -Wno-maybe-uninitialized)
                endif()
            endif()
        endif()
    endif()
else()
    message(WARNING "DaRe: Flags have not been implemented for the current compiler! (${CMAKE_CXX_COMPILER_ID})")
endif()

# locate packages
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set_target_properties(OpenMP::OpenMP_CXX
        PROPERTIES IMPORTED_GLOBAL TRUE
    )
endif()

find_package(MPI REQUIRED)
if(MPI_CXX_FOUND)
    set_target_properties(MPI::MPI_CXX
        PROPERTIES IMPORTED_GLOBAL TRUE
    )
endif()

find_package(VTK REQUIRED)
if(VTK_FOUND)
    set_target_properties(${VTK_LIBRARIES}
        PROPERTIES IMPORTED_GLOBAL TRUE
    )
endif()

find_package(Boost 1.74 REQUIRED
    COMPONENTS filesystem
               regex 
               serialization
               system
)

# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(Trilinos REQUIRED
             COMPONENTS Teuchos
                        Tpetra
                        Belos
                        MueLu
                        Kokkos
                        Ifpack2
)

if(Trilinos_FOUND)
    message(STATUS "Found Trilinos version ${Trilinos_VERSION} in ${Trilinos_DIR}")
    if(${Trilinos_VERSION} VERSION_LESS "14")
        message(ERROR "This project requires minimum Trilinos 14.0")
    endif()

    MESSAGE("\nFollowing Trilinos variables were determined:")
        MESSAGE("   Trilinos_DIR               = ${Trilinos_DIR}")
        MESSAGE("   Trilinos_VERSION           = ${Trilinos_VERSION}")
        MESSAGE("   Trilinos_PACKAGE_LIST      = ${Trilinos_PACKAGE_LIST}")
        MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
        MESSAGE("   Trilinos_INCLUDE_DIRS      = ${Trilinos_INCLUDE_DIRS}")
        MESSAGE("   Trilinos_LIBRARY_DIRS      = ${Trilinos_LIBRARY_DIRS}")
        MESSAGE("   Trilinos_MPI_LIBRARIES     = ${Trilinos_MPI_LIBRARIES}")
        MESSAGE("   Trilinos_CXX_COMPILER      = ${Trilinos_CXX_COMPILER}")
        MESSAGE("   Trilinos_CXX_COMPILER_FLAGS= ${Trilinos_CXX_COMPILER_FLAGS}")
        MESSAGE("   Trilinos_TPL_LIST          = ${Trilinos_TPL_LIST}")
        MESSAGE("   Trilinos_TPL_INCLUDE_DIRS  = ${Trilinos_TPL_INCLUDE_DIRS}")
        MESSAGE("   Trilinos_TPL_LIBRARY_DIRS  = ${Trilinos_TPL_LIBRARY_DIRS}")
        MESSAGE("   Trilinos_LIBRARIES         = ${Trilinos_LIBRARIES}")

# Interesting, yet massive output, skipping that for now
#        MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
        
        MESSAGE("\n")

    # add_library(trilinos INTERFACE)
    # set_target_properties(trilinos PROPERTIES
    #     INTERFACE_LINK_LIBRARIES "${Trilinos_LIBRARIES};${Trilinos_TPL_LIBRARIES};"
    #     INTERFACE_INCLUDE_DIRECTORIES "${Trilinos_INCLUDE_DIRS};${Trilinos_TPL_INCLUDE_DIRS};"
    # )
endif()

if(Boost_FOUND)
    set_target_properties(Boost::boost
                          Boost::filesystem
                          Boost::regex
                          Boost::serialization
                          Boost::system

        PROPERTIES IMPORTED_GLOBAL TRUE
    )
endif()

if(USE_GTEST)
    find_package(GTest REQUIRED)
    if(GTest_FOUND)
        set_target_properties(GTest::GTest
            PROPERTIES IMPORTED_GLOBAL TRUE
        )
    endif()
endif()

add_subdirectory(Algorithm                                            )
add_subdirectory(AnalyticalSolutions                                  )
add_subdirectory(Data                                                 )
add_subdirectory(Equations                                            )
add_subdirectory(FileFormats                                          )
add_subdirectory(Utilities                                            )
add_subdirectory(MPI                                                  )
add_subdirectory(Grid                                                 )
add_subdirectory(IO                                                   )
add_subdirectory(MatrixSystem                                         )
add_subdirectory(Math                                                 )
add_subdirectory(ScopeGuard                                           )
if(USE_GTEST)
    add_subdirectory(test)
endif()


if(${DEBUG})
    message("\nCompiling in DEBUG mode\n")
else()
    target_compile_definitions(dare
        PUBLIC -DDARE_NDEBUG
    )
    target_compile_definitions(exec
        PUBLIC -DDARE_NDEBUG
    )
    target_compile_definitions(unit_test
        PUBLIC -DDARE_NDEBUG
    )
endif()


if(${ENABLE_COVERAGE})
    # This function will configure the compiler flags for a target that will track the code coverage during program
    # execution.
    #
    # INPUT
    #   tgt: the target for which code coverage should be enabled
    #
    function(set_enable_coverage tgt)

        set(is_coverage_enabled OFF)
        # Verify if debug-mode is enabled. Code coverage tracking will show unexpected results if the code has been
        # optimized.
        set(is_debug OFF)
        if (NOT "${CMAKE_BUILD_TYPE}" STREQUAL "")
            string(TOUPPER ${CMAKE_BUILD_TYPE} upper_build_type)
            if (upper_build_type STREQUAL "DEBUG")
                set(is_debug ON)
            endif()
        endif()

        if (NOT is_debug)
            message(WARNING "[${tgt}]: Coverage can only be tracked for code with no optimizations enabled (Debug mode).")
        else()
            # Configure code coverage with gcov
            # Gcov is included by default with GNU GCC and Clang
            if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                message(WARNING "[${tgt}]: Coverage has not been implemented for the current compiler (${CMAKE_CXX_COMPILER_ID}).")
            else()
                set(is_coverage_enabled ON)

                target_compile_options(${tgt}
                    PUBLIC --coverage -O0
                )

                target_link_options(${tgt}
                    PUBLIC --coverage -O0
                )
            endif()
        endif()

        if (is_coverage_enabled)
            message(STATUS "[${tgt}]: Set flag to track code coverage.")
        else()
            message(STATUS "[${tgt}]: Did not set flag to track code coverage.")
        endif()

    endfunction()

    set_enable_coverage(dare)
    set_enable_coverage(exec)
    if(${USE_GTEST})
        set_enable_coverage(unit_test)
    endif()
endif()

# add target sources
target_sources(exec
    PUBLIC main.cpp
)

# Specify the include directories
target_include_directories(dare PUBLIC 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}>"
)

# link libraries
target_link_libraries(dare
    PUBLIC OpenMP::OpenMP_CXX
    PUBLIC MPI::MPI_CXX
    PUBLIC Boost::boost
    PUBLIC Boost::filesystem
    PUBLIC Boost::regex
    PUBLIC Boost::serialization
    PUBLIC Boost::system
    PUBLIC ${VTK_LIBRARIES}
    PUBLIC Trilinos::all_selected_libs
)

target_link_libraries(exec
    PRIVATE dare
)

if(USE_GTEST)
target_link_libraries(unit_test
    PRIVATE dare
    PUBLIC GTest::GTest
)
endif()


# # Create installation rules and configuration scripts
# install(
#     TARGETS trilinos
#     EXPORT "${PACKAGE_NAME_EXPORT}"
#     RUNTIME DESTINATION  "${CMAKE_INSTALL_BINDIR}"
#     LIBRARY DESTINATION  "${CMAKE_INSTALL_LIBDIR}/${PACKAGE_VERSION_DIR}"
#     ARCHIVE DESTINATION  "${CMAKE_INSTALL_LIBDIR}/${PACKAGE_VERSION_DIR}"
#     INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}"
# )

install(
    TARGETS dare
    EXPORT "${PACKAGE_NAME_EXPORT}"
    RUNTIME DESTINATION  "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION  "${CMAKE_INSTALL_LIBDIR}/${PACKAGE_VERSION_DIR}"
    ARCHIVE DESTINATION  "${CMAKE_INSTALL_LIBDIR}/${PACKAGE_VERSION_DIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_VERSION_DIR}"
)


install(EXPORT "${PACKAGE_NAME_EXPORT}"
    NAMESPACE ${CMAKE_PROJECT_NAME}::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PACKAGE_CMAKEDIR_NAME}"
)

export(EXPORT "${PACKAGE_NAME_EXPORT}"
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_BUILD_TREE_NAME}/${PACKAGE_NAME_EXPORT}.cmake"
    NAMESPACE ${CMAKE_PROJECT_NAME}::
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_BUILD_TREE_NAME}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${${CMAKE_PROJECT_NAME}_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/DaReConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_BUILD_TREE_NAME}/${CMAKE_PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION
        "${CMAKE_INSTALL_LIBDIR}/cmake/${PACKAGE_CMAKEDIR_NAME}"
    NO_SET_AND_CHECK_MACRO
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_BUILD_TREE_NAME}/${CMAKE_PROJECT_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_BUILD_TREE_NAME}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION 
        "${CMAKE_INSTALL_LIBDIR}/cmake/${PACKAGE_CMAKEDIR_NAME}"
)
