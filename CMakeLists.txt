cmake_minimum_required(VERSION 3.13)

project(DaRe
	VERSION 0.1.3
	DESCRIPTION "Dave's repository"
	LANGUAGES CXX
)

option(DEBUG
    "Enables debugging statements for the repository"
    OFF
)

option(USE_GTEST
    "Enables debugging statements for the repository"
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
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
set(CMAKE_FIND_PACKAGE_SORT_ORDER NATURAL)
set(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)

# set the properties of the target (might be moved to another file)
set_target_properties(exec
	       PROPERTIES CXX_STANDARD 17
		   CXX_STANDARD_REQUIRED TRUE
)

set_target_properties(dare
           PROPERTIES CXX_STANDARD 17
           VERSION 0.1
)

if(USE_GTEST)
    set_target_properties(unit_test
           PROPERTIES CXX_STANDARD 17
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

add_subdirectory(Data                                                 )
add_subdirectory(FileFormats                                          )
add_subdirectory(Utilities                                            )
add_subdirectory(MPI                                                  )
add_subdirectory(Grid                                                 )
add_subdirectory(MatrixSystem                                         )
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
endif()

# add target sources
target_sources(exec
    PUBLIC main.cpp
)

# Specify the include directories
target_include_directories(dare PUBLIC 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
)
target_include_directories(exec PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
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
    PUBLIC Trilinos::all_selected_libs
    # PUBLIC trilinos
)
target_link_libraries(exec
    PUBLIC OpenMP::OpenMP_CXX
    PUBLIC MPI::MPI_CXX
    PUBLIC Boost::boost
    PUBLIC Boost::filesystem
    PUBLIC Boost::regex
    PUBLIC Boost::serialization
    PUBLIC Boost::system
    PUBLIC Trilinos::all_selected_libs
    # PUBLIC trilinos
)
if(USE_GTEST)
target_link_libraries(dare
    PUBLIC GTest::GTest
)
target_link_libraries(exec
    PUBLIC GTest::GTest
)
target_link_libraries(unit_test PUBLIC Trilinos::all_selected_libs)
target_link_libraries(unit_test
    PUBLIC OpenMP::OpenMP_CXX
    PUBLIC MPI::MPI_CXX
    PUBLIC Boost::boost
    PUBLIC Boost::filesystem
    PUBLIC Boost::regex
    PUBLIC Boost::serialization
    PUBLIC Boost::system
    # PUBLIC Trilinos::all_selected_libs
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
