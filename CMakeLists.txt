
################################################################################
## PROJECT METADATA ############################################################
################################################################################

cmake_minimum_required(VERSION 4.2)

project(ATEAMS++
	VERSION 1.0.0
	LANGUAGES C CXX
	DESCRIPTION "Algebraic topology-enabled algorithms for spin systems."
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################
## LOCAL BUILD CONFIG ##########################################################
################################################################################
set(TEST TRUE)
set(PROFILE FALSE)
set(TIMING FALSE)

set(DEBUG FALSE)

include(cmake/build-config.cmake)

################################################################################
## GLOBAL DEPENDENCIES #########################################################
################################################################################

# Threading building blocks (tbb) is required to use standard library parallelism
# (and parallelism more generally on linux systems on intel hardware).
find_package(PkgConfig REQUIRED)
find_package(OpenMP REQUIRED)
find_library(TBB NAMES tbb lbb REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)

# Set the choice of memory allocator. In testing on meglTower, jemalloc does
# the best. ONLY ONE OF THESE CAN BE USED AT A TIME.
set(jemalloc TRUE)
set(mimalloc FALSE)
set(tcmalloc FALSE)

# Turn OpenMP on/off.
set(OpenMP FALSE)


################################################################################
## MAIN LIBRARY BUILD ##########################################################
################################################################################

# Add library and sources.
add_library(ATEAMS++ INTERFACE)
include(cmake/library.cmake)
include(GNUInstallDirs)

target_sources(ATEAMS++
	PUBLIC
	FILE_SET HEADERS
	BASE_DIRS ${PROJECT_SOURCE_DIR}/include/ ${PROJECT_SOURCE_DIR}/src/
	FILES ${SOURCE_HEADERS} ${SOURCE_IMPLEMENTATIONS}
)

# Special rule for building locally on macOS and using gnu screen: when screen
# creates a subshell, it doesn't remember what @rpath means, so the load path
# for libATEAMS++.dylib is uninterpretable.
if(APPLE)
	set_target_properties(ATEAMS++ PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
	pkg_check_modules(TBB REQUIRED IMPORTED_TARGET tbb)
	target_link_libraries(ATEAMS++ INTERFACE PkgConfig::TBB)
endif()

# Install dependencies.
include(cmake/dependencies.cmake)

# Install the library.
set_target_properties(ATEAMS++ PROPERTIES PUBLIC_HEADER include/ATEAMS++/ATEAMS++.h)

install(
	TARGETS ${PROJECT_NAME}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the pkg-config file.
include(cmake/pkg-config.cmake)

install(
	FILES ${CMAKE_BINARY_DIR}/ATEAMS++.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)


################################################################################
## TESTING  ####################################################################
################################################################################
set(TEST_MODELS FALSE)
set(TEST_COMPLEXES FALSE)
set(TEST_ARITHMETIC FALSE)
set(TEST_TOPOLOGY TRUE)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/testing.cmake)


################################################################################
## PROFILING  ##################################################################
################################################################################
set(PROFILE_MODELS TRUE)
set(PROFILE_ARITHMETIC FALSE)
set(PROFILE_TOPOLOGY FALSE)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/profiling.cmake)

################################################################################
## TIMING ######################################################################
################################################################################
set(TIMING_MODELS FALSE)
set(TIMING_ARITHMETIC FALSE)
set(TIMING_TOPOLOGY TRUE)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/timing.cmake)


if (PROFILE OR TIMING)
	add_executable(performance.samples performance/performance.samples.cpp)
	target_link_libraries(performance.samples PUBLIC ATEAMS++)
endif()
