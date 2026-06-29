
cmake_minimum_required(VERSION 4.2)

project(ATEAMS++
	VERSION 1.3
	LANGUAGES C CXX
	DESCRIPTION "Algebraic topology-enabled algorithms for spin systems."
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-O3)


################################################################################
## EXTERNAL DEPENDENCIES                                                      ##
################################################################################
find_package(PkgConfig REQUIRED)
find_package(OpenMP REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)

################################################################################
## TESTING/PROFILING?                                                         ##
################################################################################
set(TEST TRUE)
set(PROFILE FALSE)

if (TEST)
	# add_subdirectory(dependencies/Catch2)
	enable_testing()
endif()

if (PROFILE)
	add_compile_options(-g)
endif()



################################################################################
## MAIN LIBRARY BUILD                                                         ##
################################################################################
# Add the library.
add_library(ATEAMS++ SHARED
	src/arithmetic/kernel.cpp
	src/arithmetic/persistence.cpp
	src/complexes/Cubical.cpp
	src/models/SwendsenWang.cpp
	src/models/Bernoulli.cpp
	src/models/InvadedCluster.cpp
	src/models/Glauber.cpp
	src/models/Invasion.cpp
)


# Public headers.
target_sources(ATEAMS++
	PUBLIC
	FILE_SET HEADERS
	BASE_DIRS ${PROJECT_SOURCE_DIR}/include/
	FILES
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/common.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/util.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/arithmetic.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/kernel.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/persistence.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/options.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/complexes.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/Complex.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/Cubical.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/models.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Bernoulli.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Glauber.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/InvadedCluster.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Model.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/SwendsenWang.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Invasion.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/statistics.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/Chain.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/observables.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/schedules.h
)

set_target_properties(ATEAMS++ PROPERTIES PUBLIC_HEADER include/ATEAMS++/ATEAMS++.h)
target_link_libraries(ATEAMS++ PUBLIC PkgConfig::FLINT)

target_include_directories(ATEAMS++ PUBLIC include)
target_include_directories(ATEAMS++ PUBLIC dependencies)
target_include_directories(ATEAMS++ PUBLIC src)

include(GNUInstallDirs)

# Special rule for building locally on macOS and using gnu screen: when screen
# creates a subshell, it doesn't remember what @rpath means, so the load path
# for libATEAMS++.dylib is uninterpretable.
if(APPLE)
	set_target_properties(ATEAMS++ PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif()


# Install headers.
install(DIRECTORY "${PROJECT_SOURCE_DIR}/dependencies/phat/"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/phat"
	FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/dependencies/SparseRREF/"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/SparseRREF"
	FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

# Install the library.
install(
	TARGETS ${PROJECT_NAME}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)


# Install the pkg-config file.
configure_file(ATEAMS++.pc.in ATEAMS++.pc @ONLY)
install(
	FILES ${CMAKE_BINARY_DIR}/ATEAMS++.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)



################################################################################
## TESTING                                                                    ##
################################################################################

list(APPEND COMPLEXES
	"Cubical"
)

list(APPEND ARITHMETIC
	"kernel"
	"persistence"
)

list(APPEND MODELS
	"Bernoulli"
	"Glauber"
	"InvadedCluster"
	"Invasion"
	"SwendsenWang"
)

if (TEST)
	foreach(MODEL IN LISTS MODELS)
		add_executable(tests.${MODEL} tests/models/${MODEL}.cpp)
		target_include_directories(tests.${MODEL} PUBLIC tests)
		target_link_libraries(tests.${MODEL} PRIVATE ATEAMS++)

		add_test(NAME tests.${MODEL}.2 COMMAND tests.${MODEL} 2)
		add_test(NAME tests.${MODEL}.3 COMMAND tests.${MODEL} 3)
	endforeach()


	foreach(COMPLEX IN LISTS COMPLEXES)
		add_executable(tests.${COMPLEX} tests/complexes/${COMPLEX}.cpp)
		target_include_directories(tests.${COMPLEX} PUBLIC tests)
		target_link_libraries(tests.${COMPLEX} PRIVATE ATEAMS++)

		add_test(NAME tests.${COMPLEX} COMMAND tests.${COMPLEX})
	endforeach()


	foreach(FUNCTION IN LISTS ARITHMETIC)
		add_executable(tests.${FUNCTION} tests/arithmetic/${FUNCTION}.cpp)
		target_include_directories(tests.${FUNCTION} PUBLIC tests)
		target_link_libraries(tests.${FUNCTION} PRIVATE ATEAMS++)

		add_test(NAME tests.${FUNCTION}.2 COMMAND tests.${FUNCTION} 2)
		add_test(NAME tests.${FUNCTION}.3 COMMAND tests.${FUNCTION} 3)
	endforeach()
endif()


if (PROFILE)
	foreach(MODEL IN LISTS MODELS)
		add_executable(profiles.${MODEL} profiles/${MODEL}.cpp)
		add_compile_options(-g)
		target_link_libraries(profiles.${MODEL} PRIVATE ATEAMS++)
	endforeach()
endif()
