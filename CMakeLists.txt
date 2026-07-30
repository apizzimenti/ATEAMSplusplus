
cmake_minimum_required(VERSION 4.2)

project(ATEAMS++
	VERSION 1.8
	LANGUAGES C CXX
	DESCRIPTION "Algebraic topology-enabled algorithms for spin systems."
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-march=native -flto)


################################################################################
## EXTERNAL DEPENDENCIES                                                      ##
################################################################################
find_package(PkgConfig REQUIRED)
find_package(OpenMP REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)

# If we want to use mimalloc (which we should, by default) then set this to
# TRUE. Otherwise, set it to FALSE.
set(mimalloc FALSE)

if(mimalloc)
	find_package(mimalloc REQUIRED)
	pkg_check_modules(mimalloc REQUIRED IMPORTED_TARGET mimalloc)
endif()

################################################################################
## TESTING/PROFILING?                                                         ##
################################################################################
set(TEST TRUE)
set(PROFILE FALSE)
set(TIMING FALSE)
set(DEBUG FALSE)

set(PROFILE_COMPILE_OPTIONS
	-g
	-fno-omit-frame-pointer
	-fno-inline-functions
)

set(PROFILE_LINK_OPTIONS
	-g
)

if (TEST)
	enable_testing()
endif()

if (DEBUG)
	add_compile_options(-g -O0)
else()
	add_compile_options(-Ofast)
endif()




################################################################################
## MAIN LIBRARY BUILD                                                         ##
################################################################################
add_library(ATEAMS++ INTERFACE)

if(mimalloc)
	list(APPEND SOURCE_HEADERS
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/enable_mimalloc.h
	)
endif()


list(APPEND SOURCE_HEADERS
	# Headers
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/common.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/util.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/arithmetic.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/kernel.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/options.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/topology.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence/JIT.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence/parallel.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence/PHAT.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence/split.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence/stagger.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence/standard.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/persistence/twist.h
	# ${PROJECT_SOURCE_DIR}/include/ATEAMS++/topology/helpers.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/complexes.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/Complex.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/Cubical.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/models.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/ModelState.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/ModelParameters.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Model.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Bernoulli.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Glauber.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/InvadedCluster.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/SwendsenWang.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Invasion.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/statistics.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/Chain.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/observables.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/schedules.h
)


list(APPEND SOURCE_IMPLEMENTATIONS
	# Templated implementations
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/arithmetic/kernel.tpp
	
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/topology/persistence/JIT.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/topology/persistence/parallel.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/topology/persistence/PHAT.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/topology/persistence/split.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/topology/persistence/stagger.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/topology/persistence/standard.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/topology/persistence/twist.tpp

	${PROJECT_SOURCE_DIR}/src/ATEAMS++/complexes/Cubical.tpp

	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/SwendsenWang.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/Bernoulli.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/Glauber.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/InvadedCluster.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/Invasion.tpp
)


# Public headers.
target_sources(ATEAMS++
	PUBLIC
	FILE_SET HEADERS
	BASE_DIRS ${PROJECT_SOURCE_DIR}/include/ ${PROJECT_SOURCE_DIR}/src/
	FILES ${SOURCE_HEADERS} ${SOURCE_IMPLEMENTATIONS}
)

set_target_properties(ATEAMS++ PROPERTIES PUBLIC_HEADER include/ATEAMS++/ATEAMS++.h)
target_link_libraries(ATEAMS++ INTERFACE PkgConfig::FLINT)
target_link_libraries(ATEAMS++ INTERFACE OpenMP::OpenMP_CXX)

if(mimalloc)
	target_link_libraries(ATEAMS++ INTERFACE PkgConfig::mimalloc)
endif()

# target_include_directories(ATEAMS++ INTERFACE include)
target_include_directories(ATEAMS++ INTERFACE dependencies)
# target_include_directories(ATEAMS++ INTERFACE src)

include(GNUInstallDirs)

# Special rule for building locally on macOS and using gnu screen: when screen
# creates a subshell, it doesn't remember what @rpath means, so the load path
# for libATEAMS++.dylib is uninterpretable.
if(APPLE)
	set_target_properties(ATEAMS++ PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif()

if(LINUX)
	find_package(TBB REQUIRED)
	link_libraries(tbb)
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
if (mimalloc)
	configure_file(ATEAMS++.pc.mimalloc.in ATEAMS++.pc @ONLY)
else()
	configure_file(ATEAMS++.pc.in ATEAMS++.pc @ONLY)
endif()

install(
	FILES ${CMAKE_BINARY_DIR}/ATEAMS++.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)



################################################################################
## TESTING  ####################################################################
################################################################################

list(APPEND COMPLEXES.TESTING
	# "Cubical"
)

list(APPEND ARITHMETIC.TESTING
	# "addition"
)

list(APPEND TOPOLOGY.TESTING
	# "split"
	# "stagger"
	"JIT"
	# "standard"
	# "twist"
	# "parallel"
	# "persistence"
)

list(APPEND MODELS.TESTING
	# "Bernoulli"
	# "InvadedCluster"
	# "Invasion"
	# "SwendsenWang"
	# "Glauber"
)

list(APPEND TYPES
	# "Q"
	"Zp"
)

if (TEST)
	foreach(TYPE IN LISTS TYPES)
		# Making sure models are calling all the required subroutines appropriately
		foreach(MODEL IN LISTS MODELS.TESTING)
			add_executable(tests.${MODEL}.${TYPE} tests/models/${MODEL}.${TYPE}.cpp)
			target_include_directories(tests.${MODEL}.${TYPE} PUBLIC tests)
			target_link_libraries(tests.${MODEL}.${TYPE} PUBLIC ATEAMS++)

			if("${TYPE}" STREQUAL "Q")
				add_test(NAME ${MODEL}.${TYPE} COMMAND tests.${MODEL}.${TYPE} 0)
			else()
				add_test(NAME ${MODEL}.${TYPE}.2 COMMAND tests.${MODEL}.${TYPE} 2)
				add_test(NAME ${MODEL}.${TYPE}.3 COMMAND tests.${MODEL}.${TYPE} 3)
			endif()
		endforeach()

		# Making sure complexes work the way we want
		foreach(COMPLEX IN LISTS COMPLEXES.TESTING)
			add_executable(tests.${COMPLEX}.${TYPE} tests/complexes/${COMPLEX}.${TYPE}.cpp)
			target_include_directories(tests.${COMPLEX}.${TYPE} PUBLIC tests)
			target_link_libraries(tests.${COMPLEX}.${TYPE} PUBLIC ATEAMS++)

			add_test(NAME ${COMPLEX}.${TYPE} COMMAND tests.${COMPLEX}.${TYPE})
		endforeach()

		# Arithmetic subroutines
		foreach(FUNCTION IN LISTS ARITHMETIC.TESTING)
			add_executable(tests.${FUNCTION}.${TYPE} tests/arithmetic/${FUNCTION}.${TYPE}.cpp)
			target_include_directories(tests.${FUNCTION}.${TYPE} PUBLIC tests)
			target_link_libraries(tests.${FUNCTION}.${TYPE} PUBLIC ATEAMS++)

			if("${TYPE}" STREQUAL "Q")
				add_test(NAME ${FUNCTION}.${TYPE} COMMAND tests.${FUNCTION}.${TYPE} 0)
			else()
				# Test at 1,000,000 entries over 100 trials.
				add_test(NAME ${FUNCTION}.${TYPE}.2.e6.100 COMMAND tests.${FUNCTION}.${TYPE} 2 1000000 100)
				add_test(NAME ${FUNCTION}.${TYPE}.3.e6.100 COMMAND tests.${FUNCTION}.${TYPE} 3 1000000 100)

				# Test at 10,000,000 entries over 10 trials.
				add_test(NAME ${FUNCTION}.${TYPE}.2.e7.10 COMMAND tests.${FUNCTION}.${TYPE} 2 10000000 10)
				add_test(NAME ${FUNCTION}.${TYPE}.3.e7.10 COMMAND tests.${FUNCTION}.${TYPE} 3 10000000 10)
			endif()
		endforeach()

		# Topological support subroutines
		foreach(FUNCTION IN LISTS TOPOLOGY.TESTING)
			add_executable(tests.persistence.${FUNCTION}.${TYPE} tests/topology/persistence.${FUNCTION}.${TYPE}.cpp)
			target_include_directories(tests.persistence.${FUNCTION}.${TYPE} PUBLIC tests)
			target_link_libraries(tests.persistence.${FUNCTION}.${TYPE} PUBLIC ATEAMS++)

			if("${TYPE}" STREQUAL "Q")
				add_test(NAME topology::persistence::${FUNCTION}.${TYPE} COMMAND tests.persistence.${FUNCTION}.${TYPE} 0)
			else()
				add_test(NAME topology::persistence::${FUNCTION}.${TYPE}.2 COMMAND tests.persistence.${FUNCTION}.${TYPE} 2)
				add_test(NAME topology::persistence::${FUNCTION}.${TYPE}.3 COMMAND tests.persistence.${FUNCTION}.${TYPE} 3)
			endif()
		endforeach()
	endforeach()
endif()




################################################################################
## PROFILING  ##################################################################
################################################################################

list(APPEND ARITHMETIC.PROFILING
	# "addition"
)

list(APPEND MODELS.PROFILING
	# "Bernoulli"
	# "InvadedCluster"
	# "Invasion"
	# "SwendsenWang"
	# "Glauber"
)

list(APPEND TOPOLOGY.PROFILING
	"persistence"
)

if (PROFILE)
	foreach(MODEL IN LISTS MODELS.PROFILING)
		add_executable(profiling.${MODEL} performance/profiling.${MODEL}.cpp)
		target_link_libraries(profiling.${MODEL} PUBLIC ATEAMS++)
		target_compile_options(profiling.${MODEL} PUBLIC ${PROFILE_COMPILE_OPTIONS})
	endforeach()

	foreach(FUNCTION IN LISTS TOPOLOGY.PROFILING)
		add_executable(profiling.${FUNCTION} performance/profiling.${FUNCTION}.cpp)
		target_link_libraries(profiling.${FUNCTION} PUBLIC ATEAMS++)
		target_link_options(profiling.${FUNCTION} PUBLIC ${PROFILE_COMPILE_OPTIONS})
		target_compile_options(profiling.${FUNCTION} PUBLIC ${PROFILE_COMPILE_OPTIONS})
	endforeach()

	foreach(FUNCTION IN LISTS ARITHMETIC.PROFILING)
		add_executable(profiling.${FUNCTION} performance/profiling.${FUNCTION}.cpp)
		target_link_libraries(profiling.${FUNCTION} PUBLIC ATEAMS++)
		target_compile_options(profiling.${FUNCTION} PUBLIC ${PROFILE_COMPILE_OPTIONS})
		target_link_options(profiling.${FUNCTION} PUBLIC ${PROFILE_COMPILE_OPTIONS})
	endforeach()
endif()


if (PROFILE OR TIMING)
	add_executable(performance.samples performance/performance.samples.cpp)
	target_link_libraries(performance.samples PUBLIC ATEAMS++)
endif()


################################################################################
## TIMING ######################################################################
################################################################################


list(APPEND ARITHMETIC.TIMING
	# "addition"
)

list(APPEND MODELS.TIMING
	# "Bernoulli"
	# "InvadedCluster"
	# "Invasion"
	# "SwendsenWang"
	# "Glauber"
)

list(APPEND TOPOLOGY.TIMING
	"persistence"
)

if (TIMING)
	foreach(MODEL IN LISTS MODELS.TIMING)
		add_executable(timing.${MODEL} performance/timing.${MODEL}.cpp)
		target_link_libraries(timing.${MODEL} PUBLIC ATEAMS++)
	endforeach()

	foreach(FUNCTION IN LISTS TOPOLOGY.TIMING)
		add_executable(timing.${FUNCTION} performance/timing.${FUNCTION}.cpp)
		target_link_libraries(timing.${FUNCTION} PUBLIC ATEAMS++)
	endforeach()

	foreach(FUNCTION IN LISTS ARITHMETIC.TIMING)
		add_executable(timing.${FUNCTION} performance/timing.${FUNCTION}.cpp)
		target_link_libraries(timing.${FUNCTION} PUBLIC ATEAMS++)
	endforeach()
endif()