
cmake_minimum_required(VERSION 4.2)

project(ATEAMS++
	VERSION 1.5
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
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/persistence.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/arithmetic/options.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/complexes.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/Complex.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/complexes/Cubical.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/models.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/ModelState.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/ModelParameters.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Model.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Bernoulli.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Glauber.h
	# ${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/InvadedCluster.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/SwendsenWang.h
	# ${PROJECT_SOURCE_DIR}/include/ATEAMS++/models/Invasion.h

	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/statistics.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/Chain.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/observables.h
	${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/schedules.h
)


list(APPEND SOURCE_IMPLEMENTATIONS
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/arithmetic/kernel.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/arithmetic/persistence.tpp

	${PROJECT_SOURCE_DIR}/src/ATEAMS++/complexes/Cubical.tpp

	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/SwendsenWang.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/Bernoulli.tpp
	${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/Glauber.tpp
	# ${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/InvadedCluster.tpp
	# ${PROJECT_SOURCE_DIR}/src/ATEAMS++/models/Invasion.tpp
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
## TESTING                                                                    ##
################################################################################

list(APPEND COMPLEXES
	# "Cubical"
)

list(APPEND ARITHMETIC
	# "kernel"
	# "persistence"
)

list(APPEND MODELS
	# "Bernoulli"
	"Glauber"
	# "InvadedCluster"
	# "Invasion"
	# "SwendsenWang"
)

list(APPEND TYPES
	"Q"
	"Zp"
)

if (TEST)
	foreach(MODEL IN LISTS MODELS)
		foreach(TYPE IN LISTS TYPES)
			add_executable(tests.${MODEL}.${TYPE} tests/models/${MODEL}.${TYPE}.cpp)
			target_include_directories(tests.${MODEL}.${TYPE} PUBLIC tests)
			target_link_libraries(tests.${MODEL}.${TYPE} PUBLIC ATEAMS++)

			if("${TYPE}" STREQUAL "Q")
				add_test(NAME tests.${MODEL}.${TYPE} COMMAND tests.${MODEL}.${TYPE} 0)
			else()
				add_test(NAME tests.${MODEL}.${TYPE}.2 COMMAND tests.${MODEL}.${TYPE} 2)
				add_test(NAME tests.${MODEL}.${TYPE}.3 COMMAND tests.${MODEL}.${TYPE} 3)
			endif()
		endforeach()
	endforeach()


	foreach(COMPLEX IN LISTS COMPLEXES)
		foreach(TYPE IN LISTS TYPES)
			add_executable(tests.${COMPLEX}.${TYPE} tests/complexes/${COMPLEX}.${TYPE}.cpp)
			target_include_directories(tests.${COMPLEX}.${TYPE} PUBLIC tests)
			target_link_libraries(tests.${COMPLEX}.${TYPE} PUBLIC ATEAMS++)
			add_test(NAME tests.${COMPLEX}.${TYPE} COMMAND tests.${COMPLEX}.${TYPE})
		endforeach()
	endforeach()


	foreach(FUNCTION IN LISTS ARITHMETIC)
		foreach(TYPE IN LISTS TYPES)
			add_executable(tests.${FUNCTION}.${TYPE} tests/arithmetic/${FUNCTION}.${TYPE}.cpp)
			target_include_directories(tests.${FUNCTION}.${TYPE} PUBLIC tests)
			target_link_libraries(tests.${FUNCTION}.${TYPE} PUBLIC ATEAMS++)

			if("${TYPE}" STREQUAL "Q")
				add_test(NAME tests.${FUNCTION}.${TYPE} COMMAND tests.${FUNCTION}.${TYPE} 0)
			else()
				add_test(NAME tests.${FUNCTION}.${TYPE}.2 COMMAND tests.${FUNCTION}.${TYPE} 2)
				add_test(NAME tests.${FUNCTION}.${TYPE}.3 COMMAND tests.${FUNCTION}.${TYPE} 3)
			endif()
		endforeach()
	endforeach()
endif()


if (PROFILE)
	foreach(MODEL IN LISTS MODELS)
		add_executable(profiles.${MODEL} profiles/${MODEL}.cpp)
		add_compile_options(-g)
		target_link_libraries(profiles.${MODEL} PRIVATE ATEAMS++)
	endforeach()
endif()
