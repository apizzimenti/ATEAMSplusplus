
cmake_minimum_required(VERSION 4.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(ATEAMS++
	VERSION 1.0
	LANGUAGES C CXX
	DESCRIPTION "Algebraic topology-enabled algorithms for spin systems."
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-g -O3)

# Require OpenMP and FLINT.
find_package(PkgConfig REQUIRED)
find_package(OpenMP REQUIRED)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)


# Add the library.
add_library(ATEAMS++ SHARED
	src/arithmetic/kernel.cpp
	src/arithmetic/persistence.cpp
	src/complexes/Cubical.cpp
	src/models/SwendsenWang.cpp
	src/models/Bernoulli.cpp
	src/models/InvadedCluster.cpp
	src/models/Glauber.cpp
)

include_directories(dependencies/PHAT)
include_directories(dependencies/SparseRREF)

# Include dependency headers.
install(DIRECTORY "${PROJECT_SOURCE_DIR}/dependencies/phat"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
	FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/dependencies/SparseRREF"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
	FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
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
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/statistics.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/Chain.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/observables.h
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/statistics/schedules.h
)

target_link_libraries(ATEAMS++ PUBLIC PkgConfig::FLINT)
target_include_directories(ATEAMS++ PUBLIC include)
target_include_directories(ATEAMS++ PUBLIC dependencies)
target_include_directories(ATEAMS++ PUBLIC src)

include(GNUInstallDirs)

install(
	TARGETS ${PROJECT_NAME}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

configure_file(ATEAMS++.pc.in ATEAMS++.pc @ONLY)

install(
	FILES ${CMAKE_BINARY_DIR}/ATEAMS++.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)


# Tests.
enable_testing()

# Complexes.
add_executable(complexes.Cubical tests/complexes/Cubical.cpp)
target_link_libraries(complexes.Cubical PRIVATE ATEAMS++)

add_test(
	NAME complexes.Cubical
	COMMAND $<TARGET_FILE:complexes.Cubical>
)


# Models.
add_executable(models.SwendsenWang tests/models/SwendsenWang.cpp)
target_link_libraries(models.SwendsenWang PRIVATE ATEAMS++)

add_test(
	NAME models.SwendsenWang
	COMMAND $<TARGET_FILE:models.SwendsenWang>
)

add_executable(models.Bernoulli tests/models/Bernoulli.cpp)
target_link_libraries(models.Bernoulli PRIVATE ATEAMS++)

add_test(
	NAME models.Bernoulli
	COMMAND $<TARGET_FILE:models.Bernoulli>
)


add_executable(models.InvadedCluster tests/models/InvadedCluster.cpp)
target_link_libraries(models.InvadedCluster PRIVATE ATEAMS++)

add_test(
	NAME models.InvadedCluster
	COMMAND $<TARGET_FILE:models.InvadedCluster>
)

add_executable(models.Glauber tests/models/Glauber.cpp)
target_link_libraries(models.Glauber PRIVATE ATEAMS++)

add_test(
	NAME models.Glauber
	COMMAND $<TARGET_FILE:models.Glauber>
)


add_executable(arithmetic.persistence tests/arithmetic/persistence.cpp)
target_link_libraries(arithmetic.persistence PRIVATE ATEAMS++)

add_test(
	NAME arithmetic.persistence
	COMMAND $<TARGET_FILE:arithmetic.persistence>
)


add_executable(statistics.Chain tests/statistics/Chain.cpp)
target_link_libraries(statistics.Chain PRIVATE ATEAMS++)

add_test(
	NAME statistics.Chain
	COMMAND $<TARGET_FILE:statistics.Chain>
)

# add_executable(statistics.observables tests/statistics/observables.cpp)
# target_link_libraries(statistics.observables PRIVATE ATEAMS++)

# add_test(
# 	NAME statistics.observables
# 	COMMAND $<TARGET_FILE:statistics.observables>
# )


