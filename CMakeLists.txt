
cmake_minimum_required(VERSION 4.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(ATEAMS++
	VERSION 0.1
	LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)


# # Do this instead of requiring c++23, since we need (llvm) clang for this to work.
add_compile_options(-g -O3)

# Require a bunch of packages.
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
)

include_directories(include/libraries/PHAT)

set(ATEAMS++ PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(ATEAMS++ PROPERTIES PUBLIC_HEADER include/ATEAMS++/ATEAMS++.h)

target_link_libraries(ATEAMS++ PUBLIC PkgConfig::FLINT)
target_include_directories(ATEAMS++ PUBLIC include)
target_include_directories(ATEAMS++ PUBLIC src)


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


