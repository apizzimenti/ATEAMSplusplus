
cmake_minimum_required(VERSION 4.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(ATEAMS++
	VERSION 0.1
	LANGUAGES C CXX
)

# Do this instead of requiring c++23, since we need (llvm) clang for this to work.
add_compile_options(-std=c++23)

# Require a bunch of packages.
find_package(PkgConfig REQUIRED)
find_package(OpenMP REQUIRED)

pkg_check_modules(GIVARO REQUIRED IMPORTED_TARGET givaro)
pkg_check_modules(FFLAS_FFPACK REQUIRED IMPORTED_TARGET fflas-ffpack)
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)


# Add the library.
add_library(ATEAMS++ SHARED
	src/complexes/Cubical.cpp
)

set(ATEAMS++ PROPERTIES VERSION ${PROJECT_VERSION})
# set_target_properties(ATEAMS++ PROPERTIES PUBLIC_HEADER include/ATEAMS++/ATEAMS++.h)

target_link_libraries(ATEAMS++ PUBLIC PkgConfig::FLINT)
target_include_directories(ATEAMS++ PUBLIC include)
target_include_directories(ATEAMS++ PRIVATE src)


# Tests.
enable_testing()

# Complexes.
add_executable(complexes.Cubical tests/complexes/Cubical.cpp)
target_link_libraries(complexes.Cubical PRIVATE ATEAMS++)

add_test(
	NAME complexes.Cubical
	COMMAND $<TARGET_FILE:complexes.Cubical>
)

