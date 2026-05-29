
cmake_minimum_required(VERSION 4.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(ATEAMS++
	VERSION 0.1
	LANGUAGES C CXX
)

# Do this instead of requiring c++20, since we need (llvm) clang for this to work.
add_compile_options(-std=c++20 -fexperimental-library)

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


# Add an executable just for the fun.
add_executable(main tests/main.cpp)
target_link_libraries(main PRIVATE ATEAMS++)
