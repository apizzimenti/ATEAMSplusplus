
################################################################################
## Linking/including dependencies for local build.                            ##
################################################################################

# List core external dependencies.
list(APPEND COREDEPENDENCIES
	"PkgConfig::FLINT"
	"${TBB}"
)

if(OpenMP)
	list(APPEND COREDEPENDENCIES
		"OpenMP::OpenMP_CXX"
	)
else()
endif()

# Link external dependencies.
foreach(DEP IN LISTS COREDEPENDENCIES)
	target_link_libraries(ATEAMS++ INTERFACE ${DEP})
endforeach()

# Link memory allocator dependencies.
include(cmake/memalloc.cmake)

# Include header files for local dependencies, then link them.
target_include_directories(ATEAMS++ INTERFACE dependencies)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/dependencies/phat/"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/phat"
	FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/dependencies/SparseRREF/"
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/SparseRREF"
	FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)