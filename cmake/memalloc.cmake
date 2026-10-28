
################################################################################
## Drop-in replacement memory allocator configuration.                        ##
################################################################################

# Throw an error if any two of the allocators are used at once.
if((mimalloc AND tcmalloc) OR (mimalloc AND jemalloc) OR (tcmalloc AND jemalloc))
	message(FATAL_ERROR "Error: only one memory allocator at a time, please!")
endif()

# Send a warning if using tcmalloc on Apple.
if(APPLE AND tcmalloc)
	message(WARNING "Using tcmalloc is not recommended on macOS; aborting; switch to mimalloc or jemalloc.")
endif()

# Set library and include paths.
set(MEMALLOC_LIB_DIR )
set(MEMALLOC_LIB )

if(mimalloc)
	# If using mimalloc, find it through pkg-config.
	find_library(MIMALLOC NAMES mimalloc REQUIRED)
	pkg_check_modules(mimalloc REQUIRED IMPORTED_TARGET mimalloc)
	target_link_libraries(ATEAMS++ INTERFACE ${MIMALLOC})

	# For the main library build, enable mimalloc.
	list(APPEND SOURCE_HEADERS
		${PROJECT_SOURCE_DIR}/include/ATEAMS++/enable_mimalloc.h
	)

	# Set the memory allocator variable.
	get_filename_component(MEMALLOC_LIB_DIR "-L${MIMALLOC}" DIRECTORY)
	set(MEMALLOC_LIB "-lmimalloc")

elseif(tcmalloc AND NOT APPLE)
	# tcmalloc isn't in pkg-config (and only works on linux).
	find_library(TCMALLOC NAMES tcmalloc_minimal tcmalloc REQUIRED)
	target_link_libraries(ATEAMS++ INTERFACE ${TCMALLOC})

	get_filename_component(MEMALLOC_LIB_DIR "-L${TCMALLOC}" DIRECTORY)
	set(MEMALLOC_LIB "-ltcmalloc")

elseif(jemalloc)
	# same with jemalloc.
	find_library(JEMALLOC NAMES jemalloc jemalloc_minimal REQUIRED)
	target_link_libraries(ATEAMS++ INTERFACE ${JEMALLOC})

	get_filename_component(MEMALLOC_LIB_DIR "-L${JEMALLOC}" DIRECTORY)
	set(MEMALLOC_LIB "-ljemalloc")
endif()

set(PKG_CONFIG_MEMALLOC "${MEMALLOC_LIB_DIR} ${MEMALLOC_LIB}")
