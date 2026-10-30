
################################################################################
## Local build configuration based on user-supplied parameters.               ##
################################################################################

# If we aren't profiling or debugging, we can omit frame pointers.
if(NOT DEBUG)
	if(TIMING)
		set(TIMING_COMPILE_OPTIONS
			-Ofast
			-march=native
			-flto
			-funroll-loops
			-fprefetch-loop-arrays
			-fno-rtti
			-fomit-frame-pointer
		)

		set(TIMING_LINK_OPTIONS
			-Ofast
			-march=native
			-flto
			-funroll-loops
			-fprefetch-loop-arrays
			-fno-rtti
			-fomit-frame-pointer
		)
	endif()

	if(TEST)
		set(TESTING_COMPILE_OPTIONS
			-Ofast
			-march=native
			-flto
			-funroll-loops
			-fprefetch-loop-arrays
			-fno-rtti
			-fomit-frame-pointer
		)
	endif()


	if(PROFILE)
		set(PROFILE_COMPILE_OPTIONS
			-ggdb3
			-Ofast
			-march=native
			-flto
			-funroll-loops
			-fprefetch-loop-arrays
			-fno-rtti
			-fno-omit-frame-pointer
		)

		set(PROFILE_LINK_OPTIONS
			-ggdb3
			-Ofast
			-march=native
			-flto
			-funroll-loops
			-fprefetch-loop-arrays
			-fno-rtti
			-fno-omit-frame-pointer
		)
	endif()
endif()


if(DEBUG)
	add_compile_options(
		-O1
		-g
	)
	add_link_options(
		-O1
		-g
	)
endif()
