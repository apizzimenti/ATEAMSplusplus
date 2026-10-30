################################################################################
## Profile configurations.                                                    ##
################################################################################


list(APPEND ARITHMETIC.PROFILING
	# "addition"
)

list(APPEND MODELS.PROFILING
	"Bernoulli"
	# "InvadedCluster"
	# "Invasion"
	# "SwendsenWang"
	# "Glauber"
)

list(APPEND TOPOLOGY.PROFILING
	"persistence"
)

if (PROFILE)
	if(PROFILE_MODELS)
		foreach(MODEL IN LISTS MODELS.PROFILING)
			add_executable(profiling.${MODEL} performance/profiling.${MODEL}.cpp)
			target_link_libraries(profiling.${MODEL} PUBLIC ATEAMS++)
			target_compile_options(profiling.${MODEL} PUBLIC ${PROFILE_COMPILE_OPTIONS})
		endforeach()
	endif()

	if(PROFILE_TOPOLOGY)
		foreach(FUNCTION IN LISTS TOPOLOGY.PROFILING)
			add_executable(profiling.${FUNCTION} performance/profiling.${FUNCTION}.cpp)
			target_link_libraries(profiling.${FUNCTION} PUBLIC ATEAMS++)
			target_compile_options(profiling.${FUNCTION} PUBLIC ${PROFILE_COMPILE_OPTIONS})
		endforeach()
	endif()

	if(PROFILE_ARITHMETIC)
		foreach(FUNCTION IN LISTS ARITHMETIC.PROFILING)
			add_executable(profiling.${FUNCTION} performance/profiling.${FUNCTION}.cpp)
			target_link_libraries(profiling.${FUNCTION} PUBLIC ATEAMS++)
			target_compile_options(profiling.${FUNCTION} PUBLIC ${PROFILE_COMPILE_OPTIONS})
		endforeach()
	endif()
endif()