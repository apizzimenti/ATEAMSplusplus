################################################################################
## Test configurations.                                                       ##
################################################################################

# If we're testing, enable the testing suite.
if (TEST)
	enable_testing()
endif()

list(APPEND COMPLEXES.TESTING
	"Cubical"
)

list(APPEND ARITHMETIC.TESTING
	"matmul"
	"kernel"
)

list(APPEND TOPOLOGY.TESTING
	"JIT"
	"parallel"
	"PHAT"
	"split"
	"stagger"
	"standard"
	"twist"
)

list(APPEND MODELS.TESTING
	"Bernoulli"
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
		if(TEST_MODELS)
			foreach(MODEL IN LISTS MODELS.TESTING)
				add_executable(tests.${MODEL}.${TYPE} tests/models/${MODEL}.${TYPE}.cpp)
				target_compile_options(tests.${MODEL}.${TYPE} PUBLIC ${TESTING_COMPILE_OPTIONS})

				target_include_directories(tests.${MODEL}.${TYPE} PUBLIC tests)
				target_link_libraries(tests.${MODEL}.${TYPE} PUBLIC ATEAMS++)

				if("${TYPE}" STREQUAL "Q")
					add_test(NAME models::${MODEL}.${TYPE} COMMAND tests.${MODEL}.${TYPE} 0)
				else()
					add_test(NAME models::${MODEL}.${TYPE}.2 COMMAND tests.${MODEL}.${TYPE} 2)
					add_test(NAME models::${MODEL}.${TYPE}.3 COMMAND tests.${MODEL}.${TYPE} 3)
				endif()
			endforeach()
		endif()

		# Making sure complexes work the way we want
		if(TEST_COMPLEXES)
			foreach(COMPLEX IN LISTS COMPLEXES.TESTING)
				add_executable(tests.${COMPLEX}.${TYPE} tests/complexes/${COMPLEX}.${TYPE}.cpp)
				target_include_directories(tests.${COMPLEX}.${TYPE} PUBLIC tests)
				target_link_libraries(tests.${COMPLEX}.${TYPE} PUBLIC ATEAMS++)
				target_compile_options(tests.${COMPLEX}.${TYPE} PUBLIC ${TESTING_COMPILE_OPTIONS})

				add_test(NAME ${COMPLEX}.${TYPE} COMMAND tests.${COMPLEX}.${TYPE})
			endforeach()
		endif()

		# Arithmetic subroutines
		if(TEST_ARITHMETIC)
			foreach(FUNCTION IN LISTS ARITHMETIC.TESTING)
				add_executable(tests.${FUNCTION}.${TYPE} tests/arithmetic/${FUNCTION}.${TYPE}.cpp)
				target_include_directories(tests.${FUNCTION}.${TYPE} PUBLIC tests)
				target_link_libraries(tests.${FUNCTION}.${TYPE} PUBLIC ATEAMS++)
				target_compile_options(tests.${FUNCTION}.${TYPE} PUBLIC ${TESTING_COMPILE_OPTIONS})

				set(COPIES 1000)

				if("${TYPE}" STREQUAL "Q")
					add_test(NAME ${FUNCTION}.${TYPE} COMMAND tests.${FUNCTION}.${TYPE} 0)
				else()
					# Test for Z/2Z
					add_test(NAME ${FUNCTION}.${TYPE}.2.100x100 COMMAND tests.${FUNCTION}.${TYPE} 2 100 100 ${COPIES})
					add_test(NAME ${FUNCTION}.${TYPE}.2.100x1000 COMMAND tests.${FUNCTION}.${TYPE} 2 100 1000 ${COPIES})
					add_test(NAME ${FUNCTION}.${TYPE}.2.1000x100 COMMAND tests.${FUNCTION}.${TYPE} 2 1000 100 ${COPIES})

					# Test for Z/3Z
					add_test(NAME ${FUNCTION}.${TYPE}.3.100x100 COMMAND tests.${FUNCTION}.${TYPE} 3 100 100 ${COPIES})
					add_test(NAME ${FUNCTION}.${TYPE}.3.100x1000 COMMAND tests.${FUNCTION}.${TYPE} 3 100 1000 ${COPIES})
					add_test(NAME ${FUNCTION}.${TYPE}.3.1000x100 COMMAND tests.${FUNCTION}.${TYPE} 3 1000 100 ${COPIES})
				endif()
			endforeach()
		endif()

		# Topological support subroutines
		if(TEST_TOPOLOGY)
			foreach(FUNCTION IN LISTS TOPOLOGY.TESTING)
				add_executable(tests.persistence.${FUNCTION}.${TYPE} tests/topology/persistence.${FUNCTION}.${TYPE}.cpp)
				target_include_directories(tests.persistence.${FUNCTION}.${TYPE} PUBLIC tests)
				target_link_libraries(tests.persistence.${FUNCTION}.${TYPE} PUBLIC ATEAMS++)
				target_compile_options(tests.persistence.${FUNCTION}.${TYPE} PUBLIC ${TESTING_COMPILE_OPTIONS})

				if("${TYPE}" STREQUAL "Q")
					add_test(NAME topology::persistence::${FUNCTION}.${TYPE} COMMAND tests.persistence.${FUNCTION}.${TYPE} 0)
				else()
					add_test(NAME topology::persistence::${FUNCTION}.${TYPE}.2 COMMAND tests.persistence.${FUNCTION}.${TYPE} 2)
					add_test(NAME topology::persistence::${FUNCTION}.${TYPE}.3 COMMAND tests.persistence.${FUNCTION}.${TYPE} 3)
				endif()
			endforeach()
		endif()
	endforeach()
endif()


