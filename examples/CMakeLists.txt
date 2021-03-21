SET(${PROJECT_NAME}_EXAMPLES
	# DoubleIntegrator
	# InvertedPendulum_Robust
	# InvertedPendulum_Realizable
	# InvertedPendulum_RealizableSampled
	DoubleIntegrator_Robust
	DoubleIntegrator_RealizableSampled
	)

FOREACH(EXAMPLE ${${PROJECT_NAME}_EXAMPLES})
	set(EXAMPLE_NAME "${EXAMPLE}-example")
	add_executable(${EXAMPLE_NAME} "${EXAMPLE}.cpp")
	target_link_libraries(${EXAMPLE_NAME} PRIVATE ${PROJECT_NAME})
ENDFOREACH(EXAMPLE ${${PROJECT_NAME}_EXAMPLES})