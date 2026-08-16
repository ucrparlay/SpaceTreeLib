# Unit tests deliberately depend on PSI alone: tests/test_framework.h costs
# ~70 s and 1.7 GB per translation unit and drags in the baseline libraries,
# which have no business in a correctness gate.
set(psi_unit_tests
	test_queries
	test_updates
	test_degenerate
)

foreach(psi_test ${psi_unit_tests})
	add_executable(${psi_test} ${CMAKE_CURRENT_SOURCE_DIR}/${psi_test}.cpp)
	target_link_libraries(${psi_test} PRIVATE PSI)
	target_include_directories(${psi_test} PRIVATE ${PSI_SOURCE_DIR})
	# -Wall -Wextra so the tests themselves stay clean, minus the one
	# warning every p_tree instantiation produces from the deferred cpam
	# fork (augmented_node.h leaves half a std::pair indeterminate). The
	# examples already surface it; 21 more copies say nothing new.
	target_compile_options(${psi_test} PRIVATE ${PSI_WARNINGS}
	                       -Wno-maybe-uninitialized)
	add_test(NAME ${psi_test} COMMAND ${psi_test})
endforeach()
