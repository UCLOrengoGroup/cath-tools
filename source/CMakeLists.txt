#####################################################################
### [1/4] Preamble
#####################################################################

include( ./auto_generated_file_list.cmake )

#####################################################################

option( BUILD_EXTRA_CATH_TESTS "BUILD_EXTRA_CATH_TESTS" $ENV{BUILD_EXTRA_CATH_TESTS} )
option( BUILD_EXTRA_CATH_TOOLS "BUILD_EXTRA_CATH_TOOLS" $ENV{BUILD_EXTRA_CATH_TOOLS} )

#####################################################################

# Add INTERFACE target for rapidjson
add_library                ( cath_tools_rapidjson INTERFACE                           )
target_link_libraries      ( cath_tools_rapidjson INTERFACE RapidJSON::RapidJSON      )
target_compile_definitions ( cath_tools_rapidjson INTERFACE RAPIDJSON_HAS_STDSTRING=1 )

# Add INTERFACE target for gnuplot-iostream
add_library               ( cath_tools_gnuplot_iostream        INTERFACE                                                                         )
target_include_directories( cath_tools_gnuplot_iostream SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party_code/gnuplot> )


#####################################################################
### [2/4] Modules
#####################################################################

# 1. Specify the modules to build

set( MODULE_STEMS
	biocore
	cath_assign_domains
	cath_cluster
	cath_refine_align
	cath_score_align
	cath_superpose
	chopping
	clustagglom
	cluster
	common
	display_colour
	external_info
	options
	resolve_hits
	seq
	uni

	test
)

# 2. Set up ct_${MODULE} and testsrcs_${MODULE} for each MODULE

foreach( MODULE_STEM ${MODULE_STEMS} )
	set( MODULE "ct_${MODULE_STEM}" )
	string( TOUPPER "${MODULE}" MODULE_UPPER )

	# Create a library for the module with all the .cpp files
	# and add the directory as a public include directory
	add_library( "${MODULE}" ${NORMSOURCES_${MODULE_UPPER}} )
	
	# Add the module's include directory to its public include directories
	target_include_directories( "${MODULE}" PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}>" )

	if ( DEFINED TESTSOURCES_${MODULE_UPPER} )
		# Create the OBJECT library of test sources for this module and link it to the appropriate dependencies (determined above)
		add_library( "testsrcs_${MODULE_STEM}" OBJECT ${TESTSOURCES_${MODULE_UPPER}} )
		target_link_libraries( "testsrcs_${MODULE_STEM}" PUBLIC ct_test "${MODULE}" )

		# Record that list of test sources that are being handled (which can be checked later for omissions)
		list( APPEND ALL_MODULE_TESTSOURCES ${TESTSOURCES_${MODULE_UPPER}} )
	endif()
endforeach()

# For external_info, add the corresponding subdir of the *build* (not source) directory where configure_file() places the populated header template
target_include_directories( ct_external_info PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/source/ct_external_info> )

# 3. Specify the modules' dependencies

target_link_libraries( testsrcs_test          PUBLIC  ct_chopping                                           )
target_link_libraries( testsrcs_test          PUBLIC  ct_uni                                                ) # ct_uni for file/pdb/pdb.hpp
target_link_libraries( testsrcs_options       PUBLIC  ct_uni                                                ) # ct_uni for alignment/options_block/alignment_input_options_block.hpp

target_link_libraries( ct_common              PRIVATE cath_tools_gsl                                        )

target_link_libraries( ct_biocore             PUBLIC  ct_common                                             )
target_link_libraries( ct_cath_assign_domains PUBLIC  ct_biocore ct_options ct_uni                          ) # ct_uni for acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp
target_link_libraries( ct_cath_cluster        PUBLIC  ct_clustagglom ct_options                             )
target_link_libraries( ct_cath_refine_align   PUBLIC  ct_cath_superpose ct_common ct_uni                    ) # ct_uni for acquirer/alignment_acquirer/alignment_acquirer.hpp and alignment/options_block/alignment_input_options_block.hpp
target_link_libraries( ct_cath_score_align    PUBLIC  ct_biocore ct_cath_superpose ct_options ct_uni        ) # ct_uni for alignment/options_block/alignment_input_options_block.hpp
target_link_libraries( ct_cath_superpose      PUBLIC  ct_common ct_uni                                      ) # ct_uni for acquirer/alignment_acquirer/align_refining.hpp acquirer/alignment_acquirer/alignment_acquirer.hpp
target_link_libraries( ct_chopping            PUBLIC  ct_biocore Boost::program_options                     )
target_link_libraries( ct_clustagglom         PUBLIC  ct_common Boost::program_options                      )
target_link_libraries( ct_cluster             PUBLIC  ct_common ct_options ct_seq Boost::program_options    )
target_link_libraries( ct_common              PUBLIC  Boost::log cath_tools_rapidjson                       )
target_link_libraries( ct_display_colour      PUBLIC  ct_common                                             )
target_link_libraries( ct_options             PUBLIC  ct_chopping ct_external_info                          )
target_link_libraries( ct_resolve_hits        PUBLIC  ct_display_colour ct_options ct_seq                   )
target_link_libraries( ct_seq                 PUBLIC  ct_common                                             )
target_link_libraries( ct_test                PUBLIC  ct_common ct_external_info Boost::unit_test_framework )
target_link_libraries( ct_uni                 PUBLIC  ct_display_colour ct_external_info ct_options cath_tools_gnuplot_iostream Boost::iostreams Boost::serialization )

# Headers wanted from ct_uni by other targets:
#   1 acquirer/alignment_acquirer/align_refining.hpp
#   2 acquirer/alignment_acquirer/alignment_acquirer.hpp
#   1 acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp
#   3 alignment/options_block/alignment_input_options_block.hpp
#   3 file/pdb/pdb.hpp
#   1 file/prc_scores_file/prc_scores_file.hpp
#   1 scan/scan_tools/all_vs_all.hpp
#   1 ssap/options/cath_ssap_options.hpp


#####################################################################
### [3/4] (Non-test) executables
#####################################################################

# 1. Specify the (non-test) executables to build

list( APPEND
	NON_TEST_EXES
		cath-assign-domains
		cath-cluster
		cath-map-clusters
		cath-refine-align
		cath-resolve-hits
		cath-score-align
		cath-ssap
		cath-superpose
)
if ( BUILD_EXTRA_CATH_TOOLS )
	list( APPEND
		NON_TEST_EXES
			cath-extract-pdb
			check-pdb
			snap-judgement
	)
endif()

# 2. Set up (non-test) executables

foreach( NON_TEST_EXE ${NON_TEST_EXES} )
	string               ( TOUPPER "NORMSOURCES_EXECUTABLES_${NON_TEST_EXE}" NON_TEST_EXE_FILES_VAR )
	string               ( REPLACE "-" "_" NON_TEST_EXE_FILES_VAR "${NON_TEST_EXE_FILES_VAR}" )
	add_executable       ( "${NON_TEST_EXE}" ${${NON_TEST_EXE_FILES_VAR}} )
	set_target_properties( "${NON_TEST_EXE}" PROPERTIES RUNTIME_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/bin> )
	install              ( TARGETS "${NON_TEST_EXE}" DESTINATION bin )
endforeach()

# 3. Specify the (non-test) executables' dependencies

target_link_libraries( cath-assign-domains PRIVATE ct_cath_assign_domains ct_uni ) # ct_uni for file/prc_scores_file/prc_scores_file.hpp
target_link_libraries( cath-cluster        PRIVATE ct_cath_cluster               )
target_link_libraries( cath-map-clusters   PRIVATE ct_cluster                    )
target_link_libraries( cath-refine-align   PRIVATE ct_cath_refine_align          )
target_link_libraries( cath-resolve-hits   PRIVATE ct_resolve_hits               )
target_link_libraries( cath-score-align    PRIVATE ct_cath_score_align           )
target_link_libraries( cath-ssap           PRIVATE ct_chopping            ct_uni ) # ct_uni for ssap/options/cath_ssap_options.hpp
target_link_libraries( cath-superpose      PRIVATE ct_cath_superpose             )

IF ( BUILD_EXTRA_CATH_TOOLS )
		target_link_libraries( cath-extract-pdb PRIVATE ct_uni ) # ct_uni for file/pdb/pdb.hpp
		target_link_libraries( check-pdb        PRIVATE ct_uni ) # ct_uni for file/pdb/pdb.hpp
		target_link_libraries( snap-judgement   PRIVATE ct_uni ) # ct_uni for scan/scan_tools/all_vs_all.hpp
ENDIF()


#####################################################################
### [4/4] Tests
#####################################################################

# Check that all the TESTSOURCES are accounted for
set( TESTSOURCES_REMAIN ${TESTSOURCES} )
list( REMOVE_ITEM TESTSOURCES_REMAIN ${TESTSOURCES_EXECUTABLES_BUILD_TEST} ${ALL_MODULE_TESTSOURCES} )
if ( TESTSOURCES_REMAIN )
	message( FATAL_ERROR "Some TESTSOURCES are not being handled : \"${TESTSOURCES_REMAIN}\"" )
endif()

#####################################################################

# Add build-test, containing all the tests
add_executable( build-test ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
foreach( MODULE_STEM ${MODULE_STEMS} )
	if ( TARGET "testsrcs_${MODULE_STEM}" )
		target_link_libraries( build-test PRIVATE "testsrcs_${MODULE_STEM}" )
	endif()
endforeach()
add_test( NAME build-test COMMAND build-test )
set_target_properties( build-test PROPERTIES RUNTIME_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/bin> )

#####################################################################

# If requested, build a mod-test-XXX for each module
if ( BUILD_EXTRA_CATH_TESTS )

	foreach( MODULE_STEM ${MODULE_STEMS} )
		if ( TARGET "testsrcs_${MODULE_STEM}" )
			# Make the test name (eg resolve_hits -> mod-test-resolve-hits)
			string               ( REPLACE "_" "-" TESTED_MODULE_EXE "mod-test-${MODULE_STEM}" )

			# Make the test, link it to the testsrcs OBJECT library and register it as a test
			add_executable       ( "${TESTED_MODULE_EXE}" ${TESTSOURCES_EXECUTABLES_BUILD_TEST}  )
			target_link_libraries( "${TESTED_MODULE_EXE}" "testsrcs_${MODULE_STEM}"                   )
			add_test             ( NAME "${TESTED_MODULE_EXE}" COMMAND "${TESTED_MODULE_EXE}"    )
			set_target_properties( "${TESTED_MODULE_EXE}" PROPERTIES RUNTIME_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/bin> )
		endif()
	endforeach()

endif()
