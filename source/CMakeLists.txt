include( ./auto_generated_file_list.cmake )

#####################################################################

option( BUILD_EXTRA_CATH_TESTS "BUILD_EXTRA_CATH_TESTS" $ENV{BUILD_EXTRA_CATH_TESTS} )
option( BUILD_EXTRA_CATH_TOOLS "BUILD_EXTRA_CATH_TOOLS" $ENV{BUILD_EXTRA_CATH_TOOLS} )

#####################################################################
### Modules
#####################################################################

# 1. Specify the modules to build

set( MODULE_DIRS

	# TODO: Rename uni -> src_uni and move the other XXX inside a new src_XXX directory
	#
	# TODO: Then, remove the SRC_ from these and add it later where required rather than removing it
	CATH_REFINE_ALIGN
	CATH_SUPERPOSE
	CHOPPING
	DISPLAY_COLOUR
	UNI

	SRC_BIOCORE
	SRC_CATH_ASSIGN_DOMAINS
	SRC_CATH_CLUSTER
	SRC_CATH_SCORE_ALIGN
	SRC_CLUSTAGGLOM
	SRC_CLUSTER
	SRC_COMMON
	SRC_OPTIONS
	SRC_RESOLVE_HITS
	SRC_SEQ

	SRC_TEST
)

# 2. Set up ct_${MODULE} and testsrcs_${MODULE} for each MODULE

foreach( MODULE_DIR ${MODULE_DIRS} )

	# TODO: When all dirs begin src_, remove SRC_ from the above list and change this code to add it for the dir
	#       rather than removing it from the dir
	string( TOLOWER "${MODULE_DIR}" MODULE )
	string( REPLACE "src_" "" MODULE "${MODULE}" )
	list  ( APPEND MODULES "${MODULE}" )

	# Create a library for the module with all the .cpp files
	# and add the directory as a public include directory
	add_library( "ct_${MODULE}" ${NORMSOURCES_${MODULE_DIR}} )
	if ( MODULE_DIR MATCHES "^SRC_" )
		target_include_directories( "ct_${MODULE}" PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_${MODULE}>" )
	else()
		# TODO: When all the module dirs are prefixed with src_, nothing will use this else condition,
		#       so remove everything except the `target_include_directories( [...] /src_${MODULE}>" )`
		target_include_directories( "ct_${MODULE}" PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>" )
	endif()

	# Create the OBJECT library of test sources for this module and link it to the appropriate dependencies (determined above)
	add_library( "testsrcs_${MODULE}" OBJECT ${TESTSOURCES_${MODULE_DIR}} )
	target_link_libraries( "testsrcs_${MODULE}" PUBLIC ct_test "ct_${MODULE}" )

	# Record that list of test sources that are being handled (which can be checked later for omissions)
	list( APPEND ALL_MODULE_TESTSOURCES ${TESTSOURCES_${MODULE_DIR}} )
endforeach()

# 3. Specify the modules' dependencies

# TODO: When the uni directory has been rename to src_uni, eradicate this line
target_include_directories( ct_uni        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/uni>              )

# TODO: Can this be tightened to be PRIVATE?
# TODO: Wrap the third party dep in a target (perhaps with nothing but public target_include_directories)
target_include_directories( ct_uni SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party_code> )

target_link_libraries( testsrcs_test          PUBLIC ct_chopping )
target_link_libraries( testsrcs_test          PUBLIC ct_uni      ) # ct_uni for file/pdb/pdb.hpp
target_link_libraries( testsrcs_options       PUBLIC ct_uni      ) # ct_uni for alignment/options_block/alignment_input_options_block.hpp

target_link_libraries( ct_biocore             PUBLIC ct_common                                                                          )
target_link_libraries( ct_cath_assign_domains PUBLIC ct_biocore ct_options ct_uni                                                       ) # ct_uni for acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp
target_link_libraries( ct_cath_cluster        PUBLIC ct_clustagglom ct_options                                                          )
target_link_libraries( ct_cath_refine_align   PUBLIC ct_common ct_uni                                                                   ) # ct_uni for acquirer/alignment_acquirer/alignment_acquirer.hpp and alignment/options_block/alignment_input_options_block.hpp
target_link_libraries( ct_cath_score_align    PUBLIC ct_biocore ct_options ct_uni                                                       ) # ct_uni for alignment/options_block/alignment_input_options_block.hpp
target_link_libraries( ct_cath_superpose      PUBLIC ct_common ct_uni                                                                   ) # ct_uni for acquirer/alignment_acquirer/align_refining.hpp acquirer/alignment_acquirer/alignment_acquirer.hpp
target_link_libraries( ct_chopping            PUBLIC ct_biocore Boost::program_options                                                  )
target_link_libraries( ct_clustagglom         PUBLIC ct_common Boost::program_options                                                   )
target_link_libraries( ct_cluster             PUBLIC ct_common ct_options ct_seq Boost::program_options                                 )
target_link_libraries( ct_common              PUBLIC Boost::boost Boost::log Boost::thread Boost::timer ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( ct_display_colour      PUBLIC ct_common                                                                          )
target_link_libraries( ct_options             PUBLIC ct_chopping                                                                        )
target_link_libraries( ct_resolve_hits        PUBLIC ct_display_colour ct_options ct_seq                                                )
target_link_libraries( ct_seq                 PUBLIC ct_common                                                                          )
target_link_libraries( ct_test                PUBLIC ct_common Boost::unit_test_framework                                               )
target_link_libraries( ct_uni                 PUBLIC ct_display_colour ct_options Boost::iostreams Boost::serialization                 )

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
### (Non-test) executables
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

target_link_libraries( cath-assign-domains PRIVATE ct_cath_assign_domains )
target_link_libraries( cath-assign-domains PRIVATE ct_uni                 ) # ct_uni for file/prc_scores_file/prc_scores_file.hpp
target_link_libraries( cath-cluster        PRIVATE ct_cath_cluster        )
target_link_libraries( cath-map-clusters   PRIVATE ct_cluster             )
target_link_libraries( cath-refine-align   PRIVATE ct_cath_refine_align   )
target_link_libraries( cath-resolve-hits   PRIVATE ct_resolve_hits        )
target_link_libraries( cath-score-align    PRIVATE ct_cath_score_align    )
target_link_libraries( cath-ssap           PRIVATE ct_chopping            )
target_link_libraries( cath-ssap           PRIVATE ct_uni                 ) # ct_uni for ssap/options/cath_ssap_options.hpp
target_link_libraries( cath-superpose      PRIVATE ct_cath_superpose      )

IF ( BUILD_EXTRA_CATH_TOOLS )
	target_link_libraries( cath-extract-pdb PRIVATE ct_chopping )
	target_link_libraries( cath-extract-pdb PRIVATE ct_uni      ) # ct_uni for file/pdb/pdb.hpp
	target_link_libraries( check-pdb        PRIVATE ct_chopping )
	target_link_libraries( check-pdb        PRIVATE ct_uni      ) # ct_uni for file/pdb/pdb.hpp
	target_link_libraries( snap-judgement   PRIVATE ct_chopping )
	target_link_libraries( snap-judgement   PRIVATE ct_uni      ) # ct_uni for scan/scan_tools/all_vs_all.hpp
ENDIF()

#####################################################################

# Check that all the TESTSOURCES are accounted for
set( TESTSOURCES_REMAIN ${TESTSOURCES} )
list( REMOVE_ITEM TESTSOURCES_REMAIN ${TESTSOURCES_EXECUTABLES_BUILD_TEST} ${ALL_MODULE_TESTSOURCES} )
if ( TESTSOURCES_REMAIN )
	message( SEND_ERROR "Some TESTSOURCES are not being handled : \"${TESTSOURCES_REMAIN}\"" )
endif()

#####################################################################

# Add build-test, containing all the tests
add_executable( build-test ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
foreach( MODULE ${MODULES} )
	target_link_libraries( build-test PRIVATE "testsrcs_${MODULE}" )
endforeach()
add_test( NAME build-test COMMAND build-test )
set_target_properties( build-test PROPERTIES RUNTIME_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/bin> )

#####################################################################

# If requested, build a mod-test-XXX for each module
if ( BUILD_EXTRA_CATH_TESTS )

	foreach( MODULE ${MODULES} )
		# Make the test name (eg resolve_hits -> mod-test-resolve-hits)
		string               ( REPLACE "_" "-" TESTED_MODULE_EXE "mod-test-${MODULE}" )

		# Make the test, link it to the testsrcs OBJECT library and register it as a test
		add_executable       ( "${TESTED_MODULE_EXE}" ${TESTSOURCES_EXECUTABLES_BUILD_TEST}  )
		target_link_libraries( "${TESTED_MODULE_EXE}" "testsrcs_${MODULE}"                   )
		add_test             ( NAME "${TESTED_MODULE_EXE}" COMMAND "${TESTED_MODULE_EXE}"    )
		set_target_properties( "${TESTED_MODULE_EXE}" PROPERTIES RUNTIME_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/bin> )
	endforeach()

endif()
