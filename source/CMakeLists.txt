
include( ./auto_generated_file_list.cmake )

#####################################################################
#####################################################################

include_directories( $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/uni> )

set(
	NORMSOURCES_REMAIN
		${NORMSOURCES_UNI_ACQUIRER}
		${NORMSOURCES_UNI_ALIGNMENT}
		${NORMSOURCES_UNI_DISPLAY}
		${NORMSOURCES_UNI_FILE}
		${NORMSOURCES_UNI_OUTPUTTER}
		${NORMSOURCES_UNI_SCAN}
		${NORMSOURCES_UNI_SCORE}
		${NORMSOURCES_UNI_SSAP}
		${NORMSOURCES_UNI_STRUCTURE}
		${NORMSOURCES_UNI_SUPERPOSITION}
)

add_library( ct_biocore             ${NORMSOURCES_BIOCORE}             )
add_library( ct_cath_assign_domains ${NORMSOURCES_CATH_ASSIGN_DOMAINS} )
add_library( ct_cath_refine_align   ${NORMSOURCES_CATH_REFINE_ALIGN}   )
add_library( ct_cath_score_align    ${NORMSOURCES_CATH_SCORE_ALIGN}    )
add_library( ct_cath_superpose      ${NORMSOURCES_CATH_SUPERPOSE}      )
add_library( ct_chopping            ${NORMSOURCES_CHOPPING}            )
add_library( ct_cluster             ${NORMSOURCES_CLUSTER}             )
add_library( ct_common              ${NORMSOURCES_SRC_COMMON}          )
add_library( ct_display_colour      ${NORMSOURCES_DISPLAY_COLOUR}      )
add_library( ct_options             ${NORMSOURCES_OPTIONS}             )
add_library( ct_remain              ${NORMSOURCES_REMAIN}              )
add_library( ct_resolve_hits        ${NORMSOURCES_RESOLVE_HITS}        )
add_library( ct_seq                 ${NORMSOURCES_SEQ}                 )
add_library( ct_test                ${NORMSOURCES_SRC_TEST}            )

target_include_directories( ct_common PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common>   )
target_include_directories( ct_test   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>     )

# TODO: Try to eradicate these...
target_include_directories( ct_biocore             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_assign_domains PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_refine_align   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_score_align    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_superpose      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_chopping            PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cluster             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_display_colour      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_options             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_remain              PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_resolve_hits        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_seq                 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

target_link_libraries     ( ct_biocore             PUBLIC ct_common                                      )
target_link_libraries     ( ct_cath_assign_domains PUBLIC ct_common                                      )
target_link_libraries     ( ct_cath_refine_align   PUBLIC ct_common                                      )
target_link_libraries     ( ct_cath_score_align    PUBLIC ct_common                                      )
target_link_libraries     ( ct_cath_superpose      PUBLIC ct_common                                      )
target_link_libraries     ( ct_chopping            PUBLIC ct_common                                      )
target_link_libraries     ( ct_cluster             PUBLIC ct_common                                      )
target_link_libraries     ( ct_common              PUBLIC Boost::boost Boost::log ${RT_LIBRARY}          )
target_link_libraries     ( ct_display_colour      PUBLIC ct_common                                      )
target_link_libraries     ( ct_options             PUBLIC ct_common Boost::program_options               )
target_link_libraries     ( ct_remain              PUBLIC ct_common                                      )
target_link_libraries     ( ct_resolve_hits        PUBLIC ct_common                                      )
target_link_libraries     ( ct_seq                 PUBLIC ct_common                                      )
target_link_libraries     ( ct_test                PUBLIC ct_common                                      )

#####################################################################
#####################################################################

set(
	TESTSOURCES_REMAIN
		${TESTSOURCES_CATH_REFINE_ALIGN}
		${TESTSOURCES_CATH_SUPERPOSE}
		${TESTSOURCES_DISPLAY_COLOUR}
		${TESTSOURCES_OPTIONS}
		${TESTSOURCES_TEST_PREDICATE}
		${TESTSOURCES_UNI_ACQUIRER}
		${TESTSOURCES_UNI_ALIGNMENT}
		${TESTSOURCES_UNI_DISPLAY}
		${TESTSOURCES_UNI_FILE}
		${TESTSOURCES_UNI_OUTPUTTER}
		${TESTSOURCES_UNI_SCAN}
		${TESTSOURCES_UNI_SCORE}
		${TESTSOURCES_UNI_SSAP}
		${TESTSOURCES_UNI_STRUCTURE}
		${TESTSOURCES_UNI_SUPERPOSITION}
)

add_library( testsrcs_biocore      OBJECT ${TESTSOURCES_BIOCORE}      )
add_library( testsrcs_chopping     OBJECT ${TESTSOURCES_CHOPPING}     )
add_library( testsrcs_cluster      OBJECT ${TESTSOURCES_CLUSTER}      )
add_library( testsrcs_common       OBJECT ${TESTSOURCES_SRC_COMMON}   )
add_library( testsrcs_remain       OBJECT ${TESTSOURCES_REMAIN}       )
add_library( testsrcs_resolve_hits OBJECT ${TESTSOURCES_RESOLVE_HITS} )
add_library( testsrcs_seq          OBJECT ${TESTSOURCES_SEQ}          )

# These should really use target_link_libraries() to their corresponding library but CMake can't seem to do that
target_include_directories( testsrcs_biocore      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> )
target_include_directories( testsrcs_common       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> )
target_include_directories( testsrcs_chopping     PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> )
target_include_directories( testsrcs_cluster      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> )
target_include_directories( testsrcs_remain       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> )
target_include_directories( testsrcs_resolve_hits PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> )
target_include_directories( testsrcs_seq          PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> )

# TODO: Try to eradicate these...
target_include_directories( testsrcs_biocore      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_chopping     PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_cluster      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_common       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_remain       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_resolve_hits PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_seq          PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

#####################################################################
#####################################################################

add_executable(
	build-test
		$<TARGET_OBJECTS:testsrcs_biocore>
		$<TARGET_OBJECTS:testsrcs_chopping>
		$<TARGET_OBJECTS:testsrcs_cluster>
		$<TARGET_OBJECTS:testsrcs_common>
		$<TARGET_OBJECTS:testsrcs_remain>
		$<TARGET_OBJECTS:testsrcs_resolve_hits>
		$<TARGET_OBJECTS:testsrcs_seq>
		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
)

add_executable( cath-assign-domains ${NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS} )
add_executable( cath-map-clusters   ${NORMSOURCES_EXECUTABLES_CATH_MAP_CLUSTERS}   )
add_executable( cath-refine-align   ${NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN}   )
add_executable( cath-resolve-hits   ${NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS}   )
add_executable( cath-score-align    ${NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN}    )
add_executable( cath-ssap           ${NORMSOURCES_EXECUTABLES_CATH_SSAP}           )
add_executable( cath-superpose      ${NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE}      )

#####################################################################
#####################################################################

target_link_libraries( build-test          PRIVATE ct_cath_refine_align ct_cath_superpose ct_cluster ct_resolve_hits ct_test ct_seq ct_options ct_remain ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization Boost::unit_test_framework ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-assign-domains PRIVATE ct_cath_assign_domains                                                           ct_options ct_remain ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-map-clusters   PRIVATE ct_cluster                                                                ct_seq ct_options                                                                                                                                       ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-refine-align   PRIVATE ct_cath_refine_align                                                             ct_options ct_remain ct_biocore ct_chopping ct_display_colour                   Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-resolve-hits   PRIVATE ct_resolve_hits                                                           ct_seq ct_options           ct_biocore             ct_display_colour                                                                                                      )
target_link_libraries( cath-score-align    PRIVATE ct_cath_score_align                                                              ct_options ct_remain ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-ssap           PRIVATE                                                                                  ct_options ct_remain ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-superpose      PRIVATE ct_cath_superpose                                                                ct_options ct_remain ct_biocore ct_chopping ct_display_colour                   Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )

#####################################################################
#####################################################################

option( BUILD_EXTRA_CATH_TESTS "BUILD_EXTRA_CATH_TESTS" OFF )

IF ( BUILD_EXTRA_CATH_TESTS )

	add_executable( mod-test-biocore      $<TARGET_OBJECTS:testsrcs_biocore>      ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-chopping     $<TARGET_OBJECTS:testsrcs_chopping>     ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-cluster      $<TARGET_OBJECTS:testsrcs_cluster>      ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-common       $<TARGET_OBJECTS:testsrcs_common>       ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-resolve-hits $<TARGET_OBJECTS:testsrcs_resolve_hits> ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-seq          $<TARGET_OBJECTS:testsrcs_seq>          ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )

	target_link_libraries( mod-test-biocore      PRIVATE ct_test                                                     ct_biocore Boost::filesystem                                             Boost::unit_test_framework                   )
	target_link_libraries( mod-test-chopping     PRIVATE ct_test ct_chopping                                         ct_biocore                   Boost::program_options                      Boost::unit_test_framework                   )
	target_link_libraries( mod-test-cluster      PRIVATE ct_test ct_cluster      ct_seq                   ct_options ct_biocore                                                               Boost::unit_test_framework                   )
	target_link_libraries( mod-test-common       PRIVATE ct_test ct_common                                                      Boost::filesystem                        Boost::serialization Boost::unit_test_framework ${GSL_LIB_SUFFIX} )
	target_link_libraries( mod-test-resolve-hits PRIVATE ct_test ct_resolve_hits ct_seq ct_display_colour ct_options ct_biocore Boost::filesystem                                             Boost::unit_test_framework                   )
	target_link_libraries( mod-test-seq          PRIVATE ct_test                 ct_seq                              ct_biocore                                                               Boost::unit_test_framework                   )

ENDIF()

#####################################################################
#####################################################################

option (BUILD_EXTRA_CATH_TOOLS "BUILD_EXTRA_CATH_TOOLS" OFF )

IF ( BUILD_EXTRA_CATH_TOOLS )

	add_executable( snap-judgement ${NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT} )
	add_executable( check-pdb      ${NORMSOURCES_EXECUTABLES_CATH_CHECK_PDB} )

	target_link_libraries( check-pdb      PRIVATE                   ct_remain ct_biocore ct_chopping ct_display_colour ct_options Boost::filesystem Boost::iostreams Boost::serialization ${GSL_LIB_SUFFIX} )
	target_link_libraries( snap-judgement PRIVATE ct_cath_superpose ct_remain ct_biocore ct_chopping ct_display_colour ct_options                   Boost::iostreams Boost::serialization ${GSL_LIB_SUFFIX} )

ENDIF()

#####################################################################
#####################################################################
