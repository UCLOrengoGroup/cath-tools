
include( ./auto_generated_file_list.cmake )

#####################################################################
#####################################################################

include_directories( $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/uni> )

add_library( ct_biocore             ${NORMSOURCES_BIOCORE}             )
add_library( ct_cath_assign_domains ${NORMSOURCES_CATH_ASSIGN_DOMAINS} )
add_library( ct_cath_cluster        ${NORMSOURCES_CATH_CLUSTER}        )
add_library( ct_cath_refine_align   ${NORMSOURCES_CATH_REFINE_ALIGN}   )
add_library( ct_cath_score_align    ${NORMSOURCES_CATH_SCORE_ALIGN}    )
add_library( ct_cath_superpose      ${NORMSOURCES_CATH_SUPERPOSE}      )
add_library( ct_chopping            ${NORMSOURCES_CHOPPING}            )
add_library( ct_clustagglom         ${NORMSOURCES_SRC_CLUSTAGGLOM}     )
add_library( ct_cluster             ${NORMSOURCES_CLUSTER}             )
add_library( ct_common              ${NORMSOURCES_SRC_COMMON}          )
add_library( ct_display_colour      ${NORMSOURCES_DISPLAY_COLOUR}      )
add_library( ct_options             ${NORMSOURCES_OPTIONS}             )
add_library( ct_uni                 ${NORMSOURCES_UNI}                 )
add_library( ct_resolve_hits        ${NORMSOURCES_RESOLVE_HITS}        )
add_library( ct_seq                 ${NORMSOURCES_SEQ}                 )
add_library( ct_test                ${NORMSOURCES_SRC_TEST}            )

target_include_directories( ct_clustagglom PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_clustagglom> )
target_include_directories( ct_common      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common>      )
target_include_directories( ct_test        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>        )

# TODO: Try to eradicate these...
target_include_directories( ct_biocore             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_assign_domains PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_cluster        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_refine_align   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_score_align    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cath_superpose      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_chopping            PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_cluster             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_display_colour      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_options             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_uni                 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_resolve_hits        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( ct_seq                 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

target_include_directories( ct_uni                 SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party_code> )

target_link_libraries     ( ct_biocore             PUBLIC ct_common                                                        )
target_link_libraries     ( ct_cath_assign_domains PUBLIC ct_common                                                        )
target_link_libraries     ( ct_cath_cluster        PUBLIC ct_common ct_clustagglom                                         )
target_link_libraries     ( ct_cath_refine_align   PUBLIC ct_common                                                        )
target_link_libraries     ( ct_cath_score_align    PUBLIC ct_common                                                        )
target_link_libraries     ( ct_cath_superpose      PUBLIC ct_common                                                        )
target_link_libraries     ( ct_chopping            PUBLIC ct_common ct_biocore                                             )
target_link_libraries     ( ct_clustagglom         PUBLIC ct_common                                                        )
target_link_libraries     ( ct_cluster             PUBLIC ct_common                                                        )
target_link_libraries     ( ct_common              PUBLIC Boost::boost Boost::log Boost::thread Boost::timer ${RT_LIBRARY} )
target_link_libraries     ( ct_display_colour      PUBLIC ct_common                                                        )
target_link_libraries     ( ct_options             PUBLIC ct_common ct_chopping Boost::program_options                     )
target_link_libraries     ( ct_uni                 PUBLIC ct_common                                                        )
target_link_libraries     ( ct_resolve_hits        PUBLIC ct_common                                                        )
target_link_libraries     ( ct_seq                 PUBLIC ct_common                                                        )
target_link_libraries     ( ct_test                PUBLIC ct_common                                                        )

#####################################################################
#####################################################################

# TODO: Have ALL_NORMSOURCES and ALL_TESTSOURCES get auto-generated in auto_generated_file_list.cmake
set(
	ALL_TESTSOURCES
		${TESTSOURCES_BIOCORE}
		${TESTSOURCES_CATH_ASSIGN_DOMAINS}
		${TESTSOURCES_CATH_CLUSTER}
		${TESTSOURCES_CATH_REFINE_ALIGN}
		${TESTSOURCES_CATH_SCORE_ALIGN}
		${TESTSOURCES_CATH_SUPERPOSE}
		${TESTSOURCES_CHOPPING}
		${TESTSOURCES_CLUSTER}
		${TESTSOURCES_DISPLAY_COLOUR}
		${TESTSOURCES_EXECUTABLES}
		${TESTSOURCES_OPTIONS}
		${TESTSOURCES_RESOLVE_HITS}
		${TESTSOURCES_SEQ}
		${TESTSOURCES_SRC_CLUSTAGGLOM}
		${TESTSOURCES_SRC_COMMON}
		${TESTSOURCES_SRC_TEST}
		${TESTSOURCES_UNI}
)

set( TESTSOURCES_REMAIN ${ALL_TESTSOURCES} )
list( REMOVE_ITEM TESTSOURCES_REMAIN
		${TESTSOURCES_EXECUTABLES}

		${TESTSOURCES_BIOCORE}
		${TESTSOURCES_CHOPPING}
		${TESTSOURCES_CLUSTER}
		${TESTSOURCES_RESOLVE_HITS}
		${TESTSOURCES_SEQ}
		${TESTSOURCES_SRC_CLUSTAGGLOM}
		${TESTSOURCES_SRC_COMMON}
)

add_library( testsrcs_biocore      OBJECT ${TESTSOURCES_BIOCORE}         )
add_library( testsrcs_chopping     OBJECT ${TESTSOURCES_CHOPPING}        )
add_library( testsrcs_clustagglom  OBJECT ${TESTSOURCES_SRC_CLUSTAGGLOM} )
add_library( testsrcs_cluster      OBJECT ${TESTSOURCES_CLUSTER}         )
add_library( testsrcs_common       OBJECT ${TESTSOURCES_SRC_COMMON}      )
add_library( testsrcs_remain       OBJECT ${TESTSOURCES_REMAIN}          )
add_library( testsrcs_resolve_hits OBJECT ${TESTSOURCES_RESOLVE_HITS}    )
add_library( testsrcs_seq          OBJECT ${TESTSOURCES_SEQ}             )

# These should really use target_link_libraries() to their corresponding library but CMake can't seem to do that
target_include_directories( testsrcs_biocore      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>                                                                )
target_include_directories( testsrcs_common       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>                                                                )
target_include_directories( testsrcs_chopping     PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>                                                                )
target_include_directories( testsrcs_clustagglom  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_clustagglom> )
target_include_directories( testsrcs_cluster      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>                                                                )
target_include_directories( testsrcs_remain       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_clustagglom> )
target_include_directories( testsrcs_resolve_hits PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>                                                                )
target_include_directories( testsrcs_seq          PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>                                                                )

# TODO: Try to eradicate these...
target_include_directories( testsrcs_biocore      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_chopping     PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_cluster      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_common       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_remain       PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_resolve_hits PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_include_directories( testsrcs_seq          PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

target_include_directories( testsrcs_remain       SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party_code> )

#####################################################################
#####################################################################

add_executable(
	build-test
		$<TARGET_OBJECTS:testsrcs_biocore>
		$<TARGET_OBJECTS:testsrcs_chopping>
		$<TARGET_OBJECTS:testsrcs_clustagglom>
		$<TARGET_OBJECTS:testsrcs_cluster>
		$<TARGET_OBJECTS:testsrcs_common>
		$<TARGET_OBJECTS:testsrcs_remain>
		$<TARGET_OBJECTS:testsrcs_resolve_hits>
		$<TARGET_OBJECTS:testsrcs_seq>
		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
)

add_executable( cath-assign-domains ${NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS} )
add_executable( cath-cluster        ${NORMSOURCES_EXECUTABLES_CATH_CLUSTER}        )
add_executable( cath-map-clusters   ${NORMSOURCES_EXECUTABLES_CATH_MAP_CLUSTERS}   )
add_executable( cath-refine-align   ${NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN}   )
add_executable( cath-resolve-hits   ${NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS}   )
add_executable( cath-score-align    ${NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN}    )
add_executable( cath-ssap           ${NORMSOURCES_EXECUTABLES_CATH_SSAP}           )
add_executable( cath-superpose      ${NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE}      )

install(
	TARGETS
		cath-assign-domains
		cath-cluster
		cath-map-clusters
		cath-refine-align
		cath-resolve-hits
		cath-score-align
		cath-ssap
		cath-superpose
	DESTINATION
		bin
)

#####################################################################
#####################################################################

target_link_libraries( build-test          PRIVATE ct_cath_cluster ct_cath_refine_align ct_cath_superpose ct_clustagglom ct_cluster ct_resolve_hits ct_test ct_seq ct_uni ct_options ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization Boost::unit_test_framework ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-assign-domains PRIVATE ct_cath_assign_domains                                                                                          ct_uni ct_options ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-cluster        PRIVATE ct_cath_cluster                                                                                          ct_seq        ct_options            ct_chopping                                                                                                      ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-map-clusters   PRIVATE ct_cluster                                                                                               ct_seq        ct_options            ct_chopping                                                                                                      ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-refine-align   PRIVATE ct_cath_refine_align                                                                                            ct_uni ct_options ct_biocore ct_chopping ct_display_colour                   Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-resolve-hits   PRIVATE ct_resolve_hits                                                                                          ct_seq        ct_options ct_biocore ct_chopping ct_display_colour                                                                                                      )
target_link_libraries( cath-score-align    PRIVATE ct_cath_score_align                                                                                             ct_uni ct_options ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-ssap           PRIVATE                                                                                                                 ct_uni ct_options ct_biocore ct_chopping ct_display_colour Boost::filesystem Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-superpose      PRIVATE ct_cath_superpose                                                                                               ct_uni ct_options ct_biocore ct_chopping ct_display_colour                   Boost::iostreams Boost::serialization                            ${GSL_LIB_SUFFIX} )

#####################################################################
#####################################################################

option( BUILD_EXTRA_CATH_TESTS "BUILD_EXTRA_CATH_TESTS" OFF )

IF ( BUILD_EXTRA_CATH_TESTS )

	add_executable( mod-test-biocore      $<TARGET_OBJECTS:testsrcs_biocore>      ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-chopping     $<TARGET_OBJECTS:testsrcs_chopping>     ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-clustagglom  $<TARGET_OBJECTS:testsrcs_clustagglom>  ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-cluster      $<TARGET_OBJECTS:testsrcs_cluster>      ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-common       $<TARGET_OBJECTS:testsrcs_common>       ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-resolve-hits $<TARGET_OBJECTS:testsrcs_resolve_hits> ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable( mod-test-seq          $<TARGET_OBJECTS:testsrcs_seq>          ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )

	target_link_libraries( mod-test-biocore      PRIVATE ct_test                                                                            ct_biocore Boost::filesystem                                                          Boost::unit_test_framework                   )
	target_link_libraries( mod-test-chopping     PRIVATE ct_test ct_chopping                                                                ct_biocore                   Boost::program_options                                   Boost::unit_test_framework                   )
	target_link_libraries( mod-test-clustagglom  PRIVATE ct_test ct_clustagglom  ct_seq ct_options ct_chopping                              ct_biocore                                                                            Boost::unit_test_framework                   )
	target_link_libraries( mod-test-cluster      PRIVATE ct_test ct_cluster      ct_seq ct_options ct_chopping                              ct_biocore                                                                            Boost::unit_test_framework                   )
	target_link_libraries( mod-test-common       PRIVATE ct_test ct_common                                                                             Boost::filesystem                        Boost::serialization Boost::timer Boost::unit_test_framework ${GSL_LIB_SUFFIX} )
	target_link_libraries( mod-test-resolve-hits PRIVATE ct_test ct_resolve_hits ct_seq ct_chopping ct_options ct_display_colour ct_options ct_biocore Boost::filesystem                                                          Boost::unit_test_framework                   )
	target_link_libraries( mod-test-seq          PRIVATE ct_test                 ct_seq ct_options ct_chopping                              ct_biocore                                                                            Boost::unit_test_framework                   )

ENDIF()

#####################################################################
#####################################################################

option (BUILD_EXTRA_CATH_TOOLS "BUILD_EXTRA_CATH_TOOLS" OFF )

IF ( BUILD_EXTRA_CATH_TOOLS )

	add_executable( cath-extract-pdb ${NORMSOURCES_EXECUTABLES_CATH_EXTRACT_PDB} )
	add_executable( check-pdb        ${NORMSOURCES_EXECUTABLES_CATH_CHECK_PDB}   )
	add_executable( snap-judgement   ${NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT}   )

	install(
		TARGETS
			cath-extract-pdb
			check-pdb
			snap-judgement
		DESTINATION
			bin
	)

	target_link_libraries( cath-extract-pdb PRIVATE                   ct_uni ct_biocore ct_chopping ct_display_colour ct_options                                    Boost::iostreams Boost::serialization ${GSL_LIB_SUFFIX} )
	target_link_libraries( check-pdb        PRIVATE                   ct_uni ct_biocore ct_chopping ct_display_colour ct_options                  Boost::filesystem Boost::iostreams Boost::serialization ${GSL_LIB_SUFFIX} )
	target_link_libraries( snap-judgement   PRIVATE ct_cath_superpose ct_uni ct_biocore ct_chopping ct_display_colour ct_options                                    Boost::iostreams Boost::serialization ${GSL_LIB_SUFFIX} )

ENDIF()

#####################################################################
#####################################################################
