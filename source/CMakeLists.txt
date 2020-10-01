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
add_library( ct_resolve_hits        ${NORMSOURCES_RESOLVE_HITS}        )
add_library( ct_seq                 ${NORMSOURCES_SEQ}                 )
add_library( ct_test                ${NORMSOURCES_SRC_TEST}            )
add_library( ct_uni                 ${NORMSOURCES_UNI}                 )

target_include_directories( ct_clustagglom         PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_clustagglom>  )
target_include_directories( ct_common              PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_common>       )
target_include_directories( ct_test                PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src_test>         )

# # TODO: Try to eradicate these...
target_include_directories( ct_biocore             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_cath_assign_domains PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_cath_cluster        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_cath_refine_align   PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_cath_score_align    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_cath_superpose      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_chopping            PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_cluster             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_display_colour      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_options             PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_resolve_hits        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_seq                 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )
target_include_directories( ct_uni                 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>                  )

target_include_directories( ct_uni          SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party_code> )

target_link_libraries     ( ct_biocore             PUBLIC ct_common                                                                          )
target_link_libraries     ( ct_cath_assign_domains PUBLIC ct_common                                                                          )
target_link_libraries     ( ct_cath_cluster        PUBLIC ct_clustagglom                                                                     )
target_link_libraries     ( ct_cath_refine_align   PUBLIC ct_common ct_uni                                                                   )
target_link_libraries     ( ct_cath_score_align    PUBLIC ct_common                                                                          )
target_link_libraries     ( ct_cath_superpose      PUBLIC ct_common ct_uni                                                                   )
target_link_libraries     ( ct_chopping            PUBLIC ct_biocore Boost::program_options                                                  )
target_link_libraries     ( ct_clustagglom         PUBLIC ct_common Boost::program_options                                                   )
target_link_libraries     ( ct_cluster             PUBLIC ct_common ct_options ct_seq Boost::program_options                                 )
target_link_libraries     ( ct_common              PUBLIC Boost::boost Boost::log Boost::thread Boost::timer ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries     ( ct_display_colour      PUBLIC ct_common                                                                          )
target_link_libraries     ( ct_options             PUBLIC ct_chopping                                                                        )
target_link_libraries     ( ct_resolve_hits        PUBLIC ct_display_colour ct_options ct_seq                                                )
target_link_libraries     ( ct_seq                 PUBLIC ct_common                                                                          )
target_link_libraries     ( ct_test                PUBLIC ct_common Boost::unit_test_framework                                               )
target_link_libraries     ( ct_uni                 PUBLIC ct_display_colour ct_options Boost::iostreams Boost::serialization                 )

#####################################################################
#####################################################################

set( TESTSOURCES_REMAIN ${TESTSOURCES} )
list( REMOVE_ITEM TESTSOURCES_REMAIN
		${TESTSOURCES_EXECUTABLES}

		${TESTSOURCES_BIOCORE}
		${TESTSOURCES_CATH_ASSIGN_DOMAINS}
		${TESTSOURCES_CATH_CLUSTER}
		${TESTSOURCES_CATH_REFINE_ALIGN}
		${TESTSOURCES_CATH_SCORE_ALIGN}
		${TESTSOURCES_CATH_SUPERPOSE}
		${TESTSOURCES_CHOPPING}
		${TESTSOURCES_CLUSTER}
		${TESTSOURCES_DISPLAY_COLOUR}
		${TESTSOURCES_OPTIONS}
		${TESTSOURCES_RESOLVE_HITS}
		${TESTSOURCES_SEQ}
		${TESTSOURCES_SRC_CLUSTAGGLOM}
		${TESTSOURCES_SRC_COMMON}
		${TESTSOURCES_SRC_TEST}
		${TESTSOURCES_UNI}
)

if ( TESTSOURCES_REMAIN )
	message( SEND_ERROR "Some TESTSOURCES are not being handled : \"${TESTSOURCES_REMAIN}\"" )
endif()

add_library          ( testsrcs_biocore             OBJECT ${TESTSOURCES_BIOCORE}             )
add_library          ( testsrcs_cath_assign_domains OBJECT ${TESTSOURCES_CATH_ASSIGN_DOMAINS} )
add_library          ( testsrcs_cath_cluster        OBJECT ${TESTSOURCES_CATH_CLUSTER}        )
add_library          ( testsrcs_cath_refine_align   OBJECT ${TESTSOURCES_CATH_REFINE_ALIGN}   )
add_library          ( testsrcs_cath_score_align    OBJECT ${TESTSOURCES_CATH_SCORE_ALIGN}    )
add_library          ( testsrcs_cath_superpose      OBJECT ${TESTSOURCES_CATH_SUPERPOSE}      )
add_library          ( testsrcs_chopping            OBJECT ${TESTSOURCES_CHOPPING}            )
add_library          ( testsrcs_clustagglom         OBJECT ${TESTSOURCES_SRC_CLUSTAGGLOM}     )
add_library          ( testsrcs_cluster             OBJECT ${TESTSOURCES_CLUSTER}             )
add_library          ( testsrcs_common              OBJECT ${TESTSOURCES_SRC_COMMON}          )
add_library          ( testsrcs_display_colour      OBJECT ${TESTSOURCES_DISPLAY_COLOUR}      )
add_library          ( testsrcs_options             OBJECT ${TESTSOURCES_OPTIONS}             )
add_library          ( testsrcs_resolve_hits        OBJECT ${TESTSOURCES_RESOLVE_HITS}        )
add_library          ( testsrcs_seq                 OBJECT ${TESTSOURCES_SEQ}                 )
add_library          ( testsrcs_src_test            OBJECT ${TESTSOURCES_SRC_TEST}            )
add_library          ( testsrcs_uni                 OBJECT ${TESTSOURCES_UNI}                 )

target_link_libraries( testsrcs_biocore             PUBLIC ct_test ct_biocore             )
target_link_libraries( testsrcs_cath_assign_domains PUBLIC ct_test ct_cath_assign_domains )
target_link_libraries( testsrcs_cath_cluster        PUBLIC ct_test ct_cath_cluster        )
target_link_libraries( testsrcs_cath_refine_align   PUBLIC ct_test ct_cath_refine_align   )
target_link_libraries( testsrcs_cath_score_align    PUBLIC ct_test ct_cath_score_align    )
target_link_libraries( testsrcs_cath_superpose      PUBLIC ct_test ct_cath_superpose      )
target_link_libraries( testsrcs_chopping            PUBLIC ct_test ct_chopping            )
target_link_libraries( testsrcs_clustagglom         PUBLIC ct_test ct_clustagglom         )
target_link_libraries( testsrcs_cluster             PUBLIC ct_test ct_cluster             )
target_link_libraries( testsrcs_common              PUBLIC ct_test ct_common              )
target_link_libraries( testsrcs_display_colour      PUBLIC ct_test ct_display_colour      )
target_link_libraries( testsrcs_options             PUBLIC ct_test ct_options             )
target_link_libraries( testsrcs_resolve_hits        PUBLIC ct_test ct_resolve_hits        )
target_link_libraries( testsrcs_seq                 PUBLIC ct_test ct_seq                 )
target_link_libraries( testsrcs_src_test            PUBLIC ct_test ct_chopping            )
target_link_libraries( testsrcs_uni                 PUBLIC ct_test ct_uni                 )

#####################################################################
#####################################################################

add_executable(
	build-test

		$<TARGET_OBJECTS:testsrcs_biocore>
		$<TARGET_OBJECTS:testsrcs_cath_assign_domains>
		$<TARGET_OBJECTS:testsrcs_cath_cluster>
		$<TARGET_OBJECTS:testsrcs_cath_refine_align>
		$<TARGET_OBJECTS:testsrcs_cath_score_align>
		$<TARGET_OBJECTS:testsrcs_cath_superpose>
		$<TARGET_OBJECTS:testsrcs_chopping>
		$<TARGET_OBJECTS:testsrcs_clustagglom>
		$<TARGET_OBJECTS:testsrcs_cluster>
		$<TARGET_OBJECTS:testsrcs_common>
		$<TARGET_OBJECTS:testsrcs_display_colour>
		$<TARGET_OBJECTS:testsrcs_options>
		$<TARGET_OBJECTS:testsrcs_resolve_hits>
		$<TARGET_OBJECTS:testsrcs_seq>
		$<TARGET_OBJECTS:testsrcs_src_test>
		$<TARGET_OBJECTS:testsrcs_uni>
		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
)

target_link_libraries( build-test PRIVATE ct_cath_cluster ct_cath_refine_align ct_cath_superpose ct_cluster ct_resolve_hits ct_test Boost::filesystem )

add_test( NAME build-test COMMAND build-test )

#####################################################################
#####################################################################

add_executable( cath-assign-domains ${NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS} )
add_executable( cath-cluster        ${NORMSOURCES_EXECUTABLES_CATH_CLUSTER}        )
add_executable( cath-map-clusters   ${NORMSOURCES_EXECUTABLES_CATH_MAP_CLUSTERS}   )
add_executable( cath-refine-align   ${NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN}   )
add_executable( cath-resolve-hits   ${NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS}   )
add_executable( cath-score-align    ${NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN}    )
add_executable( cath-ssap           ${NORMSOURCES_EXECUTABLES_CATH_SSAP}           )
add_executable( cath-superpose      ${NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE}      )

target_link_libraries( cath-assign-domains PRIVATE ct_cath_assign_domains ct_uni            ct_chopping Boost::filesystem )
target_link_libraries( cath-cluster        PRIVATE ct_cath_cluster               ct_options ct_chopping                   )
target_link_libraries( cath-map-clusters   PRIVATE ct_cluster                               ct_chopping                   )
target_link_libraries( cath-refine-align   PRIVATE ct_cath_refine_align                     ct_chopping                   )
target_link_libraries( cath-resolve-hits   PRIVATE ct_resolve_hits                          ct_chopping                   )
target_link_libraries( cath-score-align    PRIVATE ct_cath_score_align    ct_uni            ct_chopping Boost::filesystem )
target_link_libraries( cath-ssap           PRIVATE                        ct_uni            ct_chopping Boost::filesystem )
target_link_libraries( cath-superpose      PRIVATE ct_cath_superpose                        ct_chopping                   )

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

option( BUILD_EXTRA_CATH_TESTS "BUILD_EXTRA_CATH_TESTS" $ENV{BUILD_EXTRA_CATH_TESTS} )

IF ( BUILD_EXTRA_CATH_TESTS )

	add_executable       ( mod-test-biocore             $<TARGET_OBJECTS:testsrcs_biocore>             ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-cath-assign-domains $<TARGET_OBJECTS:testsrcs_cath_assign_domains> ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-cath-cluster        $<TARGET_OBJECTS:testsrcs_cath_cluster>        ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-cath-refine-align   $<TARGET_OBJECTS:testsrcs_cath_refine_align>   ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-cath-score-align    $<TARGET_OBJECTS:testsrcs_cath_score_align>    ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-cath-superpose      $<TARGET_OBJECTS:testsrcs_cath_superpose>      ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-chopping            $<TARGET_OBJECTS:testsrcs_chopping>            ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-clustagglom         $<TARGET_OBJECTS:testsrcs_clustagglom>         ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-cluster             $<TARGET_OBJECTS:testsrcs_cluster>             ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-common              $<TARGET_OBJECTS:testsrcs_common>              ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-display-colour      $<TARGET_OBJECTS:testsrcs_display_colour>      ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-options             $<TARGET_OBJECTS:testsrcs_options>             ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-resolve-hits        $<TARGET_OBJECTS:testsrcs_resolve_hits>        ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-seq                 $<TARGET_OBJECTS:testsrcs_seq>                 ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-src-test            $<TARGET_OBJECTS:testsrcs_src_test>            ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )
	add_executable       ( mod-test-uni                 $<TARGET_OBJECTS:testsrcs_uni>                 ${TESTSOURCES_EXECUTABLES_BUILD_TEST} )

	target_link_libraries( mod-test-biocore             PUBLIC ct_biocore                    ct_test )
	target_link_libraries( mod-test-cath-assign-domains PUBLIC ct_cath_assign_domains        ct_test )
	target_link_libraries( mod-test-cath-cluster        PUBLIC ct_cath_cluster               ct_test )
	target_link_libraries( mod-test-cath-refine-align   PUBLIC ct_cath_refine_align          ct_test )
	target_link_libraries( mod-test-cath-score-align    PUBLIC ct_cath_score_align           ct_test )
	target_link_libraries( mod-test-cath-superpose      PUBLIC ct_cath_superpose             ct_test )
	target_link_libraries( mod-test-chopping            PUBLIC ct_chopping                   ct_test )
	target_link_libraries( mod-test-clustagglom         PUBLIC ct_clustagglom                ct_test )
	target_link_libraries( mod-test-cluster             PUBLIC ct_cluster                    ct_test )
	target_link_libraries( mod-test-common              PUBLIC ct_common                     ct_test )
	target_link_libraries( mod-test-display-colour      PUBLIC ct_display_colour             ct_test )
	target_link_libraries( mod-test-options             PUBLIC ct_options             ct_uni ct_test )
	target_link_libraries( mod-test-resolve-hits        PUBLIC ct_resolve_hits               ct_test )
	target_link_libraries( mod-test-seq                 PUBLIC ct_seq                        ct_test )
	target_link_libraries( mod-test-src-test            PUBLIC ct_src_test                   ct_test )
	target_link_libraries( mod-test-uni                 PUBLIC ct_uni                        ct_test )

	add_test             ( NAME mod-test-biocore             COMMAND mod-test-biocore             )
	add_test             ( NAME mod-test-cath-assign-domains COMMAND mod-test-cath-assign-domains )
	add_test             ( NAME mod-test-cath-cluster        COMMAND mod-test-cath-cluster        )
	add_test             ( NAME mod-test-cath-refine-align   COMMAND mod-test-cath-refine-align   )
	add_test             ( NAME mod-test-cath-score-align    COMMAND mod-test-cath-score-align    )
	add_test             ( NAME mod-test-cath-superpose      COMMAND mod-test-cath-superpose      )
	add_test             ( NAME mod-test-chopping            COMMAND mod-test-chopping            )
	add_test             ( NAME mod-test-clustagglom         COMMAND mod-test-clustagglom         )
	add_test             ( NAME mod-test-cluster             COMMAND mod-test-cluster             )
	add_test             ( NAME mod-test-common              COMMAND mod-test-common              )
	add_test             ( NAME mod-test-display-colour      COMMAND mod-test-display-colour      )
	add_test             ( NAME mod-test-options             COMMAND mod-test-options             )
	add_test             ( NAME mod-test-resolve-hits        COMMAND mod-test-resolve-hits        )
	add_test             ( NAME mod-test-seq                 COMMAND mod-test-seq                 )
	add_test             ( NAME mod-test-src-test            COMMAND mod-test-src-test            )
	add_test             ( NAME mod-test-uni                 COMMAND mod-test-uni                 )

ENDIF()

#####################################################################
#####################################################################

option (BUILD_EXTRA_CATH_TOOLS "BUILD_EXTRA_CATH_TOOLS" $ENV{BUILD_EXTRA_CATH_TOOLS} )

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

	target_link_libraries( cath-extract-pdb PRIVATE                   ct_uni ct_biocore ct_chopping                   Boost::iostreams )
	target_link_libraries( check-pdb        PRIVATE                   ct_uni ct_biocore ct_chopping Boost::filesystem Boost::iostreams )
	target_link_libraries( snap-judgement   PRIVATE ct_cath_superpose ct_uni ct_biocore ct_chopping                   Boost::iostreams )

ENDIF()

#####################################################################
#####################################################################
