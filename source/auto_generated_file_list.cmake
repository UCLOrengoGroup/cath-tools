##### DON'T EDIT THIS FILE - IT'S AUTO-GENERATED #####
set(
	NORMSOURCES_BIOCORE
		biocore/chain_label.cpp
		biocore/residue_id.cpp
		biocore/residue_name.cpp
)

set(
	NORMSOURCES_CATH_ASSIGN_DOMAINS_OPTIONS
		cath_assign_domains/options/cath_assign_domains_options.cpp
		cath_assign_domains/options/cath_assign_domains_options_block.cpp
)

set(
	NORMSOURCES_CATH_ASSIGN_DOMAINS
		${NORMSOURCES_CATH_ASSIGN_DOMAINS_OPTIONS}
)

set(
	NORMSOURCES_CATH_REFINE_ALIGN_OPTIONS
		cath_refine_align/options/cath_refine_align_options.cpp
)

set(
	NORMSOURCES_CATH_REFINE_ALIGN
		cath_refine_align/cath_align_refiner.cpp
		${NORMSOURCES_CATH_REFINE_ALIGN_OPTIONS}
)

set(
	NORMSOURCES_CATH_SCORE_ALIGN_OPTIONS
		cath_score_align/options/cath_score_align_options.cpp
)

set(
	NORMSOURCES_CATH_SCORE_ALIGN
		cath_score_align/cath_align_scorer.cpp
		${NORMSOURCES_CATH_SCORE_ALIGN_OPTIONS}
)

set(
	NORMSOURCES_CATH_SUPERPOSE_OPTIONS
		cath_superpose/options/cath_superpose_options.cpp
)

set(
	NORMSOURCES_CATH_SUPERPOSE
		cath_superpose/cath_superposer.cpp
		${NORMSOURCES_CATH_SUPERPOSE_OPTIONS}
)

set(
	NORMSOURCES_CHOPPING_CHOPPING_FORMAT
		chopping/chopping_format/chopping_format.cpp
		chopping/chopping_format/domall_chopping_format.cpp
		chopping/chopping_format/jmol_selection_chopping_format.cpp
		chopping/chopping_format/scop_chopping_format.cpp
		chopping/chopping_format/sillitoe_chopping_format.cpp
		chopping/chopping_format/simple_chopping_format.cpp
)

set(
	NORMSOURCES_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER
		chopping/chopping_io/region_io/region_reader/region_reader.cpp
		chopping/chopping_io/region_io/region_reader/std_region_reader.cpp
)

set(
	NORMSOURCES_CHOPPING_CHOPPING_IO_REGION_IO_REGION_WRITER
		chopping/chopping_io/region_io/region_writer/region_writer.cpp
		chopping/chopping_io/region_io/region_writer/std_region_writer.cpp
)

set(
	NORMSOURCES_CHOPPING_CHOPPING_IO_REGION_IO
		${NORMSOURCES_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER}
		${NORMSOURCES_CHOPPING_CHOPPING_IO_REGION_IO_REGION_WRITER}
		chopping/chopping_io/region_io/std_region_io_spec.cpp
)

set(
	NORMSOURCES_CHOPPING_CHOPPING_IO
		${NORMSOURCES_CHOPPING_CHOPPING_IO_REGION_IO}
)

set(
	NORMSOURCES_CHOPPING_DOMAIN
		chopping/domain/domain.cpp
		chopping/domain/domain_definition.cpp
)

set(
	NORMSOURCES_CHOPPING_REGION
		chopping/region/region.cpp
		chopping/region/regions_limiter.cpp
)

set(
	NORMSOURCES_CHOPPING_RESIDUE_LOCATION
		chopping/residue_location/residue_locating.cpp
		chopping/residue_location/residue_location.cpp
)

set(
	NORMSOURCES_CHOPPING
		chopping/chopping.cpp
		${NORMSOURCES_CHOPPING_CHOPPING_FORMAT}
		${NORMSOURCES_CHOPPING_CHOPPING_IO}
		${NORMSOURCES_CHOPPING_DOMAIN}
		${NORMSOURCES_CHOPPING_REGION}
		${NORMSOURCES_CHOPPING_RESIDUE_LOCATION}
)

set(
	NORMSOURCES_CLUSTER_DETAIL
		cluster/detail/mapping_job.cpp
)

set(
	NORMSOURCES_CLUSTER_FILE
		cluster/file/cluster_membership_file.cpp
)

set(
	NORMSOURCES_CLUSTER_MAP
		cluster/map/aggregate_map_results.cpp
		cluster/map/map_clusters.cpp
		cluster/map/map_results.cpp
		cluster/map/overlap_frac_distn.cpp
)

set(
	NORMSOURCES_CLUSTER_OPTIONS_OPTIONS_BLOCK
		cluster/options/options_block/clust_mapping_options_block.cpp
		cluster/options/options_block/clustmap_input_options_block.cpp
		cluster/options/options_block/clustmap_output_options_block.cpp
)

set(
	NORMSOURCES_CLUSTER_OPTIONS_SPEC
		cluster/options/spec/clust_mapping_spec.cpp
		cluster/options/spec/clustmap_input_spec.cpp
		cluster/options/spec/clustmap_output_spec.cpp
)

set(
	NORMSOURCES_CLUSTER_OPTIONS
		${NORMSOURCES_CLUSTER_OPTIONS_OPTIONS_BLOCK}
		${NORMSOURCES_CLUSTER_OPTIONS_SPEC}
)

set(
	NORMSOURCES_CLUSTER
		cluster/cath_cluster_mapper.cpp
		cluster/cluster_domains.cpp
		cluster/clustmap_options.cpp
		${NORMSOURCES_CLUSTER_DETAIL}
		cluster/domain_cluster_ids.cpp
		${NORMSOURCES_CLUSTER_FILE}
		${NORMSOURCES_CLUSTER_MAP}
		cluster/new_cluster_data.cpp
		cluster/old_cluster_data.cpp
		${NORMSOURCES_CLUSTER_OPTIONS}
)

set(
	NORMSOURCES_DISPLAY_COLOUR
		display_colour/display_colour.cpp
		display_colour/display_colour_gradient.cpp
		display_colour/display_colour_list.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS
		executables/cath_assign_domains/cath_assign_domains.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_CHECK_PDB
		executables/cath_check_pdb/cath_check_pdb.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_EXTRACT_PDB
		executables/cath_extract_pdb/cath_extract_pdb.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_MAP_CLUSTERS
		executables/cath_map_clusters/cath_map_clusters.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN
		executables/cath_refine_align/cath_refine_align.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS
		executables/cath_resolve_hits/cath_resolve_hits.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN
		executables/cath_score_align/cath_score_align.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_SSAP
		executables/cath_ssap/cath_ssap.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE
		executables/cath_superpose/cath_superpose.cpp
)

set(
	NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT
		executables/snap_judgement/snap_judgement.cpp
)

set(
	NORMSOURCES_EXECUTABLES
		${NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS}
		${NORMSOURCES_EXECUTABLES_CATH_CHECK_PDB}
		${NORMSOURCES_EXECUTABLES_CATH_EXTRACT_PDB}
		${NORMSOURCES_EXECUTABLES_CATH_MAP_CLUSTERS}
		${NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN}
		${NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS}
		${NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN}
		${NORMSOURCES_EXECUTABLES_CATH_SSAP}
		${NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE}
		${NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT}
)

set(
	NORMSOURCES_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS
		options/executable/cath_check_pdb_options/cath_check_pdb_options.cpp
)

set(
	NORMSOURCES_OPTIONS_EXECUTABLE_CATH_EXTRACT_PDB_OPTIONS
		options/executable/cath_extract_pdb_options/cath_extract_pdb_options.cpp
)

set(
	NORMSOURCES_OPTIONS_EXECUTABLE
		${NORMSOURCES_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS}
		${NORMSOURCES_OPTIONS_EXECUTABLE_CATH_EXTRACT_PDB_OPTIONS}
		options/executable/env_var_option_name_handler.cpp
		options/executable/executable_options.cpp
)

set(
	NORMSOURCES_OPTIONS_OPTIONS_BLOCK
		options/options_block/alignment_input_options_block.cpp
		options/options_block/alignment_input_spec.cpp
		options/options_block/check_pdb_options_block.cpp
		options/options_block/detail_help_options_block.cpp
		options/options_block/extract_pdb_options_block.cpp
		options/options_block/ids_options_block.cpp
		options/options_block/misc_help_version_options_block.cpp
		options/options_block/options_block.cpp
		options/options_block/options_block_tester.cpp
		options/options_block/pdb_input_options_block.cpp
		options/options_block/pdb_input_spec.cpp
		options/options_block/string_options_block.cpp
		options/options_block/superposition_input_options_block.cpp
)

set(
	NORMSOURCES_OPTIONS
		${NORMSOURCES_OPTIONS_EXECUTABLE}
		${NORMSOURCES_OPTIONS_OPTIONS_BLOCK}
)

set(
	NORMSOURCES_RESOLVE_HITS_ALGO
		resolve_hits/algo/discont_hits_index_by_start.cpp
		resolve_hits/algo/masked_bests_cacher.cpp
		resolve_hits/algo/scored_arch_proxy.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_FILE_DETAIL
		resolve_hits/file/detail/hmmer_aln.cpp
		resolve_hits/file/detail/hmmer_parser.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_FILE
		resolve_hits/file/alnd_rgn.cpp
		resolve_hits/file/cath_id_score_category.cpp
		${NORMSOURCES_RESOLVE_HITS_FILE_DETAIL}
		resolve_hits/file/hits_input_format_tag.cpp
		resolve_hits/file/parse_domain_hits_table.cpp
		resolve_hits/file/parse_hmmer_out.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_HTML_OUTPUT
		resolve_hits/html_output/html_segment.cpp
		resolve_hits/html_output/resolve_hits_html_outputter.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK
		resolve_hits/options/options_block/crh_filter_options_block.cpp
		resolve_hits/options/options_block/crh_html_options_block.cpp
		resolve_hits/options/options_block/crh_input_options_block.cpp
		resolve_hits/options/options_block/crh_output_options_block.cpp
		resolve_hits/options/options_block/crh_score_options_block.cpp
		resolve_hits/options/options_block/crh_segment_options_block.cpp
		resolve_hits/options/options_block/crh_single_output_options_block.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_OPTIONS_SPEC
		resolve_hits/options/spec/crh_filter_spec.cpp
		resolve_hits/options/spec/crh_html_spec.cpp
		resolve_hits/options/spec/crh_input_spec.cpp
		resolve_hits/options/spec/crh_output_spec.cpp
		resolve_hits/options/spec/crh_score_spec.cpp
		resolve_hits/options/spec/crh_segment_spec.cpp
		resolve_hits/options/spec/crh_single_output_spec.cpp
		resolve_hits/options/spec/crh_spec.cpp
		resolve_hits/options/spec/hit_boundary_output.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_OPTIONS
		resolve_hits/options/crh_options.cpp
		${NORMSOURCES_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK}
		${NORMSOURCES_RESOLVE_HITS_OPTIONS_SPEC}
)

set(
	NORMSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR
		resolve_hits/read_and_process_hits/hits_processor/gather_hits_processor.cpp
		resolve_hits/read_and_process_hits/hits_processor/hits_processor_list.cpp
		resolve_hits/read_and_process_hits/hits_processor/summarise_hits_processor.cpp
		resolve_hits/read_and_process_hits/hits_processor/write_html_hits_processor.cpp
		resolve_hits/read_and_process_hits/hits_processor/write_json_hits_processor.cpp
		resolve_hits/read_and_process_hits/hits_processor/write_results_hits_processor.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS
		${NORMSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR}
		resolve_hits/read_and_process_hits/read_and_process_mgr.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_RESOLVE
		resolve_hits/resolve/hit_resolver.cpp
		resolve_hits/resolve/naive_greedy_hit_resolver.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_TRIM
		resolve_hits/trim/seq_seg_boundary_fns.cpp
		resolve_hits/trim/trim_spec.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS
		${NORMSOURCES_RESOLVE_HITS_ALGO}
		resolve_hits/calc_hit.cpp
		resolve_hits/calc_hit_list.cpp
		resolve_hits/cath_hit_resolver.cpp
		${NORMSOURCES_RESOLVE_HITS_FILE}
		resolve_hits/full_hit.cpp
		resolve_hits/full_hit_fns.cpp
		resolve_hits/full_hit_list.cpp
		resolve_hits/full_hit_list_fns.cpp
		resolve_hits/hit_arch.cpp
		resolve_hits/hit_extras.cpp
		resolve_hits/hit_score_type.cpp
		${NORMSOURCES_RESOLVE_HITS_HTML_OUTPUT}
		${NORMSOURCES_RESOLVE_HITS_OPTIONS}
		${NORMSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS}
		${NORMSOURCES_RESOLVE_HITS_RESOLVE}
		resolve_hits/scored_hit_arch.cpp
		${NORMSOURCES_RESOLVE_HITS_TRIM}
)

set(
	NORMSOURCES_SEQ
		seq/seq_arrow.cpp
		seq/seq_seg.cpp
		seq/seq_seg_run.cpp
)

set(
	NORMSOURCES_SRC_COMMON_COMMON_ALGORITHM
		src_common/common/algorithm/random_split.cpp
)

set(
	NORMSOURCES_SRC_COMMON_COMMON_BATCH
		src_common/common/batch/batch_functions.cpp
)

set(
	NORMSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_GRAPH
		src_common/common/boost_addenda/graph/spanning_tree.cpp
)

set(
	NORMSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_LOG
		src_common/common/boost_addenda/log/log_to_ostream_guard.cpp
)

set(
	NORMSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA
		${NORMSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_GRAPH}
		${NORMSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_LOG}
)

set(
	NORMSOURCES_SRC_COMMON_COMMON_EXCEPTION
		src_common/common/exception/invalid_argument_exception.cpp
		src_common/common/exception/not_implemented_exception.cpp
		src_common/common/exception/out_of_range_exception.cpp
		src_common/common/exception/runtime_error_exception.cpp
)

set(
	NORMSOURCES_SRC_COMMON_COMMON_FILE
		src_common/common/file/find_file.cpp
		src_common/common/file/ofstream_list.cpp
		src_common/common/file/open_fstream.cpp
		src_common/common/file/path_or_istream.cpp
		src_common/common/file/temp_file.cpp
)

set(
	NORMSOURCES_SRC_COMMON_COMMON
		${NORMSOURCES_SRC_COMMON_COMMON_ALGORITHM}
		src_common/common/argc_argv_faker.cpp
		${NORMSOURCES_SRC_COMMON_COMMON_BATCH}
		${NORMSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA}
		src_common/common/command_executer.cpp
		${NORMSOURCES_SRC_COMMON_COMMON_EXCEPTION}
		${NORMSOURCES_SRC_COMMON_COMMON_FILE}
		src_common/common/invert_permutation.cpp
		src_common/common/logger.cpp
		src_common/common/program_exception_wrapper.cpp
)

set(
	NORMSOURCES_SRC_COMMON
		${NORMSOURCES_SRC_COMMON_COMMON}
)

set(
	NORMSOURCES_SRC_TEST_TEST_PREDICATE_DETAIL
		src_test/test/predicate/detail/strings_equal.cpp
)

set(
	NORMSOURCES_SRC_TEST_TEST_PREDICATE
		src_test/test/predicate/bootstrap_mode.cpp
		${NORMSOURCES_SRC_TEST_TEST_PREDICATE_DETAIL}
		src_test/test/predicate/files_equal.cpp
		src_test/test/predicate/istream_and_file_equal.cpp
		src_test/test/predicate/istreams_equal.cpp
		src_test/test/predicate/string_matches_file.cpp
)

set(
	NORMSOURCES_SRC_TEST_TEST
		src_test/test/global_test_constants.cpp
		${NORMSOURCES_SRC_TEST_TEST_PREDICATE}
)

set(
	NORMSOURCES_SRC_TEST
		${NORMSOURCES_SRC_TEST_TEST}
)

set(
	NORMSOURCES_UNI_ACQUIRER_ALIGNMENT_ACQUIRER
		uni/acquirer/alignment_acquirer/alignment_acquirer.cpp
		uni/acquirer/alignment_acquirer/cora_aln_file_alignment_acquirer.cpp
		uni/acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer.cpp
		uni/acquirer/alignment_acquirer/residue_name_alignment_acquirer.cpp
		uni/acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer.cpp
		uni/acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.cpp
)

set(
	NORMSOURCES_UNI_ACQUIRER_PDBS_ACQUIRER
		uni/acquirer/pdbs_acquirer/domain_defn_pdbs_acquirer.cpp
		uni/acquirer/pdbs_acquirer/file_list_pdbs_acquirer.cpp
		uni/acquirer/pdbs_acquirer/istream_pdbs_acquirer.cpp
		uni/acquirer/pdbs_acquirer/pdbs_acquirer.cpp
)

set(
	NORMSOURCES_UNI_ACQUIRER_SELECTION_POLICY_ACQUIRER
		uni/acquirer/selection_policy_acquirer/selection_policy_acquirer.cpp
)

set(
	NORMSOURCES_UNI_ACQUIRER_SUPERPOSITION_ACQUIRER
		uni/acquirer/superposition_acquirer/align_based_superposition_acquirer.cpp
		uni/acquirer/superposition_acquirer/superposition_acquirer.cpp
)

set(
	NORMSOURCES_UNI_ACQUIRER
		${NORMSOURCES_UNI_ACQUIRER_ALIGNMENT_ACQUIRER}
		${NORMSOURCES_UNI_ACQUIRER_PDBS_ACQUIRER}
		${NORMSOURCES_UNI_ACQUIRER_SELECTION_POLICY_ACQUIRER}
		${NORMSOURCES_UNI_ACQUIRER_SUPERPOSITION_ACQUIRER}
)

set(
	NORMSOURCES_UNI_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY
		uni/alignment/common_atom_selection_policy/common_atom_select_ca_policy.cpp
		uni/alignment/common_atom_selection_policy/common_atom_select_cb_policy.cpp
		uni/alignment/common_atom_selection_policy/common_atom_selection_policy.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY
		uni/alignment/common_residue_selection_policy/common_residue_score_based_selection_policy.cpp
		uni/alignment/common_residue_selection_policy/common_residue_select_all_policy.cpp
		uni/alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.cpp
		uni/alignment/common_residue_selection_policy/common_residue_select_min_score_policy.cpp
		uni/alignment/common_residue_selection_policy/common_residue_selection_policy.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_DETAIL
		uni/alignment/detail/multi_align_builder.cpp
		uni/alignment/detail/multi_align_group.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER
		uni/alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.cpp
		uni/alignment/dyn_prog_align/detail/matrix_plotter/matrix_plotter.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER
		uni/alignment/dyn_prog_align/detail/string_aligner/benchmark_dyn_prog_string_aligner.cpp
		uni/alignment/dyn_prog_align/detail/string_aligner/gen_dyn_prog_string_aligner.cpp
		uni/alignment/dyn_prog_align/detail/string_aligner/string_aligner.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL
		${NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER}
		uni/alignment/dyn_prog_align/detail/path_step.cpp
		uni/alignment/dyn_prog_align/detail/return_path_matrix.cpp
		uni/alignment/dyn_prog_align/detail/score_accumulation_matrix.cpp
		${NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER}
)

set(
	NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE
		uni/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/entry_querier_dyn_prog_score_source.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN
		${NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL}
		uni/alignment/dyn_prog_align/dyn_prog_aligner.cpp
		${NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE}
		uni/alignment/dyn_prog_align/ssap_code_dyn_prog_aligner.cpp
		uni/alignment/dyn_prog_align/std_dyn_prog_aligner.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_GAP
		uni/alignment/gap/alignment_gap.cpp
		uni/alignment/gap/gap_penalty.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_IO_OUTPUTTER
		uni/alignment/io/outputter/horiz_align_outputter.cpp
		uni/alignment/io/outputter/html_align_outputter.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_IO
		uni/alignment/io/align_scaffold.cpp
		uni/alignment/io/alignment_io.cpp
		${NORMSOURCES_UNI_ALIGNMENT_IO_OUTPUTTER}
)

set(
	NORMSOURCES_UNI_ALIGNMENT_REFINER_DETAIL
		uni/alignment/refiner/detail/alignment_split.cpp
		uni/alignment/refiner/detail/alignment_split_list.cpp
		uni/alignment/refiner/detail/alignment_split_mapping.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_REFINER
		uni/alignment/refiner/alignment_refiner.cpp
		${NORMSOURCES_UNI_ALIGNMENT_REFINER_DETAIL}
		uni/alignment/refiner/indexed_refiner.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL
		uni/alignment/residue_name_align/detail/residue_name_align_map.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN
		${NORMSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL}
		uni/alignment/residue_name_align/residue_name_aligner.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_RESIDUE_SCORE
		uni/alignment/residue_score/alignment_residue_scores.cpp
		uni/alignment/residue_score/residue_scorer.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT_TOOLS
		uni/alignment/tools/alignment_breaks.cpp
)

set(
	NORMSOURCES_UNI_ALIGNMENT
		uni/alignment/alignment.cpp
		uni/alignment/alignment_action.cpp
		uni/alignment/alignment_context.cpp
		uni/alignment/alignment_coord_extractor.cpp
		uni/alignment/alignment_row.cpp
		${NORMSOURCES_UNI_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY}
		${NORMSOURCES_UNI_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY}
		${NORMSOURCES_UNI_ALIGNMENT_DETAIL}
		${NORMSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN}
		${NORMSOURCES_UNI_ALIGNMENT_GAP}
		${NORMSOURCES_UNI_ALIGNMENT_IO}
		uni/alignment/pair_alignment.cpp
		${NORMSOURCES_UNI_ALIGNMENT_REFINER}
		${NORMSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN}
		${NORMSOURCES_UNI_ALIGNMENT_RESIDUE_SCORE}
		${NORMSOURCES_UNI_ALIGNMENT_TOOLS}
)

set(
	NORMSOURCES_UNI_DISPLAY_DISPLAY_COLOUR_SPEC
		uni/display/display_colour_spec/broad_display_colour_spec.cpp
		uni/display/display_colour_spec/display_colour_spec.cpp
)

set(
	NORMSOURCES_UNI_DISPLAY_DISPLAY_COLOURER_DETAIL
		uni/display/display_colourer/detail/score_colour_handler.cpp
)

set(
	NORMSOURCES_UNI_DISPLAY_DISPLAY_COLOURER
		uni/display/display_colourer/alignment_free_display_colourer.cpp
		${NORMSOURCES_UNI_DISPLAY_DISPLAY_COLOURER_DETAIL}
		uni/display/display_colourer/display_colourer.cpp
		uni/display/display_colourer/display_colourer_alignment.cpp
		uni/display/display_colourer/display_colourer_consecutive.cpp
		uni/display/display_colourer/display_colourer_score.cpp
)

set(
	NORMSOURCES_UNI_DISPLAY_OPTIONS
		uni/display/options/display_options_block.cpp
		uni/display/options/display_spec.cpp
)

set(
	NORMSOURCES_UNI_DISPLAY_VIEWER_PYMOL
		uni/display/viewer/pymol/pymol_tools.cpp
)

set(
	NORMSOURCES_UNI_DISPLAY_VIEWER
		uni/display/viewer/chimera_viewer.cpp
		uni/display/viewer/jmol_viewer.cpp
		${NORMSOURCES_UNI_DISPLAY_VIEWER_PYMOL}
		uni/display/viewer/pymol_viewer.cpp
		uni/display/viewer/rasmol_style_viewer.cpp
		uni/display/viewer/rasmol_viewer.cpp
		uni/display/viewer/viewer.cpp
)

set(
	NORMSOURCES_UNI_DISPLAY
		${NORMSOURCES_UNI_DISPLAY_DISPLAY_COLOUR_SPEC}
		${NORMSOURCES_UNI_DISPLAY_DISPLAY_COLOURER}
		${NORMSOURCES_UNI_DISPLAY_OPTIONS}
		${NORMSOURCES_UNI_DISPLAY_VIEWER}
)

set(
	NORMSOURCES_UNI_FILE_DOMAIN_DEFINITION_LIST
		uni/file/domain_definition_list/domain_definition_list.cpp
)

set(
	NORMSOURCES_UNI_FILE_DSSP_WOLF
		uni/file/dssp_wolf/dssp_file.cpp
		uni/file/dssp_wolf/dssp_file_io.cpp
		uni/file/dssp_wolf/tally_residue_ids.cpp
		uni/file/dssp_wolf/wolf_file.cpp
		uni/file/dssp_wolf/wolf_file_io.cpp
)

set(
	NORMSOURCES_UNI_FILE_HMMER_SCORES_FILE
		uni/file/hmmer_scores_file/hmmer_scores_entry.cpp
		uni/file/hmmer_scores_file/hmmer_scores_file.cpp
)

set(
	NORMSOURCES_UNI_FILE_NAME_SET
		uni/file/name_set/name_set.cpp
		uni/file/name_set/name_set_list.cpp
)

set(
	NORMSOURCES_UNI_FILE_OPTIONS
		uni/file/options/data_dirs_options_block.cpp
		uni/file/options/data_dirs_spec.cpp
		uni/file/options/data_option.cpp
)

set(
	NORMSOURCES_UNI_FILE_PDB
		uni/file/pdb/coarse_element_type.cpp
		uni/file/pdb/dssp_skip_policy.cpp
		uni/file/pdb/pdb.cpp
		uni/file/pdb/pdb_atom.cpp
		uni/file/pdb/pdb_atom_parse_status.cpp
		uni/file/pdb/pdb_list.cpp
		uni/file/pdb/pdb_record.cpp
		uni/file/pdb/pdb_residue.cpp
		uni/file/pdb/proximity_calculator.cpp
		uni/file/pdb/read_domain_def_from_pdb.cpp
)

set(
	NORMSOURCES_UNI_FILE_PRC_SCORES_FILE
		uni/file/prc_scores_file/prc_scores_entry.cpp
		uni/file/prc_scores_file/prc_scores_file.cpp
)

set(
	NORMSOURCES_UNI_FILE_SEC
		uni/file/sec/sec_file.cpp
		uni/file/sec/sec_file_io.cpp
		uni/file/sec/sec_file_record.cpp
)

set(
	NORMSOURCES_UNI_FILE_SSAP_SCORES_FILE
		uni/file/ssap_scores_file/ssap_scores_entry.cpp
		uni/file/ssap_scores_file/ssap_scores_file.cpp
)

set(
	NORMSOURCES_UNI_FILE
		uni/file/data_file.cpp
		${NORMSOURCES_UNI_FILE_DOMAIN_DEFINITION_LIST}
		${NORMSOURCES_UNI_FILE_DSSP_WOLF}
		${NORMSOURCES_UNI_FILE_HMMER_SCORES_FILE}
		${NORMSOURCES_UNI_FILE_NAME_SET}
		${NORMSOURCES_UNI_FILE_OPTIONS}
		${NORMSOURCES_UNI_FILE_PDB}
		${NORMSOURCES_UNI_FILE_PRC_SCORES_FILE}
		${NORMSOURCES_UNI_FILE_SEC}
		${NORMSOURCES_UNI_FILE_SSAP_SCORES_FILE}
		uni/file/strucs_context.cpp
)

set(
	NORMSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER
		uni/outputter/alignment_outputter/alignment_outputter.cpp
		uni/outputter/alignment_outputter/alignment_outputter_list.cpp
		uni/outputter/alignment_outputter/cath_aln_ostream_alignment_outputter.cpp
		uni/outputter/alignment_outputter/fasta_ostream_alignment_outputter.cpp
		uni/outputter/alignment_outputter/file_alignment_outputter.cpp
		uni/outputter/alignment_outputter/html_ostream_alignment_outputter.cpp
		uni/outputter/alignment_outputter/ssap_ostream_alignment_outputter.cpp
)

set(
	NORMSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS
		uni/outputter/alignment_outputter_options/alignment_output_options_block.cpp
)

set(
	NORMSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS
		uni/outputter/superposition_output_options/superposition_output_options_block.cpp
)

set(
	NORMSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUTTER
		uni/outputter/superposition_outputter/json_file_superposition_outputter.cpp
		uni/outputter/superposition_outputter/ostream_superposition_outputter.cpp
		uni/outputter/superposition_outputter/pdb_file_superposition_outputter.cpp
		uni/outputter/superposition_outputter/pdb_files_superposition_outputter.cpp
		uni/outputter/superposition_outputter/pymol_file_superposition_outputter.cpp
		uni/outputter/superposition_outputter/pymol_view_superposition_outputter.cpp
		uni/outputter/superposition_outputter/superposition_outputter.cpp
		uni/outputter/superposition_outputter/superposition_outputter_list.cpp
)

set(
	NORMSOURCES_UNI_OUTPUTTER
		${NORMSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER}
		${NORMSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS}
		${NORMSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS}
		${NORMSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUTTER}
)

set(
	NORMSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY
		uni/scan/detail/check_scan/test_only/alignment_scan_comparison.cpp
		uni/scan/detail/check_scan/test_only/check_scan_on_final_alignment.cpp
		uni/scan/detail/check_scan/test_only/quad_and_rep_criteria_result.cpp
		uni/scan/detail/check_scan/test_only/quad_criteria_result.cpp
)

set(
	NORMSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN
		${NORMSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY}
)

set(
	NORMSOURCES_UNI_SCAN_DETAIL_RES_PAIR
		uni/scan/detail/res_pair/multi_struc_res_rep_pair.cpp
		uni/scan/detail/res_pair/res_pair_core.cpp
		uni/scan/detail/res_pair/single_struc_res_pair.cpp
)

set(
	NORMSOURCES_UNI_SCAN_DETAIL_RES_PAIR_DIRN
		uni/scan/detail/res_pair_dirn/res_pair_dirn.cpp
)

set(
	NORMSOURCES_UNI_SCAN_DETAIL_STRIDE
		uni/scan/detail/stride/rep_strider.cpp
)

set(
	NORMSOURCES_UNI_SCAN_DETAIL
		${NORMSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN}
		${NORMSOURCES_UNI_SCAN_DETAIL_RES_PAIR}
		${NORMSOURCES_UNI_SCAN_DETAIL_RES_PAIR_DIRN}
		${NORMSOURCES_UNI_SCAN_DETAIL_STRIDE}
)

set(
	NORMSOURCES_UNI_SCAN_SCAN_ACTION
		uni/scan/scan_action/populate_matrix_scan_action.cpp
)

set(
	NORMSOURCES_UNI_SCAN_SCAN_TOOLS
		uni/scan/scan_tools/all_vs_all.cpp
		uni/scan/scan_tools/load_and_scan.cpp
		uni/scan/scan_tools/load_and_scan_metrics.cpp
		uni/scan/scan_tools/scan_metrics.cpp
		uni/scan/scan_tools/scan_type.cpp
		uni/scan/scan_tools/single_pair.cpp
)

set(
	NORMSOURCES_UNI_SCAN_SPATIAL_INDEX
		uni/scan/spatial_index/spatial_index.cpp
)

set(
	NORMSOURCES_UNI_SCAN
		${NORMSOURCES_UNI_SCAN_DETAIL}
		uni/scan/quad_criteria.cpp
		uni/scan/res_pair_index_dirn_criterion.cpp
		${NORMSOURCES_UNI_SCAN_SCAN_ACTION}
		uni/scan/scan_stride.cpp
		${NORMSOURCES_UNI_SCAN_SCAN_TOOLS}
		${NORMSOURCES_UNI_SCAN_SPATIAL_INDEX}
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_DETAIL
		uni/score/aligned_pair_score/detail/score_common_coord_handler.cpp
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE
		uni/score/aligned_pair_score/ssap_score/ssap_score_accuracy.cpp
		uni/score/aligned_pair_score/ssap_score/ssap_score_post_processing.cpp
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX
		uni/score/aligned_pair_score/substitution_matrix/blosum62_substitution_matrix.cpp
		uni/score/aligned_pair_score/substitution_matrix/identity_substitution_matrix.cpp
		uni/score/aligned_pair_score/substitution_matrix/match_substitution_matrix.cpp
		uni/score/aligned_pair_score/substitution_matrix/substitution_matrix.cpp
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE
		uni/score/aligned_pair_score/aligned_pair_score.cpp
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_DETAIL}
		uni/score/aligned_pair_score/drmsd_score.cpp
		uni/score/aligned_pair_score/gsas_score.cpp
		uni/score/aligned_pair_score/lddt_score.cpp
		uni/score/aligned_pair_score/length_score.cpp
		uni/score/aligned_pair_score/mi_score.cpp
		uni/score/aligned_pair_score/overlap_score.cpp
		uni/score/aligned_pair_score/pseudo_string_score.cpp
		uni/score/aligned_pair_score/rmsd_score.cpp
		uni/score/aligned_pair_score/sas_score.cpp
		uni/score/aligned_pair_score/sequence_similarity_score.cpp
		uni/score/aligned_pair_score/si_score.cpp
		uni/score/aligned_pair_score/ssap_score.cpp
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE}
		uni/score/aligned_pair_score/structal_score.cpp
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX}
		uni/score/aligned_pair_score/tm_score.cpp
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL
		uni/score/aligned_pair_score_list/detail/aligned_pair_score_list_append.cpp
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER
		uni/score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.cpp
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_READER
		uni/score/aligned_pair_score_list/score_value_list_reader/score_value_reader.cpp
)

set(
	NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST
		uni/score/aligned_pair_score_list/aligned_pair_score_list.cpp
		uni/score/aligned_pair_score_list/aligned_pair_score_list_factory.cpp
		uni/score/aligned_pair_score_list/aligned_pair_score_value_list.cpp
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL}
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER}
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_READER}
)

set(
	NORMSOURCES_UNI_SCORE_DETAIL
		uni/score/detail/score_name_helper.cpp
)

set(
	NORMSOURCES_UNI_SCORE_HOMCHECK_TOOLS
		uni/score/homcheck_tools/ssap_and_prc.cpp
		uni/score/homcheck_tools/ssaps_and_prcs_of_query.cpp
		uni/score/homcheck_tools/superfamily_of_domain.cpp
)

set(
	NORMSOURCES_UNI_SCORE_LENGTH_GETTER
		uni/score/length_getter/geometric_mean_length_getter.cpp
		uni/score/length_getter/length_getter.cpp
		uni/score/length_getter/length_getter_make_clone.cpp
		uni/score/length_getter/length_of_first_getter.cpp
		uni/score/length_getter/length_of_longer_getter.cpp
		uni/score/length_getter/length_of_second_getter.cpp
		uni/score/length_getter/length_of_shorter_getter.cpp
		uni/score/length_getter/mean_length_getter.cpp
		uni/score/length_getter/num_aligned_length_getter.cpp
		uni/score/length_getter/protein_only_length_getter.cpp
		uni/score/length_getter/sym_protein_only_length_getter.cpp
)

set(
	NORMSOURCES_UNI_SCORE_PAIR_SCATTER_PLOTTER
		uni/score/pair_scatter_plotter/pair_scatter_plotter.cpp
)

set(
	NORMSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_DETAIL
		uni/score/score_classification/detail/score_classn_value_list_name_less.cpp
)

set(
	NORMSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE
		uni/score/score_classification/label_pair_is_positive/label_pair_is_positive.cpp
)

set(
	NORMSOURCES_UNI_SCORE_SCORE_CLASSIFICATION
		${NORMSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_DETAIL}
		${NORMSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE}
		uni/score/score_classification/rbf_model.cpp
		uni/score/score_classification/score_classn_value.cpp
		uni/score/score_classification/score_classn_value_better_value.cpp
		uni/score/score_classification/score_classn_value_list.cpp
		uni/score/score_classification/score_classn_value_results_set.cpp
		uni/score/score_classification/value_list_scaling.cpp
)

set(
	NORMSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER
		uni/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter.cpp
		uni/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_spec.cpp
)

set(
	NORMSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG
		uni/score/true_pos_false_neg/classn_outcome.cpp
		uni/score/true_pos_false_neg/classn_rate_stat.cpp
		uni/score/true_pos_false_neg/classn_stat.cpp
		uni/score/true_pos_false_neg/classn_stat_pair_series.cpp
		uni/score/true_pos_false_neg/classn_stat_pair_series_list.cpp
		${NORMSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER}
		uni/score/true_pos_false_neg/named_true_false_pos_neg_list.cpp
		uni/score/true_pos_false_neg/named_true_false_pos_neg_list_list.cpp
		uni/score/true_pos_false_neg/true_false_pos_neg.cpp
		uni/score/true_pos_false_neg/true_false_pos_neg_list.cpp
)

set(
	NORMSOURCES_UNI_SCORE
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE}
		${NORMSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST}
		${NORMSOURCES_UNI_SCORE_DETAIL}
		${NORMSOURCES_UNI_SCORE_HOMCHECK_TOOLS}
		${NORMSOURCES_UNI_SCORE_LENGTH_GETTER}
		${NORMSOURCES_UNI_SCORE_PAIR_SCATTER_PLOTTER}
		${NORMSOURCES_UNI_SCORE_SCORE_CLASSIFICATION}
		${NORMSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG}
)

set(
	NORMSOURCES_UNI_SSAP_OPTIONS
		uni/ssap/options/cath_ssap_options.cpp
		uni/ssap/options/old_ssap_options_block.cpp
)

set(
	NORMSOURCES_UNI_SSAP
		uni/ssap/distance_score_formula.cpp
		${NORMSOURCES_UNI_SSAP_OPTIONS}
		uni/ssap/selected_pair.cpp
		uni/ssap/ssap.cpp
		uni/ssap/ssap_scores.cpp
		uni/ssap/windowed_matrix.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_ACCESSIBILITY_CALC
		uni/structure/accessibility_calc/dssp_accessibility.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_ENTRY_QUERIER
		uni/structure/entry_querier/entry_querier.cpp
		uni/structure/entry_querier/residue_querier.cpp
		uni/structure/entry_querier/sec_struc_querier.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_GEOMETRY_DETAIL
		uni/structure/geometry/detail/cross_covariance_matrix.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_GEOMETRY
		uni/structure/geometry/angle.cpp
		uni/structure/geometry/coord.cpp
		uni/structure/geometry/coord_list.cpp
		${NORMSOURCES_UNI_STRUCTURE_GEOMETRY_DETAIL}
		uni/structure/geometry/orient.cpp
		uni/structure/geometry/pca.cpp
		uni/structure/geometry/restrict_to_single_linkage_extension.cpp
		uni/structure/geometry/rotation.cpp
		uni/structure/geometry/superpose_fit.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_PROTEIN_PROTEIN_LOADER
		uni/structure/protein/protein_loader/protein_list_loader.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET
		uni/structure/protein/protein_source_file_set/protein_file_combn.cpp
		uni/structure/protein/protein_source_file_set/protein_from_pdb.cpp
		uni/structure/protein/protein_source_file_set/protein_from_pdb_and_calc.cpp
		uni/structure/protein/protein_source_file_set/protein_from_pdb_and_dssp_and_calc.cpp
		uni/structure/protein/protein_source_file_set/protein_from_pdb_dssp_and_sec.cpp
		uni/structure/protein/protein_source_file_set/protein_from_wolf_and_sec.cpp
		uni/structure/protein/protein_source_file_set/protein_source_file_set.cpp
		uni/structure/protein/protein_source_file_set/restrict_protein_source_file_set.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_PROTEIN
		uni/structure/protein/amino_acid.cpp
		uni/structure/protein/dna_atom.cpp
		uni/structure/protein/protein.cpp
		uni/structure/protein/protein_io.cpp
		uni/structure/protein/protein_list.cpp
		${NORMSOURCES_UNI_STRUCTURE_PROTEIN_PROTEIN_LOADER}
		${NORMSOURCES_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET}
		uni/structure/protein/residue.cpp
		uni/structure/protein/sec_struc.cpp
		uni/structure/protein/sec_struc_planar_angles.cpp
		uni/structure/protein/sec_struc_type.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP
		uni/structure/sec_struc_calc/dssp/bifur_hbond_list.cpp
		uni/structure/sec_struc_calc/dssp/dssp_hbond_calc.cpp
		uni/structure/sec_struc_calc/dssp/dssp_ss_calc.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_SEC
		uni/structure/sec_struc_calc/sec/sec_calc.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC
		${NORMSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP}
		${NORMSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_SEC}
)

set(
	NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER_DETAIL
		uni/structure/view_cache/filter/detail/filter_vs_full_score_less.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER
		${NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER_DETAIL}
		uni/structure/view_cache/filter/filter_vs_full_score.cpp
		uni/structure/view_cache/filter/filter_vs_full_score_list.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL
		uni/structure/view_cache/index/detail/vcie_match_criteria.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX
		${NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL}
		uni/structure/view_cache/index/quad_find_action.cpp
		uni/structure/view_cache/index/quad_find_action_check.cpp
		uni/structure/view_cache/index/view_cache_index.cpp
		uni/structure/view_cache/index/view_cache_index_entry.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE
		${NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER}
		${NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX}
		uni/structure/view_cache/view_cache.cpp
		uni/structure/view_cache/view_cache_list.cpp
)

set(
	NORMSOURCES_UNI_STRUCTURE
		${NORMSOURCES_UNI_STRUCTURE_ACCESSIBILITY_CALC}
		${NORMSOURCES_UNI_STRUCTURE_ENTRY_QUERIER}
		${NORMSOURCES_UNI_STRUCTURE_GEOMETRY}
		${NORMSOURCES_UNI_STRUCTURE_PROTEIN}
		${NORMSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC}
		${NORMSOURCES_UNI_STRUCTURE_VIEW_CACHE}
)

set(
	NORMSOURCES_UNI_SUPERPOSITION_IO
		uni/superposition/io/superposition_io.cpp
)

set(
	NORMSOURCES_UNI_SUPERPOSITION_OPTIONS
		uni/superposition/options/align_regions_options_block.cpp
		uni/superposition/options/superposition_content_options_block.cpp
)

set(
	NORMSOURCES_UNI_SUPERPOSITION
		${NORMSOURCES_UNI_SUPERPOSITION_IO}
		${NORMSOURCES_UNI_SUPERPOSITION_OPTIONS}
		uni/superposition/superpose_orient.cpp
		uni/superposition/superposition.cpp
		uni/superposition/superposition_content_spec.cpp
		uni/superposition/superposition_context.cpp
		uni/superposition/supn_regions_context.cpp
)

set(
	NORMSOURCES_UNI
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

set(
	TESTSOURCES_BIOCORE
		biocore/chain_label_test.cpp
		biocore/residue_id_test.cpp
		biocore/residue_name_test.cpp
)

set(
	TESTSOURCES_CATH_ASSIGN_DOMAINS_OPTIONS
		cath_assign_domains/options/cath_assign_domains_options_test.cpp
)

set(
	TESTSOURCES_CATH_ASSIGN_DOMAINS
		${TESTSOURCES_CATH_ASSIGN_DOMAINS_OPTIONS}
)

set(
	TESTSOURCES_CATH_REFINE_ALIGN_OPTIONS
		cath_refine_align/options/cath_refine_align_options_test.cpp
)

set(
	TESTSOURCES_CATH_REFINE_ALIGN
		cath_refine_align/cath_align_refiner_test.cpp
		${TESTSOURCES_CATH_REFINE_ALIGN_OPTIONS}
)

set(
	TESTSOURCES_CATH_SCORE_ALIGN_OPTIONS
		cath_score_align/options/cath_score_align_options_test.cpp
)

set(
	TESTSOURCES_CATH_SCORE_ALIGN
		cath_score_align/cath_align_scorer_test.cpp
		${TESTSOURCES_CATH_SCORE_ALIGN_OPTIONS}
)

set(
	TESTSOURCES_CATH_SUPERPOSE_OPTIONS
		cath_superpose/options/cath_superpose_options_test.cpp
)

set(
	TESTSOURCES_CATH_SUPERPOSE
		cath_superpose/cath_superposer_test.cpp
		${TESTSOURCES_CATH_SUPERPOSE_OPTIONS}
)

set(
	TESTSOURCES_CHOPPING_CHOPPING_FORMAT
		chopping/chopping_format/chopping_format_test.cpp
		chopping/chopping_format/sillitoe_chopping_format_test.cpp
		chopping/chopping_format/simple_chopping_format_test.cpp
)

set(
	TESTSOURCES_CHOPPING_DOMAIN
		chopping/domain/domain_test.cpp
)

set(
	TESTSOURCES_CHOPPING_REGION
		chopping/region/region_test.cpp
		chopping/region/regions_limiter_test.cpp
)

set(
	TESTSOURCES_CHOPPING
		${TESTSOURCES_CHOPPING_CHOPPING_FORMAT}
		chopping/chopping_test.cpp
		${TESTSOURCES_CHOPPING_DOMAIN}
		${TESTSOURCES_CHOPPING_REGION}
)

set(
	TESTSOURCES_CLUSTER_DETAIL
		cluster/detail/mapping_job_test.cpp
)

set(
	TESTSOURCES_CLUSTER_FILE
		cluster/file/cluster_membership_file_test.cpp
)

set(
	TESTSOURCES_CLUSTER_MAP
		cluster/map/aggregate_map_results_test.cpp
		cluster/map/map_results_test.cpp
		cluster/map/overlap_frac_distn_test.cpp
)

set(
	TESTSOURCES_CLUSTER_OPTIONS_OPTIONS_BLOCK
		cluster/options/options_block/clust_mapping_options_block_test.cpp
		cluster/options/options_block/clustmap_input_options_block_test.cpp
		cluster/options/options_block/clustmap_output_options_block_test.cpp
)

set(
	TESTSOURCES_CLUSTER_OPTIONS
		${TESTSOURCES_CLUSTER_OPTIONS_OPTIONS_BLOCK}
)

set(
	TESTSOURCES_CLUSTER_TEST
		cluster/test/map_clusters_fixture.cpp
)

set(
	TESTSOURCES_CLUSTER
		cluster/cath_cluster_mapper_test.cpp
		cluster/cluster_entry_test.cpp
		cluster/cluster_info_test.cpp
		cluster/clusters_info_test.cpp
		cluster/clustmap_options_test.cpp
		${TESTSOURCES_CLUSTER_DETAIL}
		cluster/domain_cluster_ids_by_seq_test.cpp
		cluster/domain_cluster_ids_test.cpp
		${TESTSOURCES_CLUSTER_FILE}
		${TESTSOURCES_CLUSTER_MAP}
		cluster/mapping_tool_test.cpp
		cluster/new_cluster_data_test.cpp
		cluster/old_cluster_data_test.cpp
		${TESTSOURCES_CLUSTER_OPTIONS}
		${TESTSOURCES_CLUSTER_TEST}
)

set(
	TESTSOURCES_DISPLAY_COLOUR
		display_colour/display_colour_gradient_test.cpp
		display_colour/display_colour_list_test.cpp
		display_colour/display_colour_test.cpp
)

set(
	TESTSOURCES_EXECUTABLES_BUILD_TEST
		executables/build_test/build_test.cpp
)

set(
	TESTSOURCES_EXECUTABLES
		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
)

set(
	TESTSOURCES_OPTIONS_EXECUTABLE
		options/executable/env_var_option_name_handler_test.cpp
		options/executable/executable_options_test.cpp
)

set(
	TESTSOURCES_OPTIONS_OPTIONS_BLOCK
		options/options_block/alignment_input_options_block_test.cpp
		options/options_block/check_pdb_options_block_test.cpp
		options/options_block/detail_help_options_block_test.cpp
		options/options_block/extract_pdb_options_block_test.cpp
		options/options_block/misc_help_version_options_block_test.cpp
		options/options_block/options_block_test.cpp
		options/options_block/pdb_input_options_block_test.cpp
		options/options_block/superposition_input_options_block_test.cpp
)

set(
	TESTSOURCES_OPTIONS
		${TESTSOURCES_OPTIONS_EXECUTABLE}
		${TESTSOURCES_OPTIONS_OPTIONS_BLOCK}
)

set(
	TESTSOURCES_RESOLVE_HITS_ALGO
		resolve_hits/algo/masked_bests_cache_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_FILE_DETAIL
		resolve_hits/file/detail/hmmer_parser_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_FILE
		resolve_hits/file/cath_id_score_category_test.cpp
		${TESTSOURCES_RESOLVE_HITS_FILE_DETAIL}
)

set(
	TESTSOURCES_RESOLVE_HITS_HTML_OUTPUT
		resolve_hits/html_output/resolve_hits_html_outputter_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK
		resolve_hits/options/options_block/crh_filter_options_block_test.cpp
		resolve_hits/options/options_block/crh_html_options_block_test.cpp
		resolve_hits/options/options_block/crh_input_options_block_test.cpp
		resolve_hits/options/options_block/crh_output_options_block_test.cpp
		resolve_hits/options/options_block/crh_score_options_block_test.cpp
		resolve_hits/options/options_block/crh_segment_options_block_test.cpp
		resolve_hits/options/options_block/crh_single_output_options_block_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_OPTIONS_SPEC
		resolve_hits/options/spec/crh_filter_spec_test.cpp
		resolve_hits/options/spec/crh_single_output_spec_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_OPTIONS
		resolve_hits/options/crh_options_test.cpp
		${TESTSOURCES_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK}
		${TESTSOURCES_RESOLVE_HITS_OPTIONS_SPEC}
)

set(
	TESTSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR
		resolve_hits/read_and_process_hits/hits_processor/hits_processor_list_test.cpp
		resolve_hits/read_and_process_hits/hits_processor/hits_processor_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS
		${TESTSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR}
)

set(
	TESTSOURCES_RESOLVE_HITS_RESOLVE
		resolve_hits/resolve/hit_resolver_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_TEST
		resolve_hits/test/resolve_hits_fixture.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_TRIM
		resolve_hits/trim/resolve_boundary_test.cpp
		resolve_hits/trim/seq_seg_boundary_fns_test.cpp
		resolve_hits/trim/trim_spec_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS
		${TESTSOURCES_RESOLVE_HITS_ALGO}
		resolve_hits/calc_hit_list_test.cpp
		resolve_hits/cath_hit_resolver_test.cpp
		${TESTSOURCES_RESOLVE_HITS_FILE}
		resolve_hits/first_hit_is_better_test.cpp
		resolve_hits/full_hit_list_test.cpp
		resolve_hits/full_hit_test.cpp
		resolve_hits/hit_extras_test.cpp
		resolve_hits/hit_test.cpp
		${TESTSOURCES_RESOLVE_HITS_HTML_OUTPUT}
		${TESTSOURCES_RESOLVE_HITS_OPTIONS}
		${TESTSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS}
		${TESTSOURCES_RESOLVE_HITS_RESOLVE}
		${TESTSOURCES_RESOLVE_HITS_TEST}
		${TESTSOURCES_RESOLVE_HITS_TRIM}
)

set(
	TESTSOURCES_SEQ
		seq/seq_seg_run_test.cpp
		seq/seq_seg_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_ALGORITHM
		src_common/common/algorithm/are_same_test.cpp
		src_common/common/algorithm/constexpr_find_test.cpp
		src_common/common/algorithm/constexpr_floor_test.cpp
		src_common/common/algorithm/constexpr_integer_rounding_test.cpp
		src_common/common/algorithm/constexpr_is_uniq_test.cpp
		src_common/common/algorithm/constexpr_modulo_fns_test.cpp
		src_common/common/algorithm/for_n_test.cpp
		src_common/common/algorithm/transform_build_test.cpp
		src_common/common/algorithm/variadic_and_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_GRAPH
		src_common/common/boost_addenda/graph/spanning_tree_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR
		src_common/common/boost_addenda/range/adaptor/adaptor_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR
		src_common/common/boost_addenda/range/utility/iterator/cross_itr_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_UTILITY
		${TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR}
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE
		${TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR}
		${TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_UTILITY}
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_TRIBOOL
		src_common/common/boost_addenda/tribool/tribool_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA
		${TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_GRAPH}
		${TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE}
		${TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA_TRIBOOL}
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_CLONE
		src_common/common/clone/clone_ptr_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_CONTAINER
		src_common/common/container/id_of_str_bidirnl_test.cpp
		src_common/common/container/id_of_string_ref_test.cpp
		src_common/common/container/id_of_string_test.cpp
		src_common/common/container/id_of_string_view_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_CPP17
		src_common/common/cpp17/apply_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_EXCEPTION
		src_common/common/exception/exception_is_equivalent_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_FILE
		src_common/common/file/ofstream_list_test.cpp
		src_common/common/file/open_fstream_test.cpp
		src_common/common/file/simple_file_read_write_test.cpp
		src_common/common/file/temp_file_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_FUNCTION
		src_common/common/function/multiply_args_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_GSL
		src_common/common/gsl/get_determinant_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_METAPROGRAMMING
		src_common/common/metaprogramming/append_template_params_into_first_wrapper_test.cpp
		src_common/common/metaprogramming/change_template_subwrappers_test.cpp
		src_common/common/metaprogramming/change_template_wrapper_test.cpp
		src_common/common/metaprogramming/combine_params_lists_with_template_list_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_OPTIONAL
		src_common/common/optional/make_optional_if_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_RAPIDJSON_ADDENDA
		src_common/common/rapidjson_addenda/rapidjson_writer_test.cpp
		src_common/common/rapidjson_addenda/string_of_rapidjson_write_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_STRING
		src_common/common/string/booled_to_string_test.cpp
		src_common/common/string/string_parse_tools_test.cpp
		src_common/common/string/sub_string_parser_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_TUPLE
		src_common/common/tuple/make_tuple_with_skips_test.cpp
		src_common/common/tuple/mins_maxs_tuple_pair_mins_maxs_element_test.cpp
		src_common/common/tuple/tuple_increment_test.cpp
		src_common/common/tuple/tuple_lattice_index_test.cpp
		src_common/common/tuple/tuple_mins_maxs_element_test.cpp
		src_common/common/tuple/tuple_multiply_args_test.cpp
		src_common/common/tuple/tuple_subtract_test.cpp
		src_common/common/tuple/tuple_within_range_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON_TYPE_TRAITS
		src_common/common/type_traits/is_tuple_test.cpp
)

set(
	TESTSOURCES_SRC_COMMON_COMMON
		${TESTSOURCES_SRC_COMMON_COMMON_ALGORITHM}
		src_common/common/argc_argv_faker_test.cpp
		${TESTSOURCES_SRC_COMMON_COMMON_BOOST_ADDENDA}
		${TESTSOURCES_SRC_COMMON_COMMON_CLONE}
		src_common/common/command_executer_test.cpp
		${TESTSOURCES_SRC_COMMON_COMMON_CONTAINER}
		${TESTSOURCES_SRC_COMMON_COMMON_CPP17}
		src_common/common/difference_test.cpp
		${TESTSOURCES_SRC_COMMON_COMMON_EXCEPTION}
		${TESTSOURCES_SRC_COMMON_COMMON_FILE}
		${TESTSOURCES_SRC_COMMON_COMMON_FUNCTION}
		${TESTSOURCES_SRC_COMMON_COMMON_GSL}
		src_common/common/invert_permutation_test.cpp
		src_common/common/less_than_helper_test.cpp
		src_common/common/logger_test.cpp
		${TESTSOURCES_SRC_COMMON_COMMON_METAPROGRAMMING}
		${TESTSOURCES_SRC_COMMON_COMMON_OPTIONAL}
		src_common/common/program_exception_wrapper_test.cpp
		${TESTSOURCES_SRC_COMMON_COMMON_RAPIDJSON_ADDENDA}
		${TESTSOURCES_SRC_COMMON_COMMON_STRING}
		src_common/common/temp_check_offset_1_test.cpp
		${TESTSOURCES_SRC_COMMON_COMMON_TUPLE}
		src_common/common/tuple_insertion_operator_test.cpp
		src_common/common/type_to_string_test.cpp
		${TESTSOURCES_SRC_COMMON_COMMON_TYPE_TRAITS}
)

set(
	TESTSOURCES_SRC_COMMON
		${TESTSOURCES_SRC_COMMON_COMMON}
)

set(
	TESTSOURCES_SRC_TEST_TEST_PREDICATE
		src_test/test/predicate/files_equal_test.cpp
		src_test/test/predicate/istreams_equal_test.cpp
)

set(
	TESTSOURCES_SRC_TEST_TEST
		${TESTSOURCES_SRC_TEST_TEST_PREDICATE}
		src_test/test/superposition_fixture.cpp
)

set(
	TESTSOURCES_SRC_TEST
		${TESTSOURCES_SRC_TEST_TEST}
)

set(
	TESTSOURCES_UNI_ACQUIRER_ALIGNMENT_ACQUIRER
		uni/acquirer/alignment_acquirer/alignment_acquirer_test.cpp
		uni/acquirer/alignment_acquirer/cora_aln_file_alignment_acquirer_test.cpp
		uni/acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer_test.cpp
		uni/acquirer/alignment_acquirer/residue_name_alignment_acquirer_test.cpp
		uni/acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer_test.cpp
		uni/acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer_test.cpp
)

set(
	TESTSOURCES_UNI_ACQUIRER_PDBS_ACQUIRER
		uni/acquirer/pdbs_acquirer/file_list_pdbs_acquirer_test.cpp
		uni/acquirer/pdbs_acquirer/istream_pdbs_acquirer_test.cpp
		uni/acquirer/pdbs_acquirer/pdbs_acquirer_test.cpp
)

set(
	TESTSOURCES_UNI_ACQUIRER_SELECTION_POLICY_ACQUIRER
		uni/acquirer/selection_policy_acquirer/selection_policy_acquirer_test.cpp
)

set(
	TESTSOURCES_UNI_ACQUIRER_SUPERPOSITION_ACQUIRER
		uni/acquirer/superposition_acquirer/align_based_superposition_acquirer_test.cpp
		uni/acquirer/superposition_acquirer/superposition_acquirer_test.cpp
)

set(
	TESTSOURCES_UNI_ACQUIRER
		${TESTSOURCES_UNI_ACQUIRER_ALIGNMENT_ACQUIRER}
		${TESTSOURCES_UNI_ACQUIRER_PDBS_ACQUIRER}
		${TESTSOURCES_UNI_ACQUIRER_SELECTION_POLICY_ACQUIRER}
		${TESTSOURCES_UNI_ACQUIRER_SUPERPOSITION_ACQUIRER}
)

set(
	TESTSOURCES_UNI_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY
		uni/alignment/common_residue_selection_policy/common_residue_selection_policy_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER
		uni/alignment/dyn_prog_align/detail/string_aligner/string_aligner_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL
		uni/alignment/dyn_prog_align/detail/return_path_matrix_test.cpp
		uni/alignment/dyn_prog_align/detail/score_accumulation_matrix_test.cpp
		${TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER}
)

set(
	TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE
		uni/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source_test.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/entry_querier_dyn_prog_score_source_test.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source_test.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source_test.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source_test.cpp
		uni/alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_TEST
		uni/alignment/dyn_prog_align/test/dyn_prog_score_source_fixture.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN
		${TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL}
		${TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE}
		uni/alignment/dyn_prog_align/ssap_code_dyn_prog_aligner_test.cpp
		uni/alignment/dyn_prog_align/std_dyn_prog_aligner_test.cpp
		${TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN_TEST}
)

set(
	TESTSOURCES_UNI_ALIGNMENT_GAP
		uni/alignment/gap/alignment_gap_test.cpp
		uni/alignment/gap/gap_penalty_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_IO_OUTPUTTER
		uni/alignment/io/outputter/horiz_align_outputter_test.cpp
		uni/alignment/io/outputter/html_align_outputter_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_IO
		uni/alignment/io/align_scaffold_test.cpp
		uni/alignment/io/alignment_io_test.cpp
		${TESTSOURCES_UNI_ALIGNMENT_IO_OUTPUTTER}
)

set(
	TESTSOURCES_UNI_ALIGNMENT_REFINER_DETAIL
		uni/alignment/refiner/detail/alignment_split_list_test.cpp
		uni/alignment/refiner/detail/alignment_split_mapping_test.cpp
		uni/alignment/refiner/detail/alignment_split_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_REFINER
		uni/alignment/refiner/alignment_refiner_test.cpp
		${TESTSOURCES_UNI_ALIGNMENT_REFINER_DETAIL}
		uni/alignment/refiner/indexed_refiner_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL
		uni/alignment/residue_name_align/detail/residue_name_align_map_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN
		${TESTSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL}
		uni/alignment/residue_name_align/residue_name_aligner_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_RESIDUE_SCORE
		uni/alignment/residue_score/alignment_residue_scores_test.cpp
		uni/alignment/residue_score/residue_scorer_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_TEST
		uni/alignment/test/alignment_fixture.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT_TOOLS
		uni/alignment/tools/alignment_breaks_test.cpp
)

set(
	TESTSOURCES_UNI_ALIGNMENT
		uni/alignment/alignment_action_test.cpp
		uni/alignment/alignment_context_test.cpp
		uni/alignment/alignment_coord_extractor_test.cpp
		uni/alignment/alignment_row_test.cpp
		uni/alignment/alignment_test.cpp
		${TESTSOURCES_UNI_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY}
		${TESTSOURCES_UNI_ALIGNMENT_DYN_PROG_ALIGN}
		${TESTSOURCES_UNI_ALIGNMENT_GAP}
		${TESTSOURCES_UNI_ALIGNMENT_IO}
		${TESTSOURCES_UNI_ALIGNMENT_REFINER}
		${TESTSOURCES_UNI_ALIGNMENT_RESIDUE_NAME_ALIGN}
		${TESTSOURCES_UNI_ALIGNMENT_RESIDUE_SCORE}
		${TESTSOURCES_UNI_ALIGNMENT_TEST}
		${TESTSOURCES_UNI_ALIGNMENT_TOOLS}
)

set(
	TESTSOURCES_UNI_DISPLAY_DISPLAY_COLOUR_SPEC
		uni/display/display_colour_spec/broad_display_colour_spec_test.cpp
		uni/display/display_colour_spec/display_colour_spec_test.cpp
)

set(
	TESTSOURCES_UNI_DISPLAY_DISPLAY_COLOURER
		uni/display/display_colourer/display_colourer_test.cpp
)

set(
	TESTSOURCES_UNI_DISPLAY_OPTIONS
		uni/display/options/display_options_block_test.cpp
		uni/display/options/display_spec_test.cpp
)

set(
	TESTSOURCES_UNI_DISPLAY_VIEWER_PYMOL
		uni/display/viewer/pymol/pymol_tools_test.cpp
)

set(
	TESTSOURCES_UNI_DISPLAY_VIEWER
		uni/display/viewer/jmol_viewer_test.cpp
		${TESTSOURCES_UNI_DISPLAY_VIEWER_PYMOL}
		uni/display/viewer/pymol_viewer_test.cpp
		uni/display/viewer/rasmol_style_viewer_test.cpp
		uni/display/viewer/rasmol_viewer_test.cpp
		uni/display/viewer/viewer_test.cpp
)

set(
	TESTSOURCES_UNI_DISPLAY
		${TESTSOURCES_UNI_DISPLAY_DISPLAY_COLOUR_SPEC}
		${TESTSOURCES_UNI_DISPLAY_DISPLAY_COLOURER}
		${TESTSOURCES_UNI_DISPLAY_OPTIONS}
		${TESTSOURCES_UNI_DISPLAY_VIEWER}
)

set(
	TESTSOURCES_UNI_FILE_DSSP_WOLF
		uni/file/dssp_wolf/dssp_file_test.cpp
		uni/file/dssp_wolf/tally_residue_ids_test.cpp
)

set(
	TESTSOURCES_UNI_FILE_HMMER_SCORES_FILE
		uni/file/hmmer_scores_file/hmmer_scores_entry_test.cpp
		uni/file/hmmer_scores_file/hmmer_scores_file_test.cpp
)

set(
	TESTSOURCES_UNI_FILE_NAME_SET
		uni/file/name_set/name_set_list_test.cpp
		uni/file/name_set/name_set_test.cpp
)

set(
	TESTSOURCES_UNI_FILE_OPTIONS
		uni/file/options/data_dirs_options_block_test.cpp
		uni/file/options/data_dirs_spec_test.cpp
)

set(
	TESTSOURCES_UNI_FILE_PDB
		uni/file/pdb/coarse_element_type_test.cpp
		uni/file/pdb/element_type_string_test.cpp
		uni/file/pdb/pdb_atom_test.cpp
		uni/file/pdb/pdb_list_test.cpp
		uni/file/pdb/pdb_residue_test.cpp
		uni/file/pdb/pdb_test.cpp
		uni/file/pdb/proximity_calculator_test.cpp
)

set(
	TESTSOURCES_UNI_FILE_PRC_SCORES_FILE
		uni/file/prc_scores_file/prc_scores_entry_test.cpp
		uni/file/prc_scores_file/prc_scores_file_test.cpp
)

set(
	TESTSOURCES_UNI_FILE_SEC
		uni/file/sec/sec_file_test.cpp
)

set(
	TESTSOURCES_UNI_FILE_SSAP_SCORES_FILE
		uni/file/ssap_scores_file/ssap_scores_entry_test.cpp
		uni/file/ssap_scores_file/ssap_scores_file_test.cpp
)

set(
	TESTSOURCES_UNI_FILE
		${TESTSOURCES_UNI_FILE_DSSP_WOLF}
		${TESTSOURCES_UNI_FILE_HMMER_SCORES_FILE}
		${TESTSOURCES_UNI_FILE_NAME_SET}
		${TESTSOURCES_UNI_FILE_OPTIONS}
		${TESTSOURCES_UNI_FILE_PDB}
		${TESTSOURCES_UNI_FILE_PRC_SCORES_FILE}
		${TESTSOURCES_UNI_FILE_SEC}
		${TESTSOURCES_UNI_FILE_SSAP_SCORES_FILE}
)

set(
	TESTSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER_TEST
		uni/outputter/alignment_outputter/test/alignment_outputter_fixture.cpp
)

set(
	TESTSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER
		uni/outputter/alignment_outputter/alignment_outputter_list_test.cpp
		uni/outputter/alignment_outputter/ssap_ostream_alignment_outputter_test.cpp
		${TESTSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER_TEST}
)

set(
	TESTSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS
		uni/outputter/alignment_outputter_options/alignment_output_options_block_test.cpp
)

set(
	TESTSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS
		uni/outputter/superposition_output_options/superposition_output_options_block_test.cpp
)

set(
	TESTSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUTTER
		uni/outputter/superposition_outputter/json_file_superposition_outputter_test.cpp
		uni/outputter/superposition_outputter/ostream_superposition_outputter_test.cpp
		uni/outputter/superposition_outputter/pdb_file_superposition_outputter_test.cpp
		uni/outputter/superposition_outputter/pdb_files_superposition_outputter_test.cpp
		uni/outputter/superposition_outputter/pymol_file_superposition_outputter_test.cpp
		uni/outputter/superposition_outputter/pymol_view_superposition_outputter_test.cpp
		uni/outputter/superposition_outputter/superposition_outputter_list_test.cpp
		uni/outputter/superposition_outputter/superposition_outputter_test.cpp
)

set(
	TESTSOURCES_UNI_OUTPUTTER
		${TESTSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER}
		${TESTSOURCES_UNI_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS}
		${TESTSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS}
		${TESTSOURCES_UNI_OUTPUTTER_SUPERPOSITION_OUTPUTTER}
)

set(
	TESTSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY
		uni/scan/detail/check_scan/test_only/check_scan_on_final_alignment_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN
		${TESTSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY}
)

set(
	TESTSOURCES_UNI_SCAN_DETAIL_RES_PAIR
		uni/scan/detail/res_pair/multi_struc_res_rep_pair_list_test.cpp
		uni/scan/detail/res_pair/multi_struc_res_rep_pair_test.cpp
		uni/scan/detail/res_pair/res_pair_core_test.cpp
		uni/scan/detail/res_pair/single_struc_res_pair_list_test.cpp
		uni/scan/detail/res_pair/single_struc_res_pair_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN_DETAIL_SCAN_ACTION
		uni/scan/detail/scan_action/align_scan_action_test.cpp
		uni/scan/detail/scan_action/scan_multi_action_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN_DETAIL_SCAN_INDEX_STORE
		uni/scan/detail/scan_index_store/scan_index_hash_store_test.cpp
		uni/scan/detail/scan_index_store/scan_index_lattice_store_test.cpp
		uni/scan/detail/scan_index_store/scan_index_vector_store_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN_DETAIL_STRIDE
		uni/scan/detail/stride/co_stride_test.cpp
		uni/scan/detail/stride/rep_strider_test.cpp
		uni/scan/detail/stride/roled_scan_stride_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN_DETAIL
		${TESTSOURCES_UNI_SCAN_DETAIL_CHECK_SCAN}
		${TESTSOURCES_UNI_SCAN_DETAIL_RES_PAIR}
		${TESTSOURCES_UNI_SCAN_DETAIL_SCAN_ACTION}
		${TESTSOURCES_UNI_SCAN_DETAIL_SCAN_INDEX_STORE}
		uni/scan/detail/scan_structure_data_test.cpp
		${TESTSOURCES_UNI_SCAN_DETAIL_STRIDE}
)

set(
	TESTSOURCES_UNI_SCAN_RES_PAIR_KEYER
		uni/scan/res_pair_keyer/res_pair_keyer_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN_SCAN_TOOLS
		uni/scan/scan_tools/all_vs_all_test.cpp
		uni/scan/scan_tools/load_and_scan_metrics_test.cpp
		uni/scan/scan_tools/load_and_scan_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN_SPATIAL_INDEX
		uni/scan/spatial_index/spatial_index_test.cpp
)

set(
	TESTSOURCES_UNI_SCAN
		${TESTSOURCES_UNI_SCAN_DETAIL}
		uni/scan/quad_criteria_test.cpp
		${TESTSOURCES_UNI_SCAN_RES_PAIR_KEYER}
		uni/scan/scan_index_test.cpp
		uni/scan/scan_policy_test.cpp
		uni/scan/scan_query_set_test.cpp
		uni/scan/scan_stride_test.cpp
		${TESTSOURCES_UNI_SCAN_SCAN_TOOLS}
		${TESTSOURCES_UNI_SCAN_SPATIAL_INDEX}
)

set(
	TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_DETAIL
		uni/score/aligned_pair_score/detail/score_common_coord_handler_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE
		uni/score/aligned_pair_score/ssap_score/ssap_score_accuracy_test.cpp
		uni/score/aligned_pair_score/ssap_score/ssap_score_post_processing_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX
		uni/score/aligned_pair_score/substitution_matrix/substitution_matrix_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE
		uni/score/aligned_pair_score/aligned_pair_score_test.cpp
		${TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_DETAIL}
		uni/score/aligned_pair_score/rmsd_score_test.cpp
		${TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE}
		${TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX}
)

set(
	TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST
		uni/score/aligned_pair_score_list/aligned_pair_score_list_factory_test.cpp
		uni/score/aligned_pair_score_list/aligned_pair_score_list_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_HOMCHECK_TOOLS
		uni/score/homcheck_tools/ssap_and_prc_test.cpp
		uni/score/homcheck_tools/ssaps_and_prcs_of_query_test.cpp
		uni/score/homcheck_tools/superfamily_of_domain_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_LENGTH_GETTER
		uni/score/length_getter/length_getter_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_DETAIL
		uni/score/score_classification/detail/score_classn_value_list_name_less_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE
		uni/score/score_classification/label_pair_is_positive/label_pair_is_positive_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_SCORE_CLASSIFICATION
		${TESTSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_DETAIL}
		${TESTSOURCES_UNI_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE}
		uni/score/score_classification/rbf_model_test.cpp
		uni/score/score_classification/score_classn_value_better_value_test.cpp
		uni/score/score_classification/score_classn_value_list_test.cpp
		uni/score/score_classification/score_classn_value_results_set_test.cpp
		uni/score/score_classification/score_classn_value_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER
		uni/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG
		uni/score/true_pos_false_neg/classn_rate_stat_test.cpp
		uni/score/true_pos_false_neg/classn_stat_pair_series_list_test.cpp
		uni/score/true_pos_false_neg/classn_stat_pair_series_test.cpp
		${TESTSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER}
		uni/score/true_pos_false_neg/classn_stat_test.cpp
		uni/score/true_pos_false_neg/named_true_false_pos_neg_list_test.cpp
		uni/score/true_pos_false_neg/true_false_pos_neg_list_test.cpp
		uni/score/true_pos_false_neg/true_false_pos_neg_test.cpp
)

set(
	TESTSOURCES_UNI_SCORE
		${TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE}
		${TESTSOURCES_UNI_SCORE_ALIGNED_PAIR_SCORE_LIST}
		${TESTSOURCES_UNI_SCORE_HOMCHECK_TOOLS}
		${TESTSOURCES_UNI_SCORE_LENGTH_GETTER}
		${TESTSOURCES_UNI_SCORE_SCORE_CLASSIFICATION}
		${TESTSOURCES_UNI_SCORE_TRUE_POS_FALSE_NEG}
)

set(
	TESTSOURCES_UNI_SSAP_OPTIONS
		uni/ssap/options/cath_ssap_options_test.cpp
		uni/ssap/options/old_ssap_options_block_test.cpp
)

set(
	TESTSOURCES_UNI_SSAP
		uni/ssap/distance_score_formula_test.cpp
		${TESTSOURCES_UNI_SSAP_OPTIONS}
		uni/ssap/selected_pair_test.cpp
		uni/ssap/ssap_scores_test.cpp
		uni/ssap/ssap_test.cpp
		uni/ssap/windowed_matrix_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_ACCESSIBILITY_CALC
		uni/structure/accessibility_calc/dssp_accessibility_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_ENTRY_QUERIER
		uni/structure/entry_querier/entry_querier_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_GEOMETRY
		uni/structure/geometry/angle_test.cpp
		uni/structure/geometry/coord_list_test.cpp
		uni/structure/geometry/coord_test.cpp
		uni/structure/geometry/orient_test.cpp
		uni/structure/geometry/orientation_covering_test.cpp
		uni/structure/geometry/pca_test.cpp
		uni/structure/geometry/quat_rot_test.cpp
		uni/structure/geometry/restrict_to_single_linkage_extension_test.cpp
		uni/structure/geometry/rotation_test.cpp
		uni/structure/geometry/superpose_fit_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET
		uni/structure/protein/protein_source_file_set/protein_source_file_set_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_PROTEIN
		uni/structure/protein/amino_acid_test.cpp
		uni/structure/protein/protein_list_test.cpp
		${TESTSOURCES_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET}
		uni/structure/protein/protein_test.cpp
		uni/structure/protein/residue_test.cpp
		uni/structure/protein/sec_struc_planar_angles_test.cpp
		uni/structure/protein/sec_struc_test.cpp
		uni/structure/protein/sec_struc_type_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST
		uni/structure/sec_struc_calc/dssp/test/dssp_dupl_fixture.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP
		uni/structure/sec_struc_calc/dssp/bifur_hbond_list_test.cpp
		uni/structure/sec_struc_calc/dssp/dssp_hbond_calc_test.cpp
		uni/structure/sec_struc_calc/dssp/dssp_ss_calc_test.cpp
		${TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST}
)

set(
	TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_SEC
		uni/structure/sec_struc_calc/sec/sec_calc_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC
		${TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP}
		${TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC_SEC}
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_DETAIL_PLATE
		uni/structure/view_cache/detail/plate/plate_scan_test.cpp
		uni/structure/view_cache/detail/plate/rod_cache_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_DETAIL
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_DETAIL_PLATE}
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER_DETAIL
		uni/structure/view_cache/filter/detail/filter_vs_full_score_less_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER_DETAIL}
		uni/structure/view_cache/filter/filter_vs_full_score_list_test.cpp
		uni/structure/view_cache/filter/filter_vs_full_score_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL
		uni/structure/view_cache/index/detail/dims/detail/view_cache_index_dim_linear_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL}
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS}
		uni/structure/view_cache/index/detail/vcie_match_criteria_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL}
		uni/structure/view_cache/index/view_cache_index_entry_test.cpp
		uni/structure/view_cache/index/view_cache_index_test.cpp
)

set(
	TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_DETAIL}
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_FILTER}
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE_INDEX}
)

set(
	TESTSOURCES_UNI_STRUCTURE
		${TESTSOURCES_UNI_STRUCTURE_ACCESSIBILITY_CALC}
		${TESTSOURCES_UNI_STRUCTURE_ENTRY_QUERIER}
		${TESTSOURCES_UNI_STRUCTURE_GEOMETRY}
		uni/structure/get_residue_names_test.cpp
		${TESTSOURCES_UNI_STRUCTURE_PROTEIN}
		${TESTSOURCES_UNI_STRUCTURE_SEC_STRUC_CALC}
		${TESTSOURCES_UNI_STRUCTURE_VIEW_CACHE}
)

set(
	TESTSOURCES_UNI_SUPERPOSITION_IO
		uni/superposition/io/superposition_io_test.cpp
)

set(
	TESTSOURCES_UNI_SUPERPOSITION_OPTIONS
		uni/superposition/options/align_regions_options_block_test.cpp
		uni/superposition/options/superposition_content_options_block_test.cpp
)

set(
	TESTSOURCES_UNI_SUPERPOSITION
		${TESTSOURCES_UNI_SUPERPOSITION_IO}
		${TESTSOURCES_UNI_SUPERPOSITION_OPTIONS}
		uni/superposition/superposition_context_test.cpp
		uni/superposition/superposition_test.cpp
)

set(
	TESTSOURCES_UNI
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
##### DON'T EDIT THIS FILE - IT'S AUTO-GENERATED #####
