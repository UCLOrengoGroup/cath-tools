include_directories( ${Boost_INCLUDE_DIR} )

##### BEGINNING OF AUTO-GENERATED #####
set(
	NORMSOURCES_ACQUIRER_ALIGNMENT_ACQUIRER
		acquirer/alignment_acquirer/alignment_acquirer.cpp
		acquirer/alignment_acquirer/cora_aln_file_alignment_acquirer.cpp
		acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer.cpp
		acquirer/alignment_acquirer/residue_name_alignment_acquirer.cpp
		acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer.cpp
		acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.cpp
)

set(
	NORMSOURCES_ACQUIRER_PDBS_ACQUIRER
		acquirer/pdbs_acquirer/domain_defn_pdbs_acquirer.cpp
		acquirer/pdbs_acquirer/file_list_pdbs_acquirer.cpp
		acquirer/pdbs_acquirer/istream_pdbs_acquirer.cpp
		acquirer/pdbs_acquirer/pdbs_acquirer.cpp
)

set(
	NORMSOURCES_ACQUIRER_SELECTION_POLICY_ACQUIRER
		acquirer/selection_policy_acquirer/selection_policy_acquirer.cpp
)

set(
	NORMSOURCES_ACQUIRER_SUPERPOSITION_ACQUIRER
		acquirer/superposition_acquirer/align_based_superposition_acquirer.cpp
		acquirer/superposition_acquirer/superposition_acquirer.cpp
)

set(
	NORMSOURCES_ACQUIRER
		${NORMSOURCES_ACQUIRER_ALIGNMENT_ACQUIRER}
		${NORMSOURCES_ACQUIRER_PDBS_ACQUIRER}
		${NORMSOURCES_ACQUIRER_SELECTION_POLICY_ACQUIRER}
		${NORMSOURCES_ACQUIRER_SUPERPOSITION_ACQUIRER}
)

set(
	NORMSOURCES_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY
		alignment/common_atom_selection_policy/common_atom_select_ca_policy.cpp
		alignment/common_atom_selection_policy/common_atom_select_cb_policy.cpp
		alignment/common_atom_selection_policy/common_atom_selection_policy.cpp
)

set(
	NORMSOURCES_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY
		alignment/common_residue_selection_policy/common_residue_score_based_selection_policy.cpp
		alignment/common_residue_selection_policy/common_residue_select_all_policy.cpp
		alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.cpp
		alignment/common_residue_selection_policy/common_residue_select_min_score_policy.cpp
		alignment/common_residue_selection_policy/common_residue_selection_policy.cpp
)

set(
	NORMSOURCES_ALIGNMENT_DETAIL
		alignment/detail/multi_align_builder.cpp
		alignment/detail/multi_align_group.cpp
)

set(
	NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER
		alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.cpp
		alignment/dyn_prog_align/detail/matrix_plotter/matrix_plotter.cpp
)

set(
	NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER
		alignment/dyn_prog_align/detail/string_aligner/benchmark_dyn_prog_string_aligner.cpp
		alignment/dyn_prog_align/detail/string_aligner/gen_dyn_prog_string_aligner.cpp
		alignment/dyn_prog_align/detail/string_aligner/string_aligner.cpp
)

set(
	NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL
		${NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER}
		alignment/dyn_prog_align/detail/path_step.cpp
		alignment/dyn_prog_align/detail/return_path_matrix.cpp
		alignment/dyn_prog_align/detail/score_accumulation_matrix.cpp
		${NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER}
)

set(
	NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE
		alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/entry_querier_dyn_prog_score_source.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.cpp
)

set(
	NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN
		${NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL}
		alignment/dyn_prog_align/dyn_prog_aligner.cpp
		${NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE}
		alignment/dyn_prog_align/ssap_code_dyn_prog_aligner.cpp
		alignment/dyn_prog_align/std_dyn_prog_aligner.cpp
)

set(
	NORMSOURCES_ALIGNMENT_GAP
		alignment/gap/alignment_gap.cpp
		alignment/gap/gap_penalty.cpp
)

set(
	NORMSOURCES_ALIGNMENT_IO_OUTPUTTER
		alignment/io/outputter/horiz_align_outputter.cpp
		alignment/io/outputter/html_align_outputter.cpp
)

set(
	NORMSOURCES_ALIGNMENT_IO
		alignment/io/align_scaffold.cpp
		alignment/io/alignment_io.cpp
		${NORMSOURCES_ALIGNMENT_IO_OUTPUTTER}
)

set(
	NORMSOURCES_ALIGNMENT_REFINER_DETAIL
		alignment/refiner/detail/alignment_split.cpp
		alignment/refiner/detail/alignment_split_list.cpp
		alignment/refiner/detail/alignment_split_mapping.cpp
)

set(
	NORMSOURCES_ALIGNMENT_REFINER
		alignment/refiner/alignment_refiner.cpp
		${NORMSOURCES_ALIGNMENT_REFINER_DETAIL}
		alignment/refiner/indexed_refiner.cpp
)

set(
	NORMSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL
		alignment/residue_name_align/detail/residue_name_align_map.cpp
)

set(
	NORMSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN
		${NORMSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL}
		alignment/residue_name_align/residue_name_aligner.cpp
)

set(
	NORMSOURCES_ALIGNMENT_RESIDUE_SCORE
		alignment/residue_score/alignment_residue_scores.cpp
		alignment/residue_score/residue_scorer.cpp
)

set(
	NORMSOURCES_ALIGNMENT_TOOLS
		alignment/tools/alignment_breaks.cpp
)

set(
	NORMSOURCES_ALIGNMENT
		alignment/alignment.cpp
		alignment/alignment_action.cpp
		alignment/alignment_context.cpp
		alignment/alignment_coord_extractor.cpp
		alignment/alignment_row.cpp
		${NORMSOURCES_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY}
		${NORMSOURCES_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY}
		${NORMSOURCES_ALIGNMENT_DETAIL}
		${NORMSOURCES_ALIGNMENT_DYN_PROG_ALIGN}
		${NORMSOURCES_ALIGNMENT_GAP}
		${NORMSOURCES_ALIGNMENT_IO}
		alignment/pair_alignment.cpp
		${NORMSOURCES_ALIGNMENT_REFINER}
		${NORMSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN}
		${NORMSOURCES_ALIGNMENT_RESIDUE_SCORE}
		${NORMSOURCES_ALIGNMENT_TOOLS}
)

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
	NORMSOURCES_COMMON_ALGORITHM
		common/algorithm/random_split.cpp
)

set(
	NORMSOURCES_COMMON_BATCH
		common/batch/batch_functions.cpp
)

set(
	NORMSOURCES_COMMON_BOOST_ADDENDA_LOG
		common/boost_addenda/log/log_to_ostream_guard.cpp
)

set(
	NORMSOURCES_COMMON_BOOST_ADDENDA
		${NORMSOURCES_COMMON_BOOST_ADDENDA_LOG}
)

set(
	NORMSOURCES_COMMON_FILE
		common/file/find_file.cpp
		common/file/open_fstream.cpp
		common/file/temp_file.cpp
)

set(
	NORMSOURCES_COMMON_TEST_PREDICATE
		common/test_predicate/files_equal.cpp
		common/test_predicate/istream_and_file_equal.cpp
		common/test_predicate/istreams_equal.cpp
)

set(
	NORMSOURCES_COMMON
		${NORMSOURCES_COMMON_ALGORITHM}
		common/argc_argv_faker.cpp
		${NORMSOURCES_COMMON_BATCH}
		${NORMSOURCES_COMMON_BOOST_ADDENDA}
		common/command_executer.cpp
		${NORMSOURCES_COMMON_FILE}
		common/invert_permutation.cpp
		common/logger.cpp
		common/program_exception_wrapper.cpp
		${NORMSOURCES_COMMON_TEST_PREDICATE}
)

set(
	NORMSOURCES_DISPLAY_DISPLAY_COLOUR_SPEC
		display/display_colour_spec/broad_display_colour_spec.cpp
		display/display_colour_spec/display_colour_spec.cpp
)

set(
	NORMSOURCES_DISPLAY_DISPLAY_COLOURER_DETAIL
		display/display_colourer/detail/score_colour_handler.cpp
)

set(
	NORMSOURCES_DISPLAY_DISPLAY_COLOURER
		display/display_colourer/alignment_free_display_colourer.cpp
		${NORMSOURCES_DISPLAY_DISPLAY_COLOURER_DETAIL}
		display/display_colourer/display_colourer.cpp
		display/display_colourer/display_colourer_alignment.cpp
		display/display_colourer/display_colourer_consecutive.cpp
)

set(
	NORMSOURCES_DISPLAY_OPTIONS
		display/options/display_options_block.cpp
		display/options/display_spec.cpp
)

set(
	NORMSOURCES_DISPLAY_VIEWER_PYMOL
		display/viewer/pymol/pymol_tools.cpp
)

set(
	NORMSOURCES_DISPLAY_VIEWER
		display/viewer/chimera_viewer.cpp
		display/viewer/jmol_viewer.cpp
		${NORMSOURCES_DISPLAY_VIEWER_PYMOL}
		display/viewer/pymol_viewer.cpp
		display/viewer/rasmol_style_viewer.cpp
		display/viewer/rasmol_viewer.cpp
		display/viewer/viewer.cpp
)

set(
	NORMSOURCES_DISPLAY
		${NORMSOURCES_DISPLAY_DISPLAY_COLOUR_SPEC}
		${NORMSOURCES_DISPLAY_DISPLAY_COLOURER}
		${NORMSOURCES_DISPLAY_OPTIONS}
		${NORMSOURCES_DISPLAY_VIEWER}
)

set(
	NORMSOURCES_DISPLAY_COLOUR
		display_colour/display_colour.cpp
		display_colour/display_colour_gradient.cpp
		display_colour/display_colour_list.cpp
)

set(
	NORMSOURCES_EXCEPTION
		exception/invalid_argument_exception.cpp
		exception/not_implemented_exception.cpp
		exception/out_of_range_exception.cpp
		exception/runtime_error_exception.cpp
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
	NORMSOURCES_EXECUTABLES_DEVELTEST
		executables/develtest/develtest.cpp
)

set(
	NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT
		executables/snap_judgement/snap_judgement.cpp
)

set(
	NORMSOURCES_EXECUTABLES
		${NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS}
		${NORMSOURCES_EXECUTABLES_CATH_CHECK_PDB}
		${NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN}
		${NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS}
		${NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN}
		${NORMSOURCES_EXECUTABLES_CATH_SSAP}
		${NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE}
		${NORMSOURCES_EXECUTABLES_DEVELTEST}
		${NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT}
)

set(
	NORMSOURCES_FILE_DOMAIN_DEFINITION_LIST
		file/domain_definition_list/domain_definition_list.cpp
)

set(
	NORMSOURCES_FILE_DSSP_WOLF
		file/dssp_wolf/dssp_file.cpp
		file/dssp_wolf/dssp_file_io.cpp
		file/dssp_wolf/tally_residue_ids.cpp
		file/dssp_wolf/wolf_file.cpp
		file/dssp_wolf/wolf_file_io.cpp
)

set(
	NORMSOURCES_FILE_HMMER_SCORES_FILE
		file/hmmer_scores_file/hmmer_scores_entry.cpp
		file/hmmer_scores_file/hmmer_scores_file.cpp
)

set(
	NORMSOURCES_FILE_OPTIONS
		file/options/data_dirs_options_block.cpp
		file/options/data_dirs_spec.cpp
		file/options/data_option.cpp
)

set(
	NORMSOURCES_FILE_PDB
		file/pdb/coarse_element_type.cpp
		file/pdb/dssp_skip_policy.cpp
		file/pdb/pdb.cpp
		file/pdb/pdb_atom.cpp
		file/pdb/pdb_atom_parse_status.cpp
		file/pdb/pdb_base.cpp
		file/pdb/pdb_list.cpp
		file/pdb/pdb_record.cpp
		file/pdb/pdb_residue.cpp
		file/pdb/proximity_calculator.cpp
		file/pdb/read_domain_def_from_pdb.cpp
)

set(
	NORMSOURCES_FILE_PRC_SCORES_FILE
		file/prc_scores_file/prc_scores_entry.cpp
		file/prc_scores_file/prc_scores_file.cpp
)

set(
	NORMSOURCES_FILE_SEC
		file/sec/sec_file.cpp
		file/sec/sec_file_io.cpp
		file/sec/sec_file_record.cpp
)

set(
	NORMSOURCES_FILE_SSAP_SCORES_FILE
		file/ssap_scores_file/ssap_scores_entry.cpp
		file/ssap_scores_file/ssap_scores_file.cpp
)

set(
	NORMSOURCES_FILE
		file/data_file.cpp
		${NORMSOURCES_FILE_DOMAIN_DEFINITION_LIST}
		${NORMSOURCES_FILE_DSSP_WOLF}
		${NORMSOURCES_FILE_HMMER_SCORES_FILE}
		${NORMSOURCES_FILE_OPTIONS}
		${NORMSOURCES_FILE_PDB}
		${NORMSOURCES_FILE_PRC_SCORES_FILE}
		${NORMSOURCES_FILE_SEC}
		${NORMSOURCES_FILE_SSAP_SCORES_FILE}
)

set(
	NORMSOURCES_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS
		options/executable/cath_check_pdb_options/cath_check_pdb_options.cpp
)

set(
	NORMSOURCES_OPTIONS_EXECUTABLE
		${NORMSOURCES_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS}
		options/executable/env_var_option_name_handler.cpp
		options/executable/executable_options.cpp
)

set(
	NORMSOURCES_OPTIONS_OPTIONS_BLOCK
		options/options_block/alignment_input_options_block.cpp
		options/options_block/alignment_input_spec.cpp
		options/options_block/check_pdb_options_block.cpp
		options/options_block/detail_help_options_block.cpp
		options/options_block/ids_options_block.cpp
		options/options_block/misc_help_version_options_block.cpp
		options/options_block/options_block.cpp
		options/options_block/options_block_tester.cpp
		options/options_block/pdb_input_options_block.cpp
		options/options_block/pdb_input_spec.cpp
		options/options_block/superposition_input_options_block.cpp
)

set(
	NORMSOURCES_OPTIONS
		${NORMSOURCES_OPTIONS_EXECUTABLE}
		${NORMSOURCES_OPTIONS_OPTIONS_BLOCK}
)

set(
	NORMSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER
		outputter/alignment_outputter/alignment_outputter.cpp
		outputter/alignment_outputter/alignment_outputter_list.cpp
		outputter/alignment_outputter/cath_aln_ostream_alignment_outputter.cpp
		outputter/alignment_outputter/fasta_ostream_alignment_outputter.cpp
		outputter/alignment_outputter/file_alignment_outputter.cpp
		outputter/alignment_outputter/html_ostream_alignment_outputter.cpp
		outputter/alignment_outputter/ssap_ostream_alignment_outputter.cpp
)

set(
	NORMSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS
		outputter/alignment_outputter_options/alignment_output_options_block.cpp
)

set(
	NORMSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS
		outputter/superposition_output_options/superposition_output_options_block.cpp
)

set(
	NORMSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUTTER
		outputter/superposition_outputter/json_file_superposition_outputter.cpp
		outputter/superposition_outputter/ostream_superposition_outputter.cpp
		outputter/superposition_outputter/pdb_file_superposition_outputter.cpp
		outputter/superposition_outputter/pdb_files_superposition_outputter.cpp
		outputter/superposition_outputter/pymol_file_superposition_outputter.cpp
		outputter/superposition_outputter/pymol_view_superposition_outputter.cpp
		outputter/superposition_outputter/superposition_outputter.cpp
		outputter/superposition_outputter/superposition_outputter_list.cpp
)

set(
	NORMSOURCES_OUTPUTTER
		${NORMSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER}
		${NORMSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS}
		${NORMSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS}
		${NORMSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUTTER}
)

set(
	NORMSOURCES_RESOLVE_HITS_ALGO
		resolve_hits/algo/discont_hits_index_by_start.cpp
		resolve_hits/algo/masked_bests_cacher.cpp
		resolve_hits/algo/scored_arch_proxy.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_FILE_DETAIL
		resolve_hits/file/detail/hmmsearch_aln.cpp
		resolve_hits/file/detail/hmmsearch_parser.cpp
)

set(
	NORMSOURCES_RESOLVE_HITS_FILE
		resolve_hits/file/alnd_rgn.cpp
		resolve_hits/file/cath_id_score_category.cpp
		${NORMSOURCES_RESOLVE_HITS_FILE_DETAIL}
		resolve_hits/file/hits_input_format_tag.cpp
		resolve_hits/file/hmmer_hmmsearch_domtblout.cpp
		resolve_hits/file/hmmer_hmmsearch_out.cpp
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
)

set(
	NORMSOURCES_RESOLVE_HITS_OPTIONS_SPEC
		resolve_hits/options/spec/crh_filter_spec.cpp
		resolve_hits/options/spec/crh_html_spec.cpp
		resolve_hits/options/spec/crh_input_spec.cpp
		resolve_hits/options/spec/crh_output_spec.cpp
		resolve_hits/options/spec/crh_score_spec.cpp
		resolve_hits/options/spec/crh_segment_spec.cpp
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
		resolve_hits/read_and_process_hits/hits_processor/hits_processor.cpp
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
		resolve_hits/hit_resolver.cpp
		resolve_hits/hit_score_type.cpp
		${NORMSOURCES_RESOLVE_HITS_HTML_OUTPUT}
		${NORMSOURCES_RESOLVE_HITS_OPTIONS}
		${NORMSOURCES_RESOLVE_HITS_READ_AND_PROCESS_HITS}
		resolve_hits/scored_hit_arch.cpp
		${NORMSOURCES_RESOLVE_HITS_TRIM}
)

set(
	NORMSOURCES_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY
		scan/detail/check_scan/test_only/alignment_scan_comparison.cpp
		scan/detail/check_scan/test_only/check_scan_on_final_alignment.cpp
		scan/detail/check_scan/test_only/quad_and_rep_criteria_result.cpp
		scan/detail/check_scan/test_only/quad_criteria_result.cpp
)

set(
	NORMSOURCES_SCAN_DETAIL_CHECK_SCAN
		${NORMSOURCES_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY}
)

set(
	NORMSOURCES_SCAN_DETAIL_RES_PAIR
		scan/detail/res_pair/multi_struc_res_rep_pair.cpp
		scan/detail/res_pair/res_pair_core.cpp
		scan/detail/res_pair/single_struc_res_pair.cpp
)

set(
	NORMSOURCES_SCAN_DETAIL_RES_PAIR_DIRN
		scan/detail/res_pair_dirn/res_pair_dirn.cpp
)

set(
	NORMSOURCES_SCAN_DETAIL_STRIDE
		scan/detail/stride/rep_strider.cpp
)

set(
	NORMSOURCES_SCAN_DETAIL
		${NORMSOURCES_SCAN_DETAIL_CHECK_SCAN}
		${NORMSOURCES_SCAN_DETAIL_RES_PAIR}
		${NORMSOURCES_SCAN_DETAIL_RES_PAIR_DIRN}
		${NORMSOURCES_SCAN_DETAIL_STRIDE}
)

set(
	NORMSOURCES_SCAN_SCAN_ACTION
		scan/scan_action/populate_matrix_scan_action.cpp
)

set(
	NORMSOURCES_SCAN_SCAN_TOOLS
		scan/scan_tools/all_vs_all.cpp
		scan/scan_tools/load_and_scan.cpp
		scan/scan_tools/load_and_scan_metrics.cpp
		scan/scan_tools/scan_metrics.cpp
		scan/scan_tools/scan_type.cpp
		scan/scan_tools/single_pair.cpp
)

set(
	NORMSOURCES_SCAN_SPATIAL_INDEX
		scan/spatial_index/spatial_index.cpp
)

set(
	NORMSOURCES_SCAN
		${NORMSOURCES_SCAN_DETAIL}
		scan/quad_criteria.cpp
		scan/res_pair_index_dirn_criterion.cpp
		${NORMSOURCES_SCAN_SCAN_ACTION}
		scan/scan_stride.cpp
		${NORMSOURCES_SCAN_SCAN_TOOLS}
		${NORMSOURCES_SCAN_SPATIAL_INDEX}
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_DETAIL
		score/aligned_pair_score/detail/score_common_coord_handler.cpp
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE
		score/aligned_pair_score/ssap_score/ssap_score_accuracy.cpp
		score/aligned_pair_score/ssap_score/ssap_score_post_processing.cpp
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX
		score/aligned_pair_score/substitution_matrix/blosum62_substitution_matrix.cpp
		score/aligned_pair_score/substitution_matrix/identity_substitution_matrix.cpp
		score/aligned_pair_score/substitution_matrix/match_substitution_matrix.cpp
		score/aligned_pair_score/substitution_matrix/substitution_matrix.cpp
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE
		score/aligned_pair_score/aligned_pair_score.cpp
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_DETAIL}
		score/aligned_pair_score/drmsd_score.cpp
		score/aligned_pair_score/gsas_score.cpp
		score/aligned_pair_score/lddt_score.cpp
		score/aligned_pair_score/length_score.cpp
		score/aligned_pair_score/mi_score.cpp
		score/aligned_pair_score/overlap_score.cpp
		score/aligned_pair_score/pseudo_string_score.cpp
		score/aligned_pair_score/rmsd_score.cpp
		score/aligned_pair_score/sas_score.cpp
		score/aligned_pair_score/sequence_similarity_score.cpp
		score/aligned_pair_score/si_score.cpp
		score/aligned_pair_score/ssap_score.cpp
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE}
		score/aligned_pair_score/structal_score.cpp
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX}
		score/aligned_pair_score/tm_score.cpp
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL
		score/aligned_pair_score_list/detail/aligned_pair_score_list_append.cpp
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER
		score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.cpp
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_READER
		score/aligned_pair_score_list/score_value_list_reader/score_value_reader.cpp
)

set(
	NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST
		score/aligned_pair_score_list/aligned_pair_score_list.cpp
		score/aligned_pair_score_list/aligned_pair_score_list_factory.cpp
		score/aligned_pair_score_list/aligned_pair_score_value_list.cpp
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL}
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER}
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_READER}
)

set(
	NORMSOURCES_SCORE_DETAIL
		score/detail/score_name_helper.cpp
)

set(
	NORMSOURCES_SCORE_HOMCHECK_TOOLS
		score/homcheck_tools/ssap_and_prc.cpp
		score/homcheck_tools/ssaps_and_prcs_of_query.cpp
		score/homcheck_tools/superfamily_of_domain.cpp
)

set(
	NORMSOURCES_SCORE_LENGTH_GETTER
		score/length_getter/geometric_mean_length_getter.cpp
		score/length_getter/length_getter.cpp
		score/length_getter/length_getter_make_clone.cpp
		score/length_getter/length_of_first_getter.cpp
		score/length_getter/length_of_longer_getter.cpp
		score/length_getter/length_of_second_getter.cpp
		score/length_getter/length_of_shorter_getter.cpp
		score/length_getter/mean_length_getter.cpp
		score/length_getter/num_aligned_length_getter.cpp
		score/length_getter/protein_only_length_getter.cpp
		score/length_getter/sym_protein_only_length_getter.cpp
)

set(
	NORMSOURCES_SCORE_PAIR_SCATTER_PLOTTER
		score/pair_scatter_plotter/pair_scatter_plotter.cpp
)

set(
	NORMSOURCES_SCORE_SCORE_CLASSIFICATION_DETAIL
		score/score_classification/detail/score_classn_value_list_name_less.cpp
)

set(
	NORMSOURCES_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE
		score/score_classification/label_pair_is_positive/label_pair_is_positive.cpp
)

set(
	NORMSOURCES_SCORE_SCORE_CLASSIFICATION
		${NORMSOURCES_SCORE_SCORE_CLASSIFICATION_DETAIL}
		${NORMSOURCES_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE}
		score/score_classification/rbf_model.cpp
		score/score_classification/score_classn_value.cpp
		score/score_classification/score_classn_value_better_value.cpp
		score/score_classification/score_classn_value_list.cpp
		score/score_classification/score_classn_value_results_set.cpp
		score/score_classification/value_list_scaling.cpp
)

set(
	NORMSOURCES_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER
		score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter.cpp
		score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_spec.cpp
)

set(
	NORMSOURCES_SCORE_TRUE_POS_FALSE_NEG
		score/true_pos_false_neg/classn_outcome.cpp
		score/true_pos_false_neg/classn_rate_stat.cpp
		score/true_pos_false_neg/classn_stat.cpp
		score/true_pos_false_neg/classn_stat_pair_series.cpp
		score/true_pos_false_neg/classn_stat_pair_series_list.cpp
		${NORMSOURCES_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER}
		score/true_pos_false_neg/named_true_false_pos_neg_list.cpp
		score/true_pos_false_neg/named_true_false_pos_neg_list_list.cpp
		score/true_pos_false_neg/true_false_pos_neg.cpp
		score/true_pos_false_neg/true_false_pos_neg_list.cpp
)

set(
	NORMSOURCES_SCORE
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE}
		${NORMSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST}
		${NORMSOURCES_SCORE_DETAIL}
		${NORMSOURCES_SCORE_HOMCHECK_TOOLS}
		${NORMSOURCES_SCORE_LENGTH_GETTER}
		${NORMSOURCES_SCORE_PAIR_SCATTER_PLOTTER}
		${NORMSOURCES_SCORE_SCORE_CLASSIFICATION}
		${NORMSOURCES_SCORE_TRUE_POS_FALSE_NEG}
)

set(
	NORMSOURCES_SEQ
		seq/seq_arrow.cpp
		seq/seq_seg.cpp
		seq/seq_seg_run.cpp
)

set(
	NORMSOURCES_SSAP_OPTIONS
		ssap/options/cath_ssap_options.cpp
		ssap/options/old_ssap_options_block.cpp
)

set(
	NORMSOURCES_SSAP
		ssap/distance_score_formula.cpp
		${NORMSOURCES_SSAP_OPTIONS}
		ssap/selected_pair.cpp
		ssap/ssap.cpp
		ssap/ssap_scores.cpp
		ssap/windowed_matrix.cpp
)

set(
	NORMSOURCES_STRUCTURE_ACCESSIBILITY_CALC
		structure/accessibility_calc/dssp_accessibility.cpp
)

set(
	NORMSOURCES_STRUCTURE_BIOPLIB_FACADE
		structure/bioplib_facade/bioplib_interface.cpp
		structure/bioplib_facade/bioplib_pdb.cpp
)

set(
	NORMSOURCES_STRUCTURE_ENTRY_QUERIER
		structure/entry_querier/entry_querier.cpp
		structure/entry_querier/residue_querier.cpp
		structure/entry_querier/sec_struc_querier.cpp
)

set(
	NORMSOURCES_STRUCTURE_GEOMETRY
		structure/geometry/angle.cpp
		structure/geometry/coord.cpp
		structure/geometry/coord_list.cpp
		structure/geometry/pca.cpp
		structure/geometry/restrict_to_single_linkage_extension.cpp
		structure/geometry/rotation.cpp
)

set(
	NORMSOURCES_STRUCTURE_PROTEIN_PROTEIN_LOADER
		structure/protein/protein_loader/protein_list_loader.cpp
)

set(
	NORMSOURCES_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET
		structure/protein/protein_source_file_set/protein_file_combn.cpp
		structure/protein/protein_source_file_set/protein_from_pdb.cpp
		structure/protein/protein_source_file_set/protein_from_pdb_and_calc.cpp
		structure/protein/protein_source_file_set/protein_from_pdb_and_dssp_and_calc.cpp
		structure/protein/protein_source_file_set/protein_from_pdb_dssp_and_sec.cpp
		structure/protein/protein_source_file_set/protein_from_wolf_and_sec.cpp
		structure/protein/protein_source_file_set/protein_source_file_set.cpp
)

set(
	NORMSOURCES_STRUCTURE_PROTEIN
		structure/protein/amino_acid.cpp
		structure/protein/dna_atom.cpp
		structure/protein/protein.cpp
		structure/protein/protein_io.cpp
		structure/protein/protein_list.cpp
		${NORMSOURCES_STRUCTURE_PROTEIN_PROTEIN_LOADER}
		${NORMSOURCES_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET}
		structure/protein/residue.cpp
		structure/protein/sec_struc.cpp
		structure/protein/sec_struc_planar_angles.cpp
		structure/protein/sec_struc_type.cpp
)

set(
	NORMSOURCES_STRUCTURE_SEC_STRUC_CALC_DSSP
		structure/sec_struc_calc/dssp/bifur_hbond_list.cpp
		structure/sec_struc_calc/dssp/dssp_hbond_calc.cpp
		structure/sec_struc_calc/dssp/dssp_ss_calc.cpp
)

set(
	NORMSOURCES_STRUCTURE_SEC_STRUC_CALC_SEC
		structure/sec_struc_calc/sec/sec_calc.cpp
)

set(
	NORMSOURCES_STRUCTURE_SEC_STRUC_CALC
		${NORMSOURCES_STRUCTURE_SEC_STRUC_CALC_DSSP}
		${NORMSOURCES_STRUCTURE_SEC_STRUC_CALC_SEC}
)

set(
	NORMSOURCES_STRUCTURE_VIEW_CACHE_FILTER_DETAIL
		structure/view_cache/filter/detail/filter_vs_full_score_less.cpp
)

set(
	NORMSOURCES_STRUCTURE_VIEW_CACHE_FILTER
		${NORMSOURCES_STRUCTURE_VIEW_CACHE_FILTER_DETAIL}
		structure/view_cache/filter/filter_vs_full_score.cpp
		structure/view_cache/filter/filter_vs_full_score_list.cpp
)

set(
	NORMSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL
		structure/view_cache/index/detail/vcie_match_criteria.cpp
)

set(
	NORMSOURCES_STRUCTURE_VIEW_CACHE_INDEX
		${NORMSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL}
		structure/view_cache/index/quad_find_action.cpp
		structure/view_cache/index/quad_find_action_check.cpp
		structure/view_cache/index/view_cache_index.cpp
		structure/view_cache/index/view_cache_index_entry.cpp
)

set(
	NORMSOURCES_STRUCTURE_VIEW_CACHE
		${NORMSOURCES_STRUCTURE_VIEW_CACHE_FILTER}
		${NORMSOURCES_STRUCTURE_VIEW_CACHE_INDEX}
		structure/view_cache/view_cache.cpp
		structure/view_cache/view_cache_list.cpp
)

set(
	NORMSOURCES_STRUCTURE
		${NORMSOURCES_STRUCTURE_ACCESSIBILITY_CALC}
		${NORMSOURCES_STRUCTURE_BIOPLIB_FACADE}
		${NORMSOURCES_STRUCTURE_ENTRY_QUERIER}
		${NORMSOURCES_STRUCTURE_GEOMETRY}
		${NORMSOURCES_STRUCTURE_PROTEIN}
		${NORMSOURCES_STRUCTURE_SEC_STRUC_CALC}
		${NORMSOURCES_STRUCTURE_VIEW_CACHE}
)

set(
	NORMSOURCES_SUPERPOSITION_IO
		superposition/io/superposition_io.cpp
)

set(
	NORMSOURCES_SUPERPOSITION
		${NORMSOURCES_SUPERPOSITION_IO}
		superposition/superpose_orderer.cpp
		superposition/superposition.cpp
		superposition/superposition_content_spec.cpp
		superposition/superposition_context.cpp
		superposition/supn_regions_context.cpp
)

set(
	NORMSOURCES_TEST
		test/global_test_constants.cpp
)

set(
	TESTSOURCES_ACQUIRER_ALIGNMENT_ACQUIRER
		acquirer/alignment_acquirer/alignment_acquirer_test.cpp
		acquirer/alignment_acquirer/cora_aln_file_alignment_acquirer_test.cpp
		acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer_test.cpp
		acquirer/alignment_acquirer/residue_name_alignment_acquirer_test.cpp
		acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer_test.cpp
		acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer_test.cpp
)

set(
	TESTSOURCES_ACQUIRER_PDBS_ACQUIRER
		acquirer/pdbs_acquirer/file_list_pdbs_acquirer_test.cpp
		acquirer/pdbs_acquirer/istream_pdbs_acquirer_test.cpp
		acquirer/pdbs_acquirer/pdbs_acquirer_test.cpp
)

set(
	TESTSOURCES_ACQUIRER_SELECTION_POLICY_ACQUIRER
		acquirer/selection_policy_acquirer/selection_policy_acquirer_test.cpp
)

set(
	TESTSOURCES_ACQUIRER_SUPERPOSITION_ACQUIRER
		acquirer/superposition_acquirer/align_based_superposition_acquirer_test.cpp
		acquirer/superposition_acquirer/superposition_acquirer_test.cpp
)

set(
	TESTSOURCES_ACQUIRER
		${TESTSOURCES_ACQUIRER_ALIGNMENT_ACQUIRER}
		${TESTSOURCES_ACQUIRER_PDBS_ACQUIRER}
		${TESTSOURCES_ACQUIRER_SELECTION_POLICY_ACQUIRER}
		${TESTSOURCES_ACQUIRER_SUPERPOSITION_ACQUIRER}
)

set(
	TESTSOURCES_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY
		alignment/common_residue_selection_policy/common_residue_selection_policy_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER
		alignment/dyn_prog_align/detail/string_aligner/string_aligner_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL
		alignment/dyn_prog_align/detail/return_path_matrix_test.cpp
		alignment/dyn_prog_align/detail/score_accumulation_matrix_test.cpp
		${TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER}
)

set(
	TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE
		alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source_test.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/entry_querier_dyn_prog_score_source_test.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source_test.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source_test.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source_test.cpp
		alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_TEST
		alignment/dyn_prog_align/test/dyn_prog_score_source_fixture.cpp
)

set(
	TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN
		${TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DETAIL}
		${TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE}
		alignment/dyn_prog_align/ssap_code_dyn_prog_aligner_test.cpp
		alignment/dyn_prog_align/std_dyn_prog_aligner_test.cpp
		${TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN_TEST}
)

set(
	TESTSOURCES_ALIGNMENT_GAP
		alignment/gap/alignment_gap_test.cpp
		alignment/gap/gap_penalty_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_IO_OUTPUTTER
		alignment/io/outputter/horiz_align_outputter_test.cpp
		alignment/io/outputter/html_align_outputter_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_IO
		alignment/io/align_scaffold_test.cpp
		alignment/io/alignment_io_test.cpp
		${TESTSOURCES_ALIGNMENT_IO_OUTPUTTER}
)

set(
	TESTSOURCES_ALIGNMENT_REFINER_DETAIL
		alignment/refiner/detail/alignment_split_list_test.cpp
		alignment/refiner/detail/alignment_split_mapping_test.cpp
		alignment/refiner/detail/alignment_split_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_REFINER
		alignment/refiner/alignment_refiner_test.cpp
		${TESTSOURCES_ALIGNMENT_REFINER_DETAIL}
		alignment/refiner/indexed_refiner_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL
		alignment/residue_name_align/detail/residue_name_align_map_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN
		${TESTSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL}
		alignment/residue_name_align/residue_name_aligner_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_RESIDUE_SCORE
		alignment/residue_score/alignment_residue_scores_test.cpp
		alignment/residue_score/residue_scorer_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT_TEST
		alignment/test/alignment_fixture.cpp
)

set(
	TESTSOURCES_ALIGNMENT_TOOLS
		alignment/tools/alignment_breaks_test.cpp
)

set(
	TESTSOURCES_ALIGNMENT
		alignment/alignment_action_test.cpp
		alignment/alignment_context_test.cpp
		alignment/alignment_coord_extractor_test.cpp
		alignment/alignment_row_test.cpp
		alignment/alignment_test.cpp
		${TESTSOURCES_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY}
		${TESTSOURCES_ALIGNMENT_DYN_PROG_ALIGN}
		${TESTSOURCES_ALIGNMENT_GAP}
		${TESTSOURCES_ALIGNMENT_IO}
		${TESTSOURCES_ALIGNMENT_REFINER}
		${TESTSOURCES_ALIGNMENT_RESIDUE_NAME_ALIGN}
		${TESTSOURCES_ALIGNMENT_RESIDUE_SCORE}
		${TESTSOURCES_ALIGNMENT_TEST}
		${TESTSOURCES_ALIGNMENT_TOOLS}
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
	TESTSOURCES_COMMON_ALGORITHM
		common/algorithm/constexpr_find_test.cpp
		common/algorithm/constexpr_floor_test.cpp
		common/algorithm/constexpr_integer_rounding_test.cpp
		common/algorithm/constexpr_is_uniq_test.cpp
		common/algorithm/constexpr_modulo_fns_test.cpp
		common/algorithm/for_n_test.cpp
		common/algorithm/transform_build_test.cpp
		common/algorithm/variadic_and_test.cpp
)

set(
	TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR
		common/boost_addenda/range/adaptor/adaptor_test.cpp
)

set(
	TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR
		common/boost_addenda/range/utility/iterator/cross_itr_test.cpp
)

set(
	TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE_UTILITY
		${TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR}
)

set(
	TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE
		${TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR}
		${TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE_UTILITY}
)

set(
	TESTSOURCES_COMMON_BOOST_ADDENDA_TRIBOOL
		common/boost_addenda/tribool/tribool_test.cpp
)

set(
	TESTSOURCES_COMMON_BOOST_ADDENDA
		${TESTSOURCES_COMMON_BOOST_ADDENDA_RANGE}
		${TESTSOURCES_COMMON_BOOST_ADDENDA_TRIBOOL}
)

set(
	TESTSOURCES_COMMON_CLONE
		common/clone/clone_ptr_test.cpp
)

set(
	TESTSOURCES_COMMON_CONTAINER
		common/container/id_of_string_test.cpp
)

set(
	TESTSOURCES_COMMON_CPP17
		common/cpp17/apply_test.cpp
)

set(
	TESTSOURCES_COMMON_FILE
		common/file/open_fstream_test.cpp
		common/file/simple_file_read_write_test.cpp
		common/file/temp_file_test.cpp
)

set(
	TESTSOURCES_COMMON_FUNCTION
		common/function/multiply_args_test.cpp
)

set(
	TESTSOURCES_COMMON_METAPROGRAMMING
		common/metaprogramming/append_template_params_into_first_wrapper_test.cpp
		common/metaprogramming/change_template_subwrappers_test.cpp
		common/metaprogramming/change_template_wrapper_test.cpp
		common/metaprogramming/combine_params_lists_with_template_list_test.cpp
)

set(
	TESTSOURCES_COMMON_RAPIDJSON_ADDENDA
		common/rapidjson_addenda/rapidjson_writer_test.cpp
		common/rapidjson_addenda/string_of_rapidjson_write_test.cpp
)

set(
	TESTSOURCES_COMMON_STRING
		common/string/booled_to_string_test.cpp
		common/string/string_parse_tools_test.cpp
		common/string/sub_string_parser_test.cpp
)

set(
	TESTSOURCES_COMMON_TEST_PREDICATE
		common/test_predicate/files_equal_test.cpp
		common/test_predicate/istreams_equal_test.cpp
)

set(
	TESTSOURCES_COMMON_TUPLE
		common/tuple/make_tuple_with_skips_test.cpp
		common/tuple/mins_maxs_tuple_pair_mins_maxs_element_test.cpp
		common/tuple/tuple_increment_test.cpp
		common/tuple/tuple_lattice_index_test.cpp
		common/tuple/tuple_mins_maxs_element_test.cpp
		common/tuple/tuple_multiply_args_test.cpp
		common/tuple/tuple_subtract_test.cpp
		common/tuple/tuple_within_range_test.cpp
)

set(
	TESTSOURCES_COMMON_TYPE_TRAITS
		common/type_traits/is_tuple_test.cpp
)

set(
	TESTSOURCES_COMMON
		${TESTSOURCES_COMMON_ALGORITHM}
		common/argc_argv_faker_test.cpp
		${TESTSOURCES_COMMON_BOOST_ADDENDA}
		${TESTSOURCES_COMMON_CLONE}
		common/command_executer_test.cpp
		${TESTSOURCES_COMMON_CONTAINER}
		${TESTSOURCES_COMMON_CPP17}
		common/difference_test.cpp
		${TESTSOURCES_COMMON_FILE}
		${TESTSOURCES_COMMON_FUNCTION}
		common/invert_permutation_test.cpp
		common/less_than_helper_test.cpp
		common/logger_test.cpp
		${TESTSOURCES_COMMON_METAPROGRAMMING}
		common/program_exception_wrapper_test.cpp
		${TESTSOURCES_COMMON_RAPIDJSON_ADDENDA}
		${TESTSOURCES_COMMON_STRING}
		common/temp_check_offset_1_test.cpp
		${TESTSOURCES_COMMON_TEST_PREDICATE}
		${TESTSOURCES_COMMON_TUPLE}
		common/type_to_string_test.cpp
		${TESTSOURCES_COMMON_TYPE_TRAITS}
)

set(
	TESTSOURCES_DISPLAY_DISPLAY_COLOUR_SPEC
		display/display_colour_spec/broad_display_colour_spec_test.cpp
		display/display_colour_spec/display_colour_spec_test.cpp
)

set(
	TESTSOURCES_DISPLAY_DISPLAY_COLOURER
		display/display_colourer/display_colourer_test.cpp
)

set(
	TESTSOURCES_DISPLAY_OPTIONS
		display/options/display_options_block_test.cpp
		display/options/display_spec_test.cpp
)

set(
	TESTSOURCES_DISPLAY_VIEWER_PYMOL
		display/viewer/pymol/pymol_tools_test.cpp
)

set(
	TESTSOURCES_DISPLAY_VIEWER
		display/viewer/jmol_viewer_test.cpp
		${TESTSOURCES_DISPLAY_VIEWER_PYMOL}
		display/viewer/pymol_viewer_test.cpp
		display/viewer/rasmol_style_viewer_test.cpp
		display/viewer/rasmol_viewer_test.cpp
		display/viewer/viewer_test.cpp
)

set(
	TESTSOURCES_DISPLAY
		${TESTSOURCES_DISPLAY_DISPLAY_COLOUR_SPEC}
		${TESTSOURCES_DISPLAY_DISPLAY_COLOURER}
		${TESTSOURCES_DISPLAY_OPTIONS}
		${TESTSOURCES_DISPLAY_VIEWER}
)

set(
	TESTSOURCES_DISPLAY_COLOUR
		display_colour/display_colour_gradient_test.cpp
		display_colour/display_colour_list_test.cpp
		display_colour/display_colour_test.cpp
)

set(
	TESTSOURCES_EXCEPTION
		exception/exception_is_equivalent_test.cpp
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
	TESTSOURCES_FILE_DSSP_WOLF
		file/dssp_wolf/dssp_file_test.cpp
		file/dssp_wolf/tally_residue_ids_test.cpp
)

set(
	TESTSOURCES_FILE_HMMER_SCORES_FILE
		file/hmmer_scores_file/hmmer_scores_entry_test.cpp
		file/hmmer_scores_file/hmmer_scores_file_test.cpp
)

set(
	TESTSOURCES_FILE_OPTIONS
		file/options/data_dirs_options_block_test.cpp
		file/options/data_dirs_spec_test.cpp
)

set(
	TESTSOURCES_FILE_PDB
		file/pdb/coarse_element_type_test.cpp
		file/pdb/element_type_string_test.cpp
		file/pdb/pdb_atom_test.cpp
		file/pdb/pdb_base_test.cpp
		file/pdb/pdb_list_test.cpp
		file/pdb/pdb_residue_test.cpp
		file/pdb/pdb_test.cpp
		file/pdb/proximity_calculator_test.cpp
)

set(
	TESTSOURCES_FILE_PRC_SCORES_FILE
		file/prc_scores_file/prc_scores_entry_test.cpp
		file/prc_scores_file/prc_scores_file_test.cpp
)

set(
	TESTSOURCES_FILE_SEC
		file/sec/sec_file_test.cpp
)

set(
	TESTSOURCES_FILE_SSAP_SCORES_FILE
		file/ssap_scores_file/ssap_scores_entry_test.cpp
		file/ssap_scores_file/ssap_scores_file_test.cpp
)

set(
	TESTSOURCES_FILE
		${TESTSOURCES_FILE_DSSP_WOLF}
		${TESTSOURCES_FILE_HMMER_SCORES_FILE}
		${TESTSOURCES_FILE_OPTIONS}
		${TESTSOURCES_FILE_PDB}
		${TESTSOURCES_FILE_PRC_SCORES_FILE}
		${TESTSOURCES_FILE_SEC}
		${TESTSOURCES_FILE_SSAP_SCORES_FILE}
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
	TESTSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER
		outputter/alignment_outputter/alignment_outputter_list_test.cpp
)

set(
	TESTSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS
		outputter/alignment_outputter_options/alignment_output_options_block_test.cpp
)

set(
	TESTSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS
		outputter/superposition_output_options/superposition_output_options_block_test.cpp
)

set(
	TESTSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUTTER
		outputter/superposition_outputter/json_file_superposition_outputter_test.cpp
		outputter/superposition_outputter/ostream_superposition_outputter_test.cpp
		outputter/superposition_outputter/pdb_file_superposition_outputter_test.cpp
		outputter/superposition_outputter/pdb_files_superposition_outputter_test.cpp
		outputter/superposition_outputter/pymol_file_superposition_outputter_test.cpp
		outputter/superposition_outputter/pymol_view_superposition_outputter_test.cpp
		outputter/superposition_outputter/superposition_outputter_list_test.cpp
		outputter/superposition_outputter/superposition_outputter_test.cpp
)

set(
	TESTSOURCES_OUTPUTTER
		${TESTSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER}
		${TESTSOURCES_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS}
		${TESTSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS}
		${TESTSOURCES_OUTPUTTER_SUPERPOSITION_OUTPUTTER}
)

set(
	TESTSOURCES_RESOLVE_HITS_ALGO
		resolve_hits/algo/masked_bests_cache_test.cpp
)

set(
	TESTSOURCES_RESOLVE_HITS_FILE
		resolve_hits/file/cath_id_score_category_test.cpp
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
)

set(
	TESTSOURCES_RESOLVE_HITS_OPTIONS
		resolve_hits/options/crh_options_test.cpp
		${TESTSOURCES_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK}
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
		resolve_hits/full_hit_list_test.cpp
		resolve_hits/full_hit_test.cpp
		resolve_hits/hit_extras_test.cpp
		resolve_hits/hit_resolver_test.cpp
		resolve_hits/hit_test.cpp
		${TESTSOURCES_RESOLVE_HITS_HTML_OUTPUT}
		${TESTSOURCES_RESOLVE_HITS_OPTIONS}
		${TESTSOURCES_RESOLVE_HITS_TEST}
		${TESTSOURCES_RESOLVE_HITS_TRIM}
)

set(
	TESTSOURCES_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY
		scan/detail/check_scan/test_only/check_scan_on_final_alignment_test.cpp
)

set(
	TESTSOURCES_SCAN_DETAIL_CHECK_SCAN
		${TESTSOURCES_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY}
)

set(
	TESTSOURCES_SCAN_DETAIL_RES_PAIR
		scan/detail/res_pair/multi_struc_res_rep_pair_list_test.cpp
		scan/detail/res_pair/multi_struc_res_rep_pair_test.cpp
		scan/detail/res_pair/res_pair_core_test.cpp
		scan/detail/res_pair/single_struc_res_pair_list_test.cpp
		scan/detail/res_pair/single_struc_res_pair_test.cpp
)

set(
	TESTSOURCES_SCAN_DETAIL_SCAN_ACTION
		scan/detail/scan_action/align_scan_action_test.cpp
		scan/detail/scan_action/scan_multi_action_test.cpp
)

set(
	TESTSOURCES_SCAN_DETAIL_SCAN_INDEX_STORE
		scan/detail/scan_index_store/scan_index_hash_store_test.cpp
		scan/detail/scan_index_store/scan_index_lattice_store_test.cpp
		scan/detail/scan_index_store/scan_index_vector_store_test.cpp
)

set(
	TESTSOURCES_SCAN_DETAIL_STRIDE
		scan/detail/stride/co_stride_test.cpp
		scan/detail/stride/rep_strider_test.cpp
		scan/detail/stride/roled_scan_stride_test.cpp
)

set(
	TESTSOURCES_SCAN_DETAIL
		${TESTSOURCES_SCAN_DETAIL_CHECK_SCAN}
		${TESTSOURCES_SCAN_DETAIL_RES_PAIR}
		${TESTSOURCES_SCAN_DETAIL_SCAN_ACTION}
		${TESTSOURCES_SCAN_DETAIL_SCAN_INDEX_STORE}
		scan/detail/scan_structure_data_test.cpp
		${TESTSOURCES_SCAN_DETAIL_STRIDE}
)

set(
	TESTSOURCES_SCAN_RES_PAIR_KEYER
		scan/res_pair_keyer/res_pair_keyer_test.cpp
)

set(
	TESTSOURCES_SCAN_SCAN_TOOLS
		scan/scan_tools/all_vs_all_test.cpp
		scan/scan_tools/load_and_scan_metrics_test.cpp
		scan/scan_tools/load_and_scan_test.cpp
)

set(
	TESTSOURCES_SCAN_SPATIAL_INDEX
		scan/spatial_index/spatial_index_test.cpp
)

set(
	TESTSOURCES_SCAN
		${TESTSOURCES_SCAN_DETAIL}
		scan/quad_criteria_test.cpp
		${TESTSOURCES_SCAN_RES_PAIR_KEYER}
		scan/scan_index_test.cpp
		scan/scan_policy_test.cpp
		scan/scan_query_set_test.cpp
		scan/scan_stride_test.cpp
		${TESTSOURCES_SCAN_SCAN_TOOLS}
		${TESTSOURCES_SCAN_SPATIAL_INDEX}
)

set(
	TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_DETAIL
		score/aligned_pair_score/detail/score_common_coord_handler_test.cpp
)

set(
	TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE
		score/aligned_pair_score/ssap_score/ssap_score_accuracy_test.cpp
		score/aligned_pair_score/ssap_score/ssap_score_post_processing_test.cpp
)

set(
	TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX
		score/aligned_pair_score/substitution_matrix/substitution_matrix_test.cpp
)

set(
	TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE
		score/aligned_pair_score/aligned_pair_score_test.cpp
		${TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_DETAIL}
		score/aligned_pair_score/rmsd_score_test.cpp
		${TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE}
		${TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX}
)

set(
	TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST
		score/aligned_pair_score_list/aligned_pair_score_list_factory_test.cpp
		score/aligned_pair_score_list/aligned_pair_score_list_test.cpp
)

set(
	TESTSOURCES_SCORE_HOMCHECK_TOOLS
		score/homcheck_tools/ssap_and_prc_test.cpp
		score/homcheck_tools/ssaps_and_prcs_of_query_test.cpp
		score/homcheck_tools/superfamily_of_domain_test.cpp
)

set(
	TESTSOURCES_SCORE_LENGTH_GETTER
		score/length_getter/length_getter_test.cpp
)

set(
	TESTSOURCES_SCORE_SCORE_CLASSIFICATION_DETAIL
		score/score_classification/detail/score_classn_value_list_name_less_test.cpp
)

set(
	TESTSOURCES_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE
		score/score_classification/label_pair_is_positive/label_pair_is_positive_test.cpp
)

set(
	TESTSOURCES_SCORE_SCORE_CLASSIFICATION
		${TESTSOURCES_SCORE_SCORE_CLASSIFICATION_DETAIL}
		${TESTSOURCES_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE}
		score/score_classification/rbf_model_test.cpp
		score/score_classification/score_classn_value_better_value_test.cpp
		score/score_classification/score_classn_value_list_test.cpp
		score/score_classification/score_classn_value_results_set_test.cpp
		score/score_classification/score_classn_value_test.cpp
)

set(
	TESTSOURCES_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER
		score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_test.cpp
)

set(
	TESTSOURCES_SCORE_TRUE_POS_FALSE_NEG
		score/true_pos_false_neg/classn_rate_stat_test.cpp
		score/true_pos_false_neg/classn_stat_pair_series_list_test.cpp
		score/true_pos_false_neg/classn_stat_pair_series_test.cpp
		${TESTSOURCES_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER}
		score/true_pos_false_neg/classn_stat_test.cpp
		score/true_pos_false_neg/named_true_false_pos_neg_list_test.cpp
		score/true_pos_false_neg/true_false_pos_neg_list_test.cpp
		score/true_pos_false_neg/true_false_pos_neg_test.cpp
)

set(
	TESTSOURCES_SCORE
		${TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE}
		${TESTSOURCES_SCORE_ALIGNED_PAIR_SCORE_LIST}
		${TESTSOURCES_SCORE_HOMCHECK_TOOLS}
		${TESTSOURCES_SCORE_LENGTH_GETTER}
		${TESTSOURCES_SCORE_SCORE_CLASSIFICATION}
		${TESTSOURCES_SCORE_TRUE_POS_FALSE_NEG}
)

set(
	TESTSOURCES_SEQ
		seq/seq_seg_run_test.cpp
		seq/seq_seg_test.cpp
)

set(
	TESTSOURCES_SSAP_OPTIONS
		ssap/options/cath_ssap_options_test.cpp
		ssap/options/old_ssap_options_block_test.cpp
)

set(
	TESTSOURCES_SSAP
		ssap/distance_score_formula_test.cpp
		${TESTSOURCES_SSAP_OPTIONS}
		ssap/selected_pair_test.cpp
		ssap/ssap_scores_test.cpp
		ssap/ssap_test.cpp
		ssap/windowed_matrix_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_ACCESSIBILITY_CALC
		structure/accessibility_calc/dssp_accessibility_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_BIOPLIB_FACADE
		structure/bioplib_facade/bioplib_pdb_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_ENTRY_QUERIER
		structure/entry_querier/entry_querier_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_GEOMETRY
		structure/geometry/angle_test.cpp
		structure/geometry/coord_list_test.cpp
		structure/geometry/coord_test.cpp
		structure/geometry/orientation_covering_test.cpp
		structure/geometry/pca_test.cpp
		structure/geometry/quat_rot_test.cpp
		structure/geometry/restrict_to_single_linkage_extension_test.cpp
		structure/geometry/rotation_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET
		structure/protein/protein_source_file_set/protein_source_file_set_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_PROTEIN
		structure/protein/amino_acid_test.cpp
		structure/protein/protein_list_test.cpp
		${TESTSOURCES_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET}
		structure/protein/protein_test.cpp
		structure/protein/residue_test.cpp
		structure/protein/sec_struc_planar_angles_test.cpp
		structure/protein/sec_struc_test.cpp
		structure/protein/sec_struc_type_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST
		structure/sec_struc_calc/dssp/test/dssp_dupl_fixture.cpp
)

set(
	TESTSOURCES_STRUCTURE_SEC_STRUC_CALC_DSSP
		structure/sec_struc_calc/dssp/bifur_hbond_list_test.cpp
		structure/sec_struc_calc/dssp/dssp_hbond_calc_test.cpp
		structure/sec_struc_calc/dssp/dssp_ss_calc_test.cpp
		${TESTSOURCES_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST}
)

set(
	TESTSOURCES_STRUCTURE_SEC_STRUC_CALC_SEC
		structure/sec_struc_calc/sec/sec_calc_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_SEC_STRUC_CALC
		${TESTSOURCES_STRUCTURE_SEC_STRUC_CALC_DSSP}
		${TESTSOURCES_STRUCTURE_SEC_STRUC_CALC_SEC}
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_DETAIL_PLATE
		structure/view_cache/detail/plate/plate_scan_test.cpp
		structure/view_cache/detail/plate/rod_cache_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_DETAIL
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_DETAIL_PLATE}
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_FILTER_DETAIL
		structure/view_cache/filter/detail/filter_vs_full_score_less_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_FILTER
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_FILTER_DETAIL}
		structure/view_cache/filter/filter_vs_full_score_list_test.cpp
		structure/view_cache/filter/filter_vs_full_score_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL
		structure/view_cache/index/detail/dims/detail/view_cache_index_dim_linear_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL}
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS}
		structure/view_cache/index/detail/vcie_match_criteria_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX_DETAIL}
		structure/view_cache/index/view_cache_index_entry_test.cpp
		structure/view_cache/index/view_cache_index_test.cpp
)

set(
	TESTSOURCES_STRUCTURE_VIEW_CACHE
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_DETAIL}
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_FILTER}
		${TESTSOURCES_STRUCTURE_VIEW_CACHE_INDEX}
)

set(
	TESTSOURCES_STRUCTURE
		${TESTSOURCES_STRUCTURE_ACCESSIBILITY_CALC}
		${TESTSOURCES_STRUCTURE_BIOPLIB_FACADE}
		${TESTSOURCES_STRUCTURE_ENTRY_QUERIER}
		${TESTSOURCES_STRUCTURE_GEOMETRY}
		structure/get_residue_names_test.cpp
		${TESTSOURCES_STRUCTURE_PROTEIN}
		${TESTSOURCES_STRUCTURE_SEC_STRUC_CALC}
		${TESTSOURCES_STRUCTURE_VIEW_CACHE}
)

set(
	TESTSOURCES_SUPERPOSITION_IO
		superposition/io/superposition_io_test.cpp
)

set(
	TESTSOURCES_SUPERPOSITION
		${TESTSOURCES_SUPERPOSITION_IO}
		superposition/superpose_orderer_test.cpp
		superposition/superposition_context_test.cpp
		superposition/superposition_test.cpp
)

set(
	TESTSOURCES_TEST
		test/superposition_fixture.cpp
)
##### END OF AUTO-GENERATED #####

#####################################################################
#####################################################################

set(
	STD_TESTSOURCES
		${TESTSOURCES_ACQUIRER}
		${TESTSOURCES_ALIGNMENT}
		${TESTSOURCES_CATH_REFINE_ALIGN}
		${TESTSOURCES_CATH_SUPERPOSE}
		${TESTSOURCES_CHOPPING}
		${TESTSOURCES_COMMON}
		${TESTSOURCES_DISPLAY}
		${TESTSOURCES_DISPLAY_COLOUR}
		${TESTSOURCES_EXCEPTION}
		${TESTSOURCES_FILE}
		${TESTSOURCES_OPTIONS}
		${TESTSOURCES_OUTPUTTER}
		${TESTSOURCES_RESOLVE_HITS}
		${TESTSOURCES_SCAN}
		${TESTSOURCES_SCORE}
		${TESTSOURCES_SEQ}
		${TESTSOURCES_SSAP}
		${TESTSOURCES_STRUCTURE}
		${TESTSOURCES_SUPERPOSITION}
)

#####################################################################

add_library(
	std_testsources
	OBJECT
		${STD_TESTSOURCES}
)

#####################################################################
#####################################################################

add_library(
	normsources_cath_assign_domains
		${NORMSOURCES_CATH_ASSIGN_DOMAINS}
)

add_library(
	normsources_cath_refine_align
		${NORMSOURCES_CATH_REFINE_ALIGN}
)

add_library(
	normsources_cath_score_align
		${NORMSOURCES_CATH_SCORE_ALIGN}
)

add_library(
	normsources_cath_superpose
		${NORMSOURCES_CATH_SUPERPOSE}
)

add_library(
	normsources_chopping
		${NORMSOURCES_CHOPPING}
)

add_library(
	normsources_common
		${NORMSOURCES_COMMON}
)

add_library(
	normsources_exception
		${NORMSOURCES_EXCEPTION}
)

add_library(
	normsources_resolve_hits
		${NORMSOURCES_RESOLVE_HITS}
)

add_library(
	normsources_seq
		${NORMSOURCES_SEQ}
)

add_library(
	normsources_test
		${NORMSOURCES_TEST}
)

#####################################################################

add_library(
	std_normsources_remain
		${NORMSOURCES_ACQUIRER}
		${NORMSOURCES_ALIGNMENT}
		${NORMSOURCES_BIOCORE}
		${NORMSOURCES_CHOPPING}
		${NORMSOURCES_DISPLAY}
		${NORMSOURCES_DISPLAY_COLOUR}
		${NORMSOURCES_FILE}
		${NORMSOURCES_OPTIONS}
		${NORMSOURCES_OUTPUTTER}
		${NORMSOURCES_SCAN}
		${NORMSOURCES_SCORE}
		${NORMSOURCES_SSAP}
		${NORMSOURCES_STRUCTURE}
		${NORMSOURCES_SUPERPOSITION}
)

#####################################################################



#####################################################################
#####################################################################

add_executable(
	build-test
		$<TARGET_OBJECTS:std_testsources>
		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
)

add_executable(
	cath-assign-domains
		${NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS}
)

add_executable(
	cath-refine-align
		${NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN}
)

add_executable(
	cath-resolve-hits
		${NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS}
		${NORMSOURCES_DISPLAY_COLOUR}
		${NORMSOURCES_EXCEPTION}
		${NORMSOURCES_OPTIONS}
)

add_executable(
	cath-score-align
		${NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN}
)

add_executable(
	cath-ssap
		${NORMSOURCES_EXECUTABLES_CATH_SSAP}
)

add_executable(
	cath-superpose
		${NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE}
)

add_executable(
	mod-test-chopping
		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
		${TESTSOURCES_CHOPPING}
		${NORMSOURCES_BIOCORE}
)

add_executable(
	mod-test-seq
		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
		${TESTSOURCES_SEQ}
		#${NORMSOURCES_BIOCORE}
)

# add_executable(
# 	mod-test-common
# 		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
# 		${TESTSOURCES_COMMON}
# )
#
# add_executable(
# 	mod-test-exception
# 		${TESTSOURCES_EXECUTABLES_BUILD_TEST}
# 		${TESTSOURCES_EXCEPTION}
# )

#####################################################################
#####################################################################

target_link_libraries( build-test          normsources_cath_refine_align normsources_cath_superpose normsources_resolve_hits normsources_test normsources_seq std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-assign-domains normsources_cath_assign_domains                                                                                    std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-refine-align   normsources_cath_refine_align                                                                                      std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-resolve-hits   normsources_resolve_hits                                                                           normsources_seq                        normsources_common normsources_exception          pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-score-align    normsources_cath_score_align                                                                                       std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-ssap                                                                                                                              std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( cath-superpose      normsources_cath_superpose                                                                                         std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )

target_link_libraries( mod-test-chopping    normsources_test normsources_chopping normsources_common normsources_exception pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
target_link_libraries( mod-test-seq         normsources_test normsources_seq      normsources_common normsources_exception pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
# target_link_libraries( mod-test-common    normsources_test                      normsources_common normsources_exception pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
# target_link_libraries( mod-test-exception normsources_test                                         normsources_exception pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )

add_dependencies( build-test          built_bioplib_libraries )
add_dependencies( cath-assign-domains built_bioplib_libraries )
add_dependencies( cath-refine-align   built_bioplib_libraries )
add_dependencies( cath-score-align    built_bioplib_libraries )
add_dependencies( cath-ssap           built_bioplib_libraries )
add_dependencies( cath-superpose      built_bioplib_libraries )

#####################################################################
#####################################################################

if ( DEFINED ENV{BUILD_EXTRA_CATH_TOOLS} )
	
	add_executable(
		snap-judgement
			${NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT}
	)

	add_executable(
		check-pdb
			${NORMSOURCES_EXECUTABLES_CATH_CHECK_PDB}
	)

#	set(
#		DEVELTEST_ALIGNMENT_SOURCES
#			alignment/alignment.cpp
#			alignment/alignment_row.cpp
#			alignment/residue_score/alignment_residue_scores.cpp
#			exception/invalid_argument_exception.cpp
#	)

#	add_executable(
#		develtest
#			executables/develtest/develtest.cpp
#
#			# superposition/superposition_test.cpp // requires complex spider of stuff
#
#			# alignment/alignment_test.cpp
#			# alignment/io/align_scaffold.cpp
#			# alignment/pair_alignment.cpp
#			# test/alignment_fixture.cpp
#			# ${DEVELTEST_ALIGNMENT_SOURCES}
#
#			alignment/tools/alignment_breaks_test.cpp
#			alignment/tools/alignment_breaks.cpp
#			alignment/io/align_scaffold.cpp
#			${DEVELTEST_ALIGNMENT_SOURCES}
#
#			# alignment/io/align_scaffold_test.cpp
#			# alignment/io/align_scaffold.cpp
#			# test/alignment_fixture.cpp
#			# ${DEVELTEST_ALIGNMENT_SOURCES}
#
#			# common/boost_addenda/range/adaptor/adaptor_test.cpp
#			# exception/invalid_argument_exception.cpp
#
#			# common/file/open_fstream.cpp
#			# exception/invalid_argument_exception.cpp
#			# exception/runtime_error_exception.cpp
#			# structure/geometry/orientation_covering.cpp
#			# structure/geometry/orientation_covering_test.cpp
#			# structure/geometry/quat_rot.cpp
#	)

	target_link_libraries( snap-judgement normsources_cath_superpose std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
	target_link_libraries( check-pdb                                 std_normsources_remain normsources_common normsources_exception biop gen pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )
#	target_link_libraries( develtest                                                                                                                  pthread ${Boost_LIBRARIES} ${RT_LIBRARY} ${GSL_LIB_SUFFIX} )

	add_dependencies( snap-judgement built_bioplib_libraries )
	add_dependencies( check-pdb      built_bioplib_libraries )

endif ( DEFINED ENV{BUILD_EXTRA_CATH_TOOLS} )

#####################################################################
#####################################################################

