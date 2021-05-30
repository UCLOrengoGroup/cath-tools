##### DON'T EDIT THIS FILE - IT'S AUTO-GENERATED #####

set(
	NORMSOURCES_CT_BIOCORE_CATH_BIOCORE
		ct_biocore/cath/biocore/chain_label.cpp
		ct_biocore/cath/biocore/residue_id.cpp
		ct_biocore/cath/biocore/residue_name.cpp
)

set(
	NORMSOURCES_CT_BIOCORE_CATH
		${NORMSOURCES_CT_BIOCORE_CATH_BIOCORE}
)

set(
	NORMSOURCES_CT_BIOCORE
		${NORMSOURCES_CT_BIOCORE_CATH}
)

set(
	NORMSOURCES_CT_CATH_ASSIGN_DOMAINS_CATH_CATH_ASSIGN_DOMAINS_OPTIONS
		ct_cath_assign_domains/cath/cath_assign_domains/options/cath_assign_domains_options.cpp
		ct_cath_assign_domains/cath/cath_assign_domains/options/cath_assign_domains_options_block.cpp
)

set(
	NORMSOURCES_CT_CATH_ASSIGN_DOMAINS_CATH_CATH_ASSIGN_DOMAINS
		${NORMSOURCES_CT_CATH_ASSIGN_DOMAINS_CATH_CATH_ASSIGN_DOMAINS_OPTIONS}
)

set(
	NORMSOURCES_CT_CATH_ASSIGN_DOMAINS_CATH
		${NORMSOURCES_CT_CATH_ASSIGN_DOMAINS_CATH_CATH_ASSIGN_DOMAINS}
)

set(
	NORMSOURCES_CT_CATH_ASSIGN_DOMAINS
		${NORMSOURCES_CT_CATH_ASSIGN_DOMAINS_CATH}
)

set(
	NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK
		ct_cath_cluster/cath/cath_cluster/options/options_block/cath_cluster_clustering_options_block.cpp
		ct_cath_cluster/cath/cath_cluster/options/options_block/cath_cluster_input_options_block.cpp
		ct_cath_cluster/cath/cath_cluster/options/options_block/cath_cluster_output_options_block.cpp
)

set(
	NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC
		ct_cath_cluster/cath/cath_cluster/options/spec/cath_cluster_clustering_spec.cpp
		ct_cath_cluster/cath/cath_cluster/options/spec/cath_cluster_input_spec.cpp
		ct_cath_cluster/cath/cath_cluster/options/spec/cath_cluster_output_spec.cpp
		ct_cath_cluster/cath/cath_cluster/options/spec/clustering_levels.cpp
)

set(
	NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS
		ct_cath_cluster/cath/cath_cluster/options/cath_cluster_options.cpp
		${NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK}
		${NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC}
)

set(
	NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER
		ct_cath_cluster/cath/cath_cluster/cath_clusterer.cpp
		${NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS}
)

set(
	NORMSOURCES_CT_CATH_CLUSTER_CATH
		${NORMSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER}
)

set(
	NORMSOURCES_CT_CATH_CLUSTER
		${NORMSOURCES_CT_CATH_CLUSTER_CATH}
)

set(
	NORMSOURCES_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN_OPTIONS
		ct_cath_refine_align/cath/cath_refine_align/options/cath_refine_align_options.cpp
)

set(
	NORMSOURCES_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN
		ct_cath_refine_align/cath/cath_refine_align/cath_align_refiner.cpp
		${NORMSOURCES_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN_OPTIONS}
)

set(
	NORMSOURCES_CT_CATH_REFINE_ALIGN_CATH
		${NORMSOURCES_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN}
)

set(
	NORMSOURCES_CT_CATH_REFINE_ALIGN
		${NORMSOURCES_CT_CATH_REFINE_ALIGN_CATH}
)

set(
	NORMSOURCES_CT_CATH_SCORE_ALIGN_CATH_CATH_SCORE_ALIGN_OPTIONS
		ct_cath_score_align/cath/cath_score_align/options/cath_score_align_options.cpp
)

set(
	NORMSOURCES_CT_CATH_SCORE_ALIGN_CATH_CATH_SCORE_ALIGN
		ct_cath_score_align/cath/cath_score_align/cath_align_scorer.cpp
		${NORMSOURCES_CT_CATH_SCORE_ALIGN_CATH_CATH_SCORE_ALIGN_OPTIONS}
)

set(
	NORMSOURCES_CT_CATH_SCORE_ALIGN_CATH
		${NORMSOURCES_CT_CATH_SCORE_ALIGN_CATH_CATH_SCORE_ALIGN}
)

set(
	NORMSOURCES_CT_CATH_SCORE_ALIGN
		${NORMSOURCES_CT_CATH_SCORE_ALIGN_CATH}
)

set(
	NORMSOURCES_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE_OPTIONS
		ct_cath_superpose/cath/cath_superpose/options/cath_superpose_options.cpp
)

set(
	NORMSOURCES_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE
		ct_cath_superpose/cath/cath_superpose/cath_superposer.cpp
		${NORMSOURCES_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE_OPTIONS}
)

set(
	NORMSOURCES_CT_CATH_SUPERPOSE_CATH
		${NORMSOURCES_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE}
)

set(
	NORMSOURCES_CT_CATH_SUPERPOSE
		${NORMSOURCES_CT_CATH_SUPERPOSE_CATH}
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT
		ct_chopping/cath/chopping/chopping_format/chopping_format.cpp
		ct_chopping/cath/chopping/chopping_format/domall_chopping_format.cpp
		ct_chopping/cath/chopping/chopping_format/jmol_selection_chopping_format.cpp
		ct_chopping/cath/chopping/chopping_format/scop_chopping_format.cpp
		ct_chopping/cath/chopping/chopping_format/sillitoe_chopping_format.cpp
		ct_chopping/cath/chopping/chopping_format/simple_chopping_format.cpp
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER
		ct_chopping/cath/chopping/chopping_io/region_io/region_reader/region_reader.cpp
		ct_chopping/cath/chopping/chopping_io/region_io/region_reader/std_region_reader.cpp
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO_REGION_WRITER
		ct_chopping/cath/chopping/chopping_io/region_io/region_writer/region_writer.cpp
		ct_chopping/cath/chopping/chopping_io/region_io/region_writer/std_region_writer.cpp
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER}
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO_REGION_WRITER}
		ct_chopping/cath/chopping/chopping_io/region_io/std_region_io_spec.cpp
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO_REGION_IO}
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_DOMAIN
		ct_chopping/cath/chopping/domain/domain.cpp
		ct_chopping/cath/chopping/domain/domain_definition.cpp
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_REGION
		ct_chopping/cath/chopping/region/region.cpp
		ct_chopping/cath/chopping/region/regions_limiter.cpp
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_RESIDUE_LOCATION
		ct_chopping/cath/chopping/residue_location/residue_locating.cpp
		ct_chopping/cath/chopping/residue_location/residue_location.cpp
)

set(
	NORMSOURCES_CT_CHOPPING_CATH_CHOPPING
		ct_chopping/cath/chopping/chopping.cpp
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT}
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_IO}
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_DOMAIN}
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_REGION}
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING_RESIDUE_LOCATION}
)

set(
	NORMSOURCES_CT_CHOPPING_CATH
		${NORMSOURCES_CT_CHOPPING_CATH_CHOPPING}
)

set(
	NORMSOURCES_CT_CHOPPING
		${NORMSOURCES_CT_CHOPPING_CATH}
)

set(
	NORMSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE
		ct_clustagglom/cath/clustagglom/file/dissimilarities_file.cpp
		ct_clustagglom/cath/clustagglom/file/names_file.cpp
)

set(
	NORMSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY
		ct_clustagglom/cath/clustagglom/hierarchy/hierarchy_group.cpp
		ct_clustagglom/cath/clustagglom/hierarchy/hierarchy_layer.cpp
		ct_clustagglom/cath/clustagglom/hierarchy/hierarchy_value.cpp
)

set(
	NORMSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM
		ct_clustagglom/cath/clustagglom/calc_complete_linkage_merge_list.cpp
		${NORMSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE}
		ct_clustagglom/cath/clustagglom/get_sorting_scores.cpp
		${NORMSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY}
		ct_clustagglom/cath/clustagglom/hierarchy.cpp
		ct_clustagglom/cath/clustagglom/link_dirn.cpp
		ct_clustagglom/cath/clustagglom/link_list.cpp
		ct_clustagglom/cath/clustagglom/links.cpp
		ct_clustagglom/cath/clustagglom/make_clusters_from_merges.cpp
		ct_clustagglom/cath/clustagglom/merge.cpp
)

set(
	NORMSOURCES_CT_CLUSTAGGLOM_CATH
		${NORMSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM}
)

set(
	NORMSOURCES_CT_CLUSTAGGLOM
		${NORMSOURCES_CT_CLUSTAGGLOM_CATH}
)

set(
	NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_DETAIL
		ct_cluster/cath/cluster/detail/mapping_job.cpp
)

set(
	NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_FILE
		ct_cluster/cath/cluster/file/cluster_membership_file.cpp
)

set(
	NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_MAP
		ct_cluster/cath/cluster/map/aggregate_map_results.cpp
		ct_cluster/cath/cluster/map/map_clusters.cpp
		ct_cluster/cath/cluster/map/map_results.cpp
		ct_cluster/cath/cluster/map/overlap_frac_distn.cpp
)

set(
	NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK
		ct_cluster/cath/cluster/options/options_block/clust_mapping_options_block.cpp
		ct_cluster/cath/cluster/options/options_block/clustmap_input_options_block.cpp
		ct_cluster/cath/cluster/options/options_block/clustmap_output_options_block.cpp
)

set(
	NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC
		ct_cluster/cath/cluster/options/spec/clustmap_input_spec.cpp
)

set(
	NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_OPTIONS
		${NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK}
		${NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC}
)

set(
	NORMSOURCES_CT_CLUSTER_CATH_CLUSTER
		ct_cluster/cath/cluster/cath_cluster_mapper.cpp
		ct_cluster/cath/cluster/cluster_domains.cpp
		ct_cluster/cath/cluster/clustmap_options.cpp
		${NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_DETAIL}
		ct_cluster/cath/cluster/domain_cluster_ids.cpp
		${NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_FILE}
		${NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_MAP}
		ct_cluster/cath/cluster/new_cluster_data.cpp
		ct_cluster/cath/cluster/old_cluster_data.cpp
		${NORMSOURCES_CT_CLUSTER_CATH_CLUSTER_OPTIONS}
)

set(
	NORMSOURCES_CT_CLUSTER_CATH
		${NORMSOURCES_CT_CLUSTER_CATH_CLUSTER}
)

set(
	NORMSOURCES_CT_CLUSTER
		${NORMSOURCES_CT_CLUSTER_CATH}
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON_ALGORITHM
		ct_common/cath/common/algorithm/random_split.cpp
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON_BATCH
		ct_common/cath/common/batch/batch_functions.cpp
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_GRAPH
		ct_common/cath/common/boost_addenda/graph/spanning_tree.cpp
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_LOG
		ct_common/cath/common/boost_addenda/log/log_to_ostream_guard.cpp
		ct_common/cath/common/boost_addenda/log/stringstream_log_sink.cpp
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA
		${NORMSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_GRAPH}
		${NORMSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_LOG}
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON_EXCEPTION
		ct_common/cath/common/exception/invalid_argument_exception.cpp
		ct_common/cath/common/exception/not_implemented_exception.cpp
		ct_common/cath/common/exception/out_of_range_exception.cpp
		ct_common/cath/common/exception/runtime_error_exception.cpp
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON_FILE
		ct_common/cath/common/file/find_file.cpp
		ct_common/cath/common/file/ofstream_list.cpp
		ct_common/cath/common/file/open_fstream.cpp
		ct_common/cath/common/file/path_or_istream.cpp
		ct_common/cath/common/file/temp_file.cpp
)

set(
	NORMSOURCES_CT_COMMON_CATH_COMMON
		${NORMSOURCES_CT_COMMON_CATH_COMMON_ALGORITHM}
		ct_common/cath/common/argc_argv_faker.cpp
		${NORMSOURCES_CT_COMMON_CATH_COMMON_BATCH}
		${NORMSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA}
		ct_common/cath/common/command_executer.cpp
		${NORMSOURCES_CT_COMMON_CATH_COMMON_EXCEPTION}
		${NORMSOURCES_CT_COMMON_CATH_COMMON_FILE}
		ct_common/cath/common/logger.cpp
		ct_common/cath/common/program_exception_wrapper.cpp
		ct_common/cath/common/test_or_exe_run_mode.cpp
)

set(
	NORMSOURCES_CT_COMMON_CATH
		${NORMSOURCES_CT_COMMON_CATH_COMMON}
)

set(
	NORMSOURCES_CT_COMMON
		${NORMSOURCES_CT_COMMON_CATH}
)

set(
	NORMSOURCES_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR
		ct_display_colour/cath/display_colour/display_colour.cpp
		ct_display_colour/cath/display_colour/display_colour_gradient.cpp
		ct_display_colour/cath/display_colour/display_colour_list.cpp
)

set(
	NORMSOURCES_CT_DISPLAY_COLOUR_CATH
		${NORMSOURCES_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR}
)

set(
	NORMSOURCES_CT_DISPLAY_COLOUR
		${NORMSOURCES_CT_DISPLAY_COLOUR_CATH}
)

set(
	NORMSOURCES_CT_EXTERNAL_INFO_CATH_EXTERNAL_INFO
		ct_external_info/cath/external_info/cath_tools_cmake_dirs.cpp
		ct_external_info/cath/external_info/cath_tools_git_version.cpp
)

set(
	NORMSOURCES_CT_EXTERNAL_INFO_CATH
		${NORMSOURCES_CT_EXTERNAL_INFO_CATH_EXTERNAL_INFO}
)

set(
	NORMSOURCES_CT_EXTERNAL_INFO
		${NORMSOURCES_CT_EXTERNAL_INFO_CATH}
)

set(
	NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS
		ct_options/cath/options/executable/cath_check_pdb_options/cath_check_pdb_options.cpp
)

set(
	NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_CATH_EXTRACT_PDB_OPTIONS
		ct_options/cath/options/executable/cath_extract_pdb_options/cath_extract_pdb_options.cpp
)

set(
	NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE
		${NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_CATH_CHECK_PDB_OPTIONS}
		${NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_CATH_EXTRACT_PDB_OPTIONS}
		ct_options/cath/options/executable/env_var_option_name_handler.cpp
		ct_options/cath/options/executable/executable_options.cpp
)

set(
	NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK
		ct_options/cath/options/options_block/check_pdb_options_block.cpp
		ct_options/cath/options/options_block/detail_help_options_block.cpp
		ct_options/cath/options/options_block/extract_pdb_options_block.cpp
		ct_options/cath/options/options_block/ids_options_block.cpp
		ct_options/cath/options/options_block/misc_help_version_options_block.cpp
		ct_options/cath/options/options_block/options_block.cpp
		ct_options/cath/options/options_block/options_block_tester.cpp
		ct_options/cath/options/options_block/pdb_input_options_block.cpp
		ct_options/cath/options/options_block/pdb_input_spec.cpp
		ct_options/cath/options/options_block/string_options_block.cpp
		ct_options/cath/options/options_block/superposition_input_options_block.cpp
)

set(
	NORMSOURCES_CT_OPTIONS_CATH_OPTIONS
		${NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE}
		${NORMSOURCES_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK}
)

set(
	NORMSOURCES_CT_OPTIONS_CATH
		${NORMSOURCES_CT_OPTIONS_CATH_OPTIONS}
)

set(
	NORMSOURCES_CT_OPTIONS
		${NORMSOURCES_CT_OPTIONS_CATH}
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO
		ct_resolve_hits/cath/resolve_hits/algo/discont_hits_index_by_start.cpp
		ct_resolve_hits/cath/resolve_hits/algo/masked_bests_cacher.cpp
		ct_resolve_hits/cath/resolve_hits/algo/scored_arch_proxy.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE
		ct_resolve_hits/cath/resolve_hits/file/alnd_rgn.cpp
		ct_resolve_hits/cath/resolve_hits/file/cath_id_score_category.cpp
		ct_resolve_hits/cath/resolve_hits/file/hits_input_format_tag.cpp
		ct_resolve_hits/cath/resolve_hits/file/parse_domain_hits_table.cpp
		ct_resolve_hits/cath/resolve_hits/file/parse_hmmer_out.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT
		ct_resolve_hits/cath/resolve_hits/html_output/html_segment.cpp
		ct_resolve_hits/cath/resolve_hits/html_output/resolve_hits_html_outputter.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK
		ct_resolve_hits/cath/resolve_hits/options/options_block/crh_filter_options_block.cpp
		ct_resolve_hits/cath/resolve_hits/options/options_block/crh_html_options_block.cpp
		ct_resolve_hits/cath/resolve_hits/options/options_block/crh_input_options_block.cpp
		ct_resolve_hits/cath/resolve_hits/options/options_block/crh_output_options_block.cpp
		ct_resolve_hits/cath/resolve_hits/options/options_block/crh_score_options_block.cpp
		ct_resolve_hits/cath/resolve_hits/options/options_block/crh_segment_options_block.cpp
		ct_resolve_hits/cath/resolve_hits/options/options_block/crh_single_output_options_block.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_filter_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_html_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_input_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_output_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_score_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_segment_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_single_output_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_spec.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/hit_boundary_output.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS
		ct_resolve_hits/cath/resolve_hits/options/crh_options.cpp
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_OPTIONS_BLOCK}
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC}
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/gather_hits_processor.cpp
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor_list.cpp
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/summarise_hits_processor.cpp
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/write_html_hits_processor.cpp
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/write_json_hits_processor.cpp
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/write_results_hits_processor.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR}
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/read_and_process_mgr.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_RESOLVE
		ct_resolve_hits/cath/resolve_hits/resolve/hit_resolver.cpp
		ct_resolve_hits/cath/resolve_hits/resolve/naive_greedy_hit_resolver.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TRIM
		ct_resolve_hits/cath/resolve_hits/trim/seq_seg_boundary_fns.cpp
		ct_resolve_hits/cath/resolve_hits/trim/trim_spec.cpp
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO}
		ct_resolve_hits/cath/resolve_hits/calc_hit.cpp
		ct_resolve_hits/cath/resolve_hits/calc_hit_list.cpp
		ct_resolve_hits/cath/resolve_hits/cath_hit_resolver.cpp
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE}
		ct_resolve_hits/cath/resolve_hits/full_hit.cpp
		ct_resolve_hits/cath/resolve_hits/full_hit_fns.cpp
		ct_resolve_hits/cath/resolve_hits/full_hit_list.cpp
		ct_resolve_hits/cath/resolve_hits/full_hit_list_fns.cpp
		ct_resolve_hits/cath/resolve_hits/hit_arch.cpp
		ct_resolve_hits/cath/resolve_hits/hit_extras.cpp
		ct_resolve_hits/cath/resolve_hits/hit_score_type.cpp
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT}
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS}
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS}
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_RESOLVE}
		ct_resolve_hits/cath/resolve_hits/scored_hit_arch.cpp
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TRIM}
)

set(
	NORMSOURCES_CT_RESOLVE_HITS_CATH
		${NORMSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS}
)

set(
	NORMSOURCES_CT_RESOLVE_HITS
		${NORMSOURCES_CT_RESOLVE_HITS_CATH}
)

set(
	NORMSOURCES_CT_SEQ_CATH_SEQ
		ct_seq/cath/seq/seq_arrow.cpp
		ct_seq/cath/seq/seq_seg.cpp
		ct_seq/cath/seq/seq_seg_run.cpp
)

set(
	NORMSOURCES_CT_SEQ_CATH
		${NORMSOURCES_CT_SEQ_CATH_SEQ}
)

set(
	NORMSOURCES_CT_SEQ
		${NORMSOURCES_CT_SEQ_CATH}
)

set(
	NORMSOURCES_CT_TEST_CATH_TEST_PREDICATE_DETAIL
		ct_test/cath/test/predicate/detail/strings_equal.cpp
)

set(
	NORMSOURCES_CT_TEST_CATH_TEST_PREDICATE
		ct_test/cath/test/predicate/bootstrap_mode.cpp
		${NORMSOURCES_CT_TEST_CATH_TEST_PREDICATE_DETAIL}
		ct_test/cath/test/predicate/files_equal.cpp
		ct_test/cath/test/predicate/istream_and_file_equal.cpp
		ct_test/cath/test/predicate/istreams_equal.cpp
		ct_test/cath/test/predicate/string_matches_file.cpp
)

set(
	NORMSOURCES_CT_TEST_CATH_TEST
		ct_test/cath/test/global_test_constants.cpp
		${NORMSOURCES_CT_TEST_CATH_TEST_PREDICATE}
)

set(
	NORMSOURCES_CT_TEST_CATH
		${NORMSOURCES_CT_TEST_CATH_TEST}
)

set(
	NORMSOURCES_CT_TEST
		${NORMSOURCES_CT_TEST_CATH}
)

set(
	NORMSOURCES_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER
		ct_uni/cath/acquirer/alignment_acquirer/align_refining.cpp
		ct_uni/cath/acquirer/alignment_acquirer/alignment_acquirer.cpp
		ct_uni/cath/acquirer/alignment_acquirer/cora_aln_file_alignment_acquirer.cpp
		ct_uni/cath/acquirer/alignment_acquirer/do_the_ssaps_alignment_acquirer.cpp
		ct_uni/cath/acquirer/alignment_acquirer/fasta_aln_file_alignment_acquirer.cpp
		ct_uni/cath/acquirer/alignment_acquirer/post_refine_alignment_acquirer.cpp
		ct_uni/cath/acquirer/alignment_acquirer/residue_name_alignment_acquirer.cpp
		ct_uni/cath/acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer.cpp
		ct_uni/cath/acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER
		ct_uni/cath/acquirer/pdbs_acquirer/domain_defn_pdbs_acquirer.cpp
		ct_uni/cath/acquirer/pdbs_acquirer/file_list_pdbs_acquirer.cpp
		ct_uni/cath/acquirer/pdbs_acquirer/istream_pdbs_acquirer.cpp
		ct_uni/cath/acquirer/pdbs_acquirer/pdbs_acquirer.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ACQUIRER_SELECTION_POLICY_ACQUIRER
		ct_uni/cath/acquirer/selection_policy_acquirer/selection_policy_acquirer.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER
		ct_uni/cath/acquirer/superposition_acquirer/align_based_superposition_acquirer.cpp
		ct_uni/cath/acquirer/superposition_acquirer/superposition_acquirer.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ACQUIRER
		${NORMSOURCES_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER}
		${NORMSOURCES_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER}
		${NORMSOURCES_CT_UNI_CATH_ACQUIRER_SELECTION_POLICY_ACQUIRER}
		${NORMSOURCES_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER}
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY
		ct_uni/cath/alignment/common_atom_selection_policy/common_atom_select_ca_policy.cpp
		ct_uni/cath/alignment/common_atom_selection_policy/common_atom_select_cb_policy.cpp
		ct_uni/cath/alignment/common_atom_selection_policy/common_atom_selection_policy.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY
		ct_uni/cath/alignment/common_residue_selection_policy/common_residue_score_based_selection_policy.cpp
		ct_uni/cath/alignment/common_residue_selection_policy/common_residue_select_all_policy.cpp
		ct_uni/cath/alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.cpp
		ct_uni/cath/alignment/common_residue_selection_policy/common_residue_select_min_score_policy.cpp
		ct_uni/cath/alignment/common_residue_selection_policy/common_residue_selection_policy.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DETAIL
		ct_uni/cath/alignment/detail/multi_align_builder.cpp
		ct_uni/cath/alignment/detail/multi_align_group.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER
		ct_uni/cath/alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.cpp
		ct_uni/cath/alignment/dyn_prog_align/detail/matrix_plotter/matrix_plotter.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER
		ct_uni/cath/alignment/dyn_prog_align/detail/string_aligner/benchmark_dyn_prog_string_aligner.cpp
		ct_uni/cath/alignment/dyn_prog_align/detail/string_aligner/gen_dyn_prog_string_aligner.cpp
		ct_uni/cath/alignment/dyn_prog_align/detail/string_aligner/string_aligner.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_MATRIX_PLOTTER}
		ct_uni/cath/alignment/dyn_prog_align/detail/path_step.cpp
		ct_uni/cath/alignment/dyn_prog_align/detail/return_path_matrix.cpp
		ct_uni/cath/alignment/dyn_prog_align/detail/score_accumulation_matrix.cpp
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER}
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.cpp
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_score_source/entry_querier_dyn_prog_score_source.cpp
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source.cpp
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source.cpp
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source.cpp
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL}
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_aligner.cpp
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE}
		ct_uni/cath/alignment/dyn_prog_align/ssap_code_dyn_prog_aligner.cpp
		ct_uni/cath/alignment/dyn_prog_align/std_dyn_prog_aligner.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_GAP
		ct_uni/cath/alignment/gap/alignment_gap.cpp
		ct_uni/cath/alignment/gap/gap_penalty.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_IO_OUTPUTTER
		ct_uni/cath/alignment/io/outputter/horiz_align_outputter.cpp
		ct_uni/cath/alignment/io/outputter/html_align_outputter.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_IO
		ct_uni/cath/alignment/io/align_scaffold.cpp
		ct_uni/cath/alignment/io/alignment_io.cpp
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_IO_OUTPUTTER}
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK
		ct_uni/cath/alignment/options_block/alignment_input_options_block.cpp
		ct_uni/cath/alignment/options_block/alignment_input_spec.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER_DETAIL
		ct_uni/cath/alignment/refiner/detail/alignment_split.cpp
		ct_uni/cath/alignment/refiner/detail/alignment_split_list.cpp
		ct_uni/cath/alignment/refiner/detail/alignment_split_mapping.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER
		ct_uni/cath/alignment/refiner/alignment_refiner.cpp
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER_DETAIL}
		ct_uni/cath/alignment/refiner/indexed_refiner.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL
		ct_uni/cath/alignment/residue_name_align/detail/residue_name_align_map.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL}
		ct_uni/cath/alignment/residue_name_align/residue_name_aligner.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_SCORE
		ct_uni/cath/alignment/residue_score/alignment_residue_scores.cpp
		ct_uni/cath/alignment/residue_score/residue_scorer.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT_TOOLS
		ct_uni/cath/alignment/tools/alignment_breaks.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_ALIGNMENT
		ct_uni/cath/alignment/alignment.cpp
		ct_uni/cath/alignment/alignment_action.cpp
		ct_uni/cath/alignment/alignment_context.cpp
		ct_uni/cath/alignment/alignment_coord_extractor.cpp
		ct_uni/cath/alignment/alignment_row.cpp
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DETAIL}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_GAP}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_IO}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_OPTIONS_BLOCK}
		ct_uni/cath/alignment/pair_alignment.cpp
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_SCORE}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT_TOOLS}
)

set(
	NORMSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC
		ct_uni/cath/display/display_colour_spec/broad_display_colour_spec.cpp
		ct_uni/cath/display/display_colour_spec/display_colour_spec.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOURER_DETAIL
		ct_uni/cath/display/display_colourer/detail/score_colour_handler.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOURER
		ct_uni/cath/display/display_colourer/alignment_free_display_colourer.cpp
		${NORMSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOURER_DETAIL}
		ct_uni/cath/display/display_colourer/display_colourer.cpp
		ct_uni/cath/display/display_colourer/display_colourer_alignment.cpp
		ct_uni/cath/display/display_colourer/display_colourer_consecutive.cpp
		ct_uni/cath/display/display_colourer/display_colourer_score.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_DISPLAY_OPTIONS
		ct_uni/cath/display/options/display_options_block.cpp
		ct_uni/cath/display/options/display_spec.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL
		ct_uni/cath/display/viewer/pymol/pymol_tools.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_DISPLAY_VIEWER
		ct_uni/cath/display/viewer/chimera_viewer.cpp
		ct_uni/cath/display/viewer/jmol_viewer.cpp
		${NORMSOURCES_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL}
		ct_uni/cath/display/viewer/pymol_viewer.cpp
		ct_uni/cath/display/viewer/rasmol_style_viewer.cpp
		ct_uni/cath/display/viewer/rasmol_viewer.cpp
		ct_uni/cath/display/viewer/viewer.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_DISPLAY
		${NORMSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC}
		${NORMSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOURER}
		${NORMSOURCES_CT_UNI_CATH_DISPLAY_OPTIONS}
		${NORMSOURCES_CT_UNI_CATH_DISPLAY_VIEWER}
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_DOMAIN_DEFINITION_LIST
		ct_uni/cath/file/domain_definition_list/domain_definition_list.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_DSSP_WOLF
		ct_uni/cath/file/dssp_wolf/dssp_file.cpp
		ct_uni/cath/file/dssp_wolf/dssp_file_io.cpp
		ct_uni/cath/file/dssp_wolf/tally_residue_ids.cpp
		ct_uni/cath/file/dssp_wolf/wolf_file.cpp
		ct_uni/cath/file/dssp_wolf/wolf_file_io.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_HMMER_SCORES_FILE
		ct_uni/cath/file/hmmer_scores_file/hmmer_scores_entry.cpp
		ct_uni/cath/file/hmmer_scores_file/hmmer_scores_file.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_NAME_SET
		ct_uni/cath/file/name_set/name_set.cpp
		ct_uni/cath/file/name_set/name_set_list.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_OPTIONS
		ct_uni/cath/file/options/data_dirs_options_block.cpp
		ct_uni/cath/file/options/data_dirs_spec.cpp
		ct_uni/cath/file/options/data_option.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_PDB
		ct_uni/cath/file/pdb/coarse_element_type.cpp
		ct_uni/cath/file/pdb/dssp_skip_policy.cpp
		ct_uni/cath/file/pdb/pdb.cpp
		ct_uni/cath/file/pdb/pdb_atom.cpp
		ct_uni/cath/file/pdb/pdb_atom_parse_status.cpp
		ct_uni/cath/file/pdb/pdb_list.cpp
		ct_uni/cath/file/pdb/pdb_record.cpp
		ct_uni/cath/file/pdb/pdb_residue.cpp
		ct_uni/cath/file/pdb/proximity_calculator.cpp
		ct_uni/cath/file/pdb/read_domain_def_from_pdb.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_PRC_SCORES_FILE
		ct_uni/cath/file/prc_scores_file/prc_scores_entry.cpp
		ct_uni/cath/file/prc_scores_file/prc_scores_file.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_SEC
		ct_uni/cath/file/sec/sec_file.cpp
		ct_uni/cath/file/sec/sec_file_io.cpp
		ct_uni/cath/file/sec/sec_file_record.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE_SSAP_SCORES_FILE
		ct_uni/cath/file/ssap_scores_file/ssap_scores_entry.cpp
		ct_uni/cath/file/ssap_scores_file/ssap_scores_file.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_FILE
		ct_uni/cath/file/data_file.cpp
		${NORMSOURCES_CT_UNI_CATH_FILE_DOMAIN_DEFINITION_LIST}
		${NORMSOURCES_CT_UNI_CATH_FILE_DSSP_WOLF}
		${NORMSOURCES_CT_UNI_CATH_FILE_HMMER_SCORES_FILE}
		${NORMSOURCES_CT_UNI_CATH_FILE_NAME_SET}
		${NORMSOURCES_CT_UNI_CATH_FILE_OPTIONS}
		${NORMSOURCES_CT_UNI_CATH_FILE_PDB}
		${NORMSOURCES_CT_UNI_CATH_FILE_PRC_SCORES_FILE}
		${NORMSOURCES_CT_UNI_CATH_FILE_SEC}
		${NORMSOURCES_CT_UNI_CATH_FILE_SSAP_SCORES_FILE}
		ct_uni/cath/file/strucs_context.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER
		ct_uni/cath/outputter/alignment_outputter/alignment_outputter.cpp
		ct_uni/cath/outputter/alignment_outputter/alignment_outputter_list.cpp
		ct_uni/cath/outputter/alignment_outputter/cath_aln_ostream_alignment_outputter.cpp
		ct_uni/cath/outputter/alignment_outputter/fasta_ostream_alignment_outputter.cpp
		ct_uni/cath/outputter/alignment_outputter/file_alignment_outputter.cpp
		ct_uni/cath/outputter/alignment_outputter/html_ostream_alignment_outputter.cpp
		ct_uni/cath/outputter/alignment_outputter/ssap_ostream_alignment_outputter.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS
		ct_uni/cath/outputter/alignment_outputter_options/alignment_output_options_block.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS
		ct_uni/cath/outputter/superposition_output_options/superposition_output_options_block.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUTTER
		ct_uni/cath/outputter/superposition_outputter/json_file_superposition_outputter.cpp
		ct_uni/cath/outputter/superposition_outputter/ostream_superposition_outputter.cpp
		ct_uni/cath/outputter/superposition_outputter/pdb_file_superposition_outputter.cpp
		ct_uni/cath/outputter/superposition_outputter/pdb_files_superposition_outputter.cpp
		ct_uni/cath/outputter/superposition_outputter/pymol_file_superposition_outputter.cpp
		ct_uni/cath/outputter/superposition_outputter/pymol_view_superposition_outputter.cpp
		ct_uni/cath/outputter/superposition_outputter/superposition_outputter.cpp
		ct_uni/cath/outputter/superposition_outputter/superposition_outputter_list.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_OUTPUTTER
		${NORMSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER}
		${NORMSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_OPTIONS}
		ct_uni/cath/outputter/run_pymol.cpp
		${NORMSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS}
		${NORMSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUTTER}
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY
		ct_uni/cath/scan/detail/check_scan/test_only/alignment_scan_comparison.cpp
		ct_uni/cath/scan/detail/check_scan/test_only/check_scan_on_final_alignment.cpp
		ct_uni/cath/scan/detail/check_scan/test_only/quad_and_rep_criteria_result.cpp
		ct_uni/cath/scan/detail/check_scan/test_only/quad_criteria_result.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN
		${NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY}
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR
		ct_uni/cath/scan/detail/res_pair/multi_struc_res_rep_pair.cpp
		ct_uni/cath/scan/detail/res_pair/res_pair_core.cpp
		ct_uni/cath/scan/detail/res_pair/single_struc_res_pair.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_DIRN
		ct_uni/cath/scan/detail/res_pair_dirn/res_pair_dirn.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_STRIDE
		ct_uni/cath/scan/detail/stride/rep_strider.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL
		${NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN}
		${NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR}
		${NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_DIRN}
		${NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL_STRIDE}
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_SCAN_ACTION
		ct_uni/cath/scan/scan_action/populate_matrix_scan_action.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_SCAN_TOOLS
		ct_uni/cath/scan/scan_tools/all_vs_all.cpp
		ct_uni/cath/scan/scan_tools/load_and_scan.cpp
		ct_uni/cath/scan/scan_tools/load_and_scan_metrics.cpp
		ct_uni/cath/scan/scan_tools/scan_metrics.cpp
		ct_uni/cath/scan/scan_tools/scan_type.cpp
		ct_uni/cath/scan/scan_tools/single_pair.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN_SPATIAL_INDEX
		ct_uni/cath/scan/spatial_index/spatial_index.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCAN
		${NORMSOURCES_CT_UNI_CATH_SCAN_DETAIL}
		ct_uni/cath/scan/quad_criteria.cpp
		ct_uni/cath/scan/res_pair_index_dirn_criterion.cpp
		${NORMSOURCES_CT_UNI_CATH_SCAN_SCAN_ACTION}
		ct_uni/cath/scan/scan_stride.cpp
		${NORMSOURCES_CT_UNI_CATH_SCAN_SCAN_TOOLS}
		${NORMSOURCES_CT_UNI_CATH_SCAN_SPATIAL_INDEX}
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_DETAIL
		ct_uni/cath/score/aligned_pair_score/detail/score_common_coord_handler.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE
		ct_uni/cath/score/aligned_pair_score/ssap_score/ssap_score_accuracy.cpp
		ct_uni/cath/score/aligned_pair_score/ssap_score/ssap_score_post_processing.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX
		ct_uni/cath/score/aligned_pair_score/substitution_matrix/blosum62_substitution_matrix.cpp
		ct_uni/cath/score/aligned_pair_score/substitution_matrix/identity_substitution_matrix.cpp
		ct_uni/cath/score/aligned_pair_score/substitution_matrix/match_substitution_matrix.cpp
		ct_uni/cath/score/aligned_pair_score/substitution_matrix/substitution_matrix.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE
		ct_uni/cath/score/aligned_pair_score/aligned_pair_score.cpp
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_DETAIL}
		ct_uni/cath/score/aligned_pair_score/drmsd_score.cpp
		ct_uni/cath/score/aligned_pair_score/gsas_score.cpp
		ct_uni/cath/score/aligned_pair_score/lddt_score.cpp
		ct_uni/cath/score/aligned_pair_score/length_score.cpp
		ct_uni/cath/score/aligned_pair_score/mi_score.cpp
		ct_uni/cath/score/aligned_pair_score/overlap_score.cpp
		ct_uni/cath/score/aligned_pair_score/pseudo_string_score.cpp
		ct_uni/cath/score/aligned_pair_score/rmsd_score.cpp
		ct_uni/cath/score/aligned_pair_score/sas_score.cpp
		ct_uni/cath/score/aligned_pair_score/sequence_similarity_score.cpp
		ct_uni/cath/score/aligned_pair_score/si_score.cpp
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE}
		ct_uni/cath/score/aligned_pair_score/ssap_score.cpp
		ct_uni/cath/score/aligned_pair_score/structal_score.cpp
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX}
		ct_uni/cath/score/aligned_pair_score/tm_score.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL
		ct_uni/cath/score/aligned_pair_score_list/detail/aligned_pair_score_list_append.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER
		ct_uni/cath/score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_READER
		ct_uni/cath/score/aligned_pair_score_list/score_value_list_reader/score_value_reader.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST
		ct_uni/cath/score/aligned_pair_score_list/aligned_pair_score_list.cpp
		ct_uni/cath/score/aligned_pair_score_list/aligned_pair_score_list_factory.cpp
		ct_uni/cath/score/aligned_pair_score_list/aligned_pair_score_value_list.cpp
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL}
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_OUTPUTTER}
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_SCORE_VALUE_LIST_READER}
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_DETAIL
		ct_uni/cath/score/detail/score_name_helper.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS
		ct_uni/cath/score/homcheck_tools/ssap_and_prc.cpp
		ct_uni/cath/score/homcheck_tools/ssaps_and_prcs_of_query.cpp
		ct_uni/cath/score/homcheck_tools/superfamily_of_domain.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_LENGTH_GETTER
		ct_uni/cath/score/length_getter/geometric_mean_length_getter.cpp
		ct_uni/cath/score/length_getter/length_getter.cpp
		ct_uni/cath/score/length_getter/length_getter_make_clone.cpp
		ct_uni/cath/score/length_getter/length_of_first_getter.cpp
		ct_uni/cath/score/length_getter/length_of_longer_getter.cpp
		ct_uni/cath/score/length_getter/length_of_second_getter.cpp
		ct_uni/cath/score/length_getter/length_of_shorter_getter.cpp
		ct_uni/cath/score/length_getter/mean_length_getter.cpp
		ct_uni/cath/score/length_getter/num_aligned_length_getter.cpp
		ct_uni/cath/score/length_getter/protein_only_length_getter.cpp
		ct_uni/cath/score/length_getter/sym_protein_only_length_getter.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_PAIR_SCATTER_PLOTTER
		ct_uni/cath/score/pair_scatter_plotter/pair_scatter_plotter.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_DETAIL
		ct_uni/cath/score/score_classification/detail/score_classn_value_list_name_less.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE
		ct_uni/cath/score/score_classification/label_pair_is_positive/label_pair_is_positive.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION
		${NORMSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_DETAIL}
		${NORMSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE}
		ct_uni/cath/score/score_classification/rbf_model.cpp
		ct_uni/cath/score/score_classification/score_classn_value.cpp
		ct_uni/cath/score/score_classification/score_classn_value_better_value.cpp
		ct_uni/cath/score/score_classification/score_classn_value_list.cpp
		ct_uni/cath/score/score_classification/score_classn_value_results_set.cpp
		ct_uni/cath/score/score_classification/value_list_scaling.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER
		ct_uni/cath/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter.cpp
		ct_uni/cath/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_spec.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG
		ct_uni/cath/score/true_pos_false_neg/classn_outcome.cpp
		ct_uni/cath/score/true_pos_false_neg/classn_stat.cpp
		ct_uni/cath/score/true_pos_false_neg/classn_stat_pair_series.cpp
		ct_uni/cath/score/true_pos_false_neg/classn_stat_pair_series_list.cpp
		${NORMSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER}
		ct_uni/cath/score/true_pos_false_neg/named_true_false_pos_neg_list.cpp
		ct_uni/cath/score/true_pos_false_neg/named_true_false_pos_neg_list_list.cpp
		ct_uni/cath/score/true_pos_false_neg/true_false_pos_neg.cpp
		ct_uni/cath/score/true_pos_false_neg/true_false_pos_neg_list.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SCORE
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE}
		${NORMSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST}
		${NORMSOURCES_CT_UNI_CATH_SCORE_DETAIL}
		${NORMSOURCES_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS}
		${NORMSOURCES_CT_UNI_CATH_SCORE_LENGTH_GETTER}
		${NORMSOURCES_CT_UNI_CATH_SCORE_PAIR_SCATTER_PLOTTER}
		${NORMSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION}
		${NORMSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG}
)

set(
	NORMSOURCES_CT_UNI_CATH_SSAP_OPTIONS
		ct_uni/cath/ssap/options/cath_ssap_options.cpp
		ct_uni/cath/ssap/options/old_ssap_options_block.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SSAP
		ct_uni/cath/ssap/distance_score_formula.cpp
		${NORMSOURCES_CT_UNI_CATH_SSAP_OPTIONS}
		ct_uni/cath/ssap/selected_pair.cpp
		ct_uni/cath/ssap/ssap.cpp
		ct_uni/cath/ssap/ssap_scores.cpp
		ct_uni/cath/ssap/windowed_matrix.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_ACCESSIBILITY_CALC
		ct_uni/cath/structure/accessibility_calc/dssp_accessibility.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_ENTRY_QUERIER
		ct_uni/cath/structure/entry_querier/entry_querier.cpp
		ct_uni/cath/structure/entry_querier/residue_querier.cpp
		ct_uni/cath/structure/entry_querier/sec_struc_querier.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_GEOMETRY_DETAIL
		ct_uni/cath/structure/geometry/detail/cross_covariance_matrix.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_GEOMETRY
		ct_uni/cath/structure/geometry/angle.cpp
		ct_uni/cath/structure/geometry/coord.cpp
		ct_uni/cath/structure/geometry/coord_list.cpp
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_GEOMETRY_DETAIL}
		ct_uni/cath/structure/geometry/orient.cpp
		ct_uni/cath/structure/geometry/pca.cpp
		ct_uni/cath/structure/geometry/restrict_to_single_linkage_extension.cpp
		ct_uni/cath/structure/geometry/rotation.cpp
		ct_uni/cath/structure/geometry/superpose_fit.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_LOADER
		ct_uni/cath/structure/protein/protein_loader/protein_list_loader.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET
		ct_uni/cath/structure/protein/protein_source_file_set/protein_file_combn.cpp
		ct_uni/cath/structure/protein/protein_source_file_set/protein_from_pdb.cpp
		ct_uni/cath/structure/protein/protein_source_file_set/protein_from_pdb_and_calc.cpp
		ct_uni/cath/structure/protein/protein_source_file_set/protein_from_pdb_and_dssp_and_calc.cpp
		ct_uni/cath/structure/protein/protein_source_file_set/protein_from_pdb_dssp_and_sec.cpp
		ct_uni/cath/structure/protein/protein_source_file_set/protein_from_wolf_and_sec.cpp
		ct_uni/cath/structure/protein/protein_source_file_set/protein_source_file_set.cpp
		ct_uni/cath/structure/protein/protein_source_file_set/restrict_protein_source_file_set.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN
		ct_uni/cath/structure/protein/amino_acid.cpp
		ct_uni/cath/structure/protein/protein.cpp
		ct_uni/cath/structure/protein/protein_io.cpp
		ct_uni/cath/structure/protein/protein_list.cpp
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_LOADER}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET}
		ct_uni/cath/structure/protein/residue.cpp
		ct_uni/cath/structure/protein/sec_struc.cpp
		ct_uni/cath/structure/protein/sec_struc_planar_angles.cpp
		ct_uni/cath/structure/protein/sec_struc_type.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP
		ct_uni/cath/structure/sec_struc_calc/dssp/bifur_hbond_list.cpp
		ct_uni/cath/structure/sec_struc_calc/dssp/dssp_hbond_calc.cpp
		ct_uni/cath/structure/sec_struc_calc/dssp/dssp_ss_calc.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_SEC
		ct_uni/cath/structure/sec_struc_calc/sec/sec_calc.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_SEC}
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_DETAIL
		ct_uni/cath/structure/view_cache/filter/detail/filter_vs_full_score_less.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_DETAIL}
		ct_uni/cath/structure/view_cache/filter/filter_vs_full_score.cpp
		ct_uni/cath/structure/view_cache/filter/filter_vs_full_score_list.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL
		ct_uni/cath/structure/view_cache/index/detail/vcie_match_criteria.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL}
		ct_uni/cath/structure/view_cache/index/quad_find_action.cpp
		ct_uni/cath/structure/view_cache/index/quad_find_action_check.cpp
		ct_uni/cath/structure/view_cache/index/view_cache_index.cpp
		ct_uni/cath/structure/view_cache/index/view_cache_index_entry.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX}
		ct_uni/cath/structure/view_cache/view_cache.cpp
		ct_uni/cath/structure/view_cache/view_cache_list.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_STRUCTURE
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_ACCESSIBILITY_CALC}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_ENTRY_QUERIER}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_GEOMETRY}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE}
)

set(
	NORMSOURCES_CT_UNI_CATH_SUPERPOSITION_IO
		ct_uni/cath/superposition/io/superposition_io.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SUPERPOSITION_OPTIONS
		ct_uni/cath/superposition/options/align_regions_options_block.cpp
		ct_uni/cath/superposition/options/superposition_content_options_block.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH_SUPERPOSITION
		${NORMSOURCES_CT_UNI_CATH_SUPERPOSITION_IO}
		${NORMSOURCES_CT_UNI_CATH_SUPERPOSITION_OPTIONS}
		ct_uni/cath/superposition/superpose_orient.cpp
		ct_uni/cath/superposition/superposition.cpp
		ct_uni/cath/superposition/superposition_content_spec.cpp
		ct_uni/cath/superposition/superposition_context.cpp
		ct_uni/cath/superposition/supn_regions_context.cpp
)

set(
	NORMSOURCES_CT_UNI_CATH
		${NORMSOURCES_CT_UNI_CATH_ACQUIRER}
		${NORMSOURCES_CT_UNI_CATH_ALIGNMENT}
		${NORMSOURCES_CT_UNI_CATH_DISPLAY}
		${NORMSOURCES_CT_UNI_CATH_FILE}
		${NORMSOURCES_CT_UNI_CATH_OUTPUTTER}
		${NORMSOURCES_CT_UNI_CATH_SCAN}
		${NORMSOURCES_CT_UNI_CATH_SCORE}
		${NORMSOURCES_CT_UNI_CATH_SSAP}
		${NORMSOURCES_CT_UNI_CATH_STRUCTURE}
		${NORMSOURCES_CT_UNI_CATH_SUPERPOSITION}
)

set(
	NORMSOURCES_CT_UNI
		${NORMSOURCES_CT_UNI_CATH}
)

set(
	NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS
		executables/cath_assign_domains/cath_assign_domains.cpp
)

set(
	NORMSOURCES_EXECUTABLES_CATH_CLUSTER
		executables/cath_cluster/cath_cluster.cpp
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
	NORMSOURCES_EXECUTABLES_CHECK_PDB
		executables/check_pdb/cath_check_pdb.cpp
)

set(
	NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT
		executables/snap_judgement/snap_judgement.cpp
)

set(
	NORMSOURCES_EXECUTABLES
		${NORMSOURCES_EXECUTABLES_CATH_ASSIGN_DOMAINS}
		${NORMSOURCES_EXECUTABLES_CATH_CLUSTER}
		${NORMSOURCES_EXECUTABLES_CATH_EXTRACT_PDB}
		${NORMSOURCES_EXECUTABLES_CATH_MAP_CLUSTERS}
		${NORMSOURCES_EXECUTABLES_CATH_REFINE_ALIGN}
		${NORMSOURCES_EXECUTABLES_CATH_RESOLVE_HITS}
		${NORMSOURCES_EXECUTABLES_CATH_SCORE_ALIGN}
		${NORMSOURCES_EXECUTABLES_CATH_SSAP}
		${NORMSOURCES_EXECUTABLES_CATH_SUPERPOSE}
		${NORMSOURCES_EXECUTABLES_CHECK_PDB}
		${NORMSOURCES_EXECUTABLES_SNAP_JUDGEMENT}
)

set(
	NORMSOURCES
		${NORMSOURCES_CT_BIOCORE}
		${NORMSOURCES_CT_CATH_ASSIGN_DOMAINS}
		${NORMSOURCES_CT_CATH_CLUSTER}
		${NORMSOURCES_CT_CATH_REFINE_ALIGN}
		${NORMSOURCES_CT_CATH_SCORE_ALIGN}
		${NORMSOURCES_CT_CATH_SUPERPOSE}
		${NORMSOURCES_CT_CHOPPING}
		${NORMSOURCES_CT_CLUSTAGGLOM}
		${NORMSOURCES_CT_CLUSTER}
		${NORMSOURCES_CT_COMMON}
		${NORMSOURCES_CT_DISPLAY_COLOUR}
		${NORMSOURCES_CT_EXTERNAL_INFO}
		${NORMSOURCES_CT_OPTIONS}
		${NORMSOURCES_CT_RESOLVE_HITS}
		${NORMSOURCES_CT_SEQ}
		${NORMSOURCES_CT_TEST}
		${NORMSOURCES_CT_UNI}
		${NORMSOURCES_EXECUTABLES}
)

set(
	TESTSOURCES_CT_BIOCORE_CATH_BIOCORE
		ct_biocore/cath/biocore/chain_label_test.cpp
		ct_biocore/cath/biocore/residue_id_test.cpp
		ct_biocore/cath/biocore/residue_name_test.cpp
)

set(
	TESTSOURCES_CT_BIOCORE_CATH
		${TESTSOURCES_CT_BIOCORE_CATH_BIOCORE}
)

set(
	TESTSOURCES_CT_BIOCORE
		${TESTSOURCES_CT_BIOCORE_CATH}
)

set(
	TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK
		ct_cath_cluster/cath/cath_cluster/options/options_block/cath_cluster_input_options_block_test.cpp
)

set(
	TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC
		ct_cath_cluster/cath/cath_cluster/options/spec/cath_cluster_clustering_spec_test.cpp
)

set(
	TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS
		${TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK}
		${TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC}
)

set(
	TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER
		${TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS}
)

set(
	TESTSOURCES_CT_CATH_CLUSTER_CATH
		${TESTSOURCES_CT_CATH_CLUSTER_CATH_CATH_CLUSTER}
)

set(
	TESTSOURCES_CT_CATH_CLUSTER
		${TESTSOURCES_CT_CATH_CLUSTER_CATH}
)

set(
	TESTSOURCES_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN
		ct_cath_refine_align/cath/cath_refine_align/cath_align_refiner_test.cpp
)

set(
	TESTSOURCES_CT_CATH_REFINE_ALIGN_CATH
		${TESTSOURCES_CT_CATH_REFINE_ALIGN_CATH_CATH_REFINE_ALIGN}
)

set(
	TESTSOURCES_CT_CATH_REFINE_ALIGN
		${TESTSOURCES_CT_CATH_REFINE_ALIGN_CATH}
)

set(
	TESTSOURCES_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE
		ct_cath_superpose/cath/cath_superpose/cath_superposer_test.cpp
)

set(
	TESTSOURCES_CT_CATH_SUPERPOSE_CATH
		${TESTSOURCES_CT_CATH_SUPERPOSE_CATH_CATH_SUPERPOSE}
)

set(
	TESTSOURCES_CT_CATH_SUPERPOSE
		${TESTSOURCES_CT_CATH_SUPERPOSE_CATH}
)

set(
	TESTSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT
		ct_chopping/cath/chopping/chopping_format/sillitoe_chopping_format_test.cpp
		ct_chopping/cath/chopping/chopping_format/simple_chopping_format_test.cpp
)

set(
	TESTSOURCES_CT_CHOPPING_CATH_CHOPPING_DOMAIN
		ct_chopping/cath/chopping/domain/domain_test.cpp
)

set(
	TESTSOURCES_CT_CHOPPING_CATH_CHOPPING_REGION
		ct_chopping/cath/chopping/region/region_test.cpp
		ct_chopping/cath/chopping/region/regions_limiter_test.cpp
)

set(
	TESTSOURCES_CT_CHOPPING_CATH_CHOPPING
		${TESTSOURCES_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT}
		ct_chopping/cath/chopping/chopping_test.cpp
		${TESTSOURCES_CT_CHOPPING_CATH_CHOPPING_DOMAIN}
		${TESTSOURCES_CT_CHOPPING_CATH_CHOPPING_REGION}
)

set(
	TESTSOURCES_CT_CHOPPING_CATH
		${TESTSOURCES_CT_CHOPPING_CATH_CHOPPING}
)

set(
	TESTSOURCES_CT_CHOPPING
		${TESTSOURCES_CT_CHOPPING_CATH}
)

set(
	TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_DETAIL
		ct_clustagglom/cath/clustagglom/detail/clust_id_pot_test.cpp
)

set(
	TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE
		ct_clustagglom/cath/clustagglom/file/dissimilarities_file_test.cpp
		ct_clustagglom/cath/clustagglom/file/names_file_test.cpp
)

set(
	TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY
		ct_clustagglom/cath/clustagglom/hierarchy/hierarchy_group_test.cpp
		ct_clustagglom/cath/clustagglom/hierarchy/hierarchy_layer_test.cpp
		ct_clustagglom/cath/clustagglom/hierarchy/hierarchy_value_test.cpp
)

set(
	TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM
		ct_clustagglom/cath/clustagglom/calc_complete_linkage_merge_list_test.cpp
		ct_clustagglom/cath/clustagglom/clustagglom_fixture.cpp
		${TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_DETAIL}
		${TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE}
		${TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY}
		ct_clustagglom/cath/clustagglom/hierarchy_test.cpp
		ct_clustagglom/cath/clustagglom/link_list_test.cpp
		ct_clustagglom/cath/clustagglom/link_test.cpp
		ct_clustagglom/cath/clustagglom/links_test.cpp
		ct_clustagglom/cath/clustagglom/make_clusters_from_merges_test.cpp
		ct_clustagglom/cath/clustagglom/merge_test.cpp
)

set(
	TESTSOURCES_CT_CLUSTAGGLOM_CATH
		${TESTSOURCES_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM}
)

set(
	TESTSOURCES_CT_CLUSTAGGLOM
		${TESTSOURCES_CT_CLUSTAGGLOM_CATH}
)

set(
	TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_DETAIL
		ct_cluster/cath/cluster/detail/mapping_job_test.cpp
)

set(
	TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_FILE
		ct_cluster/cath/cluster/file/cluster_membership_file_test.cpp
)

set(
	TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_MAP
		ct_cluster/cath/cluster/map/aggregate_map_results_test.cpp
		ct_cluster/cath/cluster/map/map_results_test.cpp
		ct_cluster/cath/cluster/map/overlap_frac_distn_test.cpp
)

set(
	TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_TEST
		ct_cluster/cath/cluster/test/map_clusters_fixture.cpp
)

set(
	TESTSOURCES_CT_CLUSTER_CATH_CLUSTER
		ct_cluster/cath/cluster/cath_cluster_mapper_test.cpp
		ct_cluster/cath/cluster/cluster_entry_test.cpp
		ct_cluster/cath/cluster/cluster_info_test.cpp
		ct_cluster/cath/cluster/clusters_info_test.cpp
		${TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_DETAIL}
		ct_cluster/cath/cluster/domain_cluster_ids_by_seq_test.cpp
		ct_cluster/cath/cluster/domain_cluster_ids_test.cpp
		${TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_FILE}
		${TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_MAP}
		ct_cluster/cath/cluster/mapping_tool_test.cpp
		ct_cluster/cath/cluster/new_cluster_data_test.cpp
		ct_cluster/cath/cluster/old_cluster_data_test.cpp
		${TESTSOURCES_CT_CLUSTER_CATH_CLUSTER_TEST}
)

set(
	TESTSOURCES_CT_CLUSTER_CATH
		${TESTSOURCES_CT_CLUSTER_CATH_CLUSTER}
)

set(
	TESTSOURCES_CT_CLUSTER
		${TESTSOURCES_CT_CLUSTER_CATH}
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_ALGORITHM
		ct_common/cath/common/algorithm/are_same_test.cpp
		ct_common/cath/common/algorithm/constexpr_find_test.cpp
		ct_common/cath/common/algorithm/constexpr_floor_test.cpp
		ct_common/cath/common/algorithm/constexpr_integer_rounding_test.cpp
		ct_common/cath/common/algorithm/constexpr_is_uniq_test.cpp
		ct_common/cath/common/algorithm/constexpr_modulo_fns_test.cpp
		ct_common/cath/common/algorithm/for_n_test.cpp
		ct_common/cath/common/algorithm/transform_build_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_GRAPH
		ct_common/cath/common/boost_addenda/graph/spanning_tree_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR
		ct_common/cath/common/boost_addenda/range/adaptor/adaptor_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR
		ct_common/cath/common/boost_addenda/range/utility/iterator/cross_itr_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_UTILITY
		${TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR}
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE
		${TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_UTILITY}
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_TRIBOOL
		ct_common/cath/common/boost_addenda/tribool/tribool_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA
		${TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_GRAPH}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_TRIBOOL}
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_CLONE
		ct_common/cath/common/clone/clone_ptr_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_CONTAINER
		ct_common/cath/common/container/id_of_str_bidirnl_test.cpp
		ct_common/cath/common/container/id_of_string_test.cpp
		ct_common/cath/common/container/id_of_string_view_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_EXCEPTION
		ct_common/cath/common/exception/exception_is_equivalent_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_FILE
		ct_common/cath/common/file/ofstream_list_test.cpp
		ct_common/cath/common/file/open_fstream_test.cpp
		ct_common/cath/common/file/simple_file_read_write_test.cpp
		ct_common/cath/common/file/temp_file_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_GSL
		ct_common/cath/common/gsl/get_determinant_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_MATRIX
		ct_common/cath/common/matrix/matrix_index_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_METAPROGRAMMING
		ct_common/cath/common/metaprogramming/append_template_params_into_first_wrapper_test.cpp
		ct_common/cath/common/metaprogramming/change_template_subwrappers_test.cpp
		ct_common/cath/common/metaprogramming/change_template_wrapper_test.cpp
		ct_common/cath/common/metaprogramming/combine_params_lists_with_template_list_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_OPTIONAL
		ct_common/cath/common/optional/make_optional_if_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA
		ct_common/cath/common/rapidjson_addenda/rapidjson_writer_test.cpp
		ct_common/cath/common/rapidjson_addenda/string_of_rapidjson_write_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_STRING
		ct_common/cath/common/string/booled_to_string_test.cpp
		ct_common/cath/common/string/cath_to_string_test.cpp
		ct_common/cath/common/string/string_parse_tools_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_TUPLE
		ct_common/cath/common/tuple/make_tuple_with_skips_test.cpp
		ct_common/cath/common/tuple/mins_maxs_tuple_pair_mins_maxs_element_test.cpp
		ct_common/cath/common/tuple/tuple_increment_test.cpp
		ct_common/cath/common/tuple/tuple_lattice_index_test.cpp
		ct_common/cath/common/tuple/tuple_mins_maxs_element_test.cpp
		ct_common/cath/common/tuple/tuple_multiply_args_test.cpp
		ct_common/cath/common/tuple/tuple_subtract_test.cpp
		ct_common/cath/common/tuple/tuple_within_range_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON_TYPE_TRAITS
		ct_common/cath/common/type_traits/is_tuple_test.cpp
)

set(
	TESTSOURCES_CT_COMMON_CATH_COMMON
		${TESTSOURCES_CT_COMMON_CATH_COMMON_ALGORITHM}
		ct_common/cath/common/argc_argv_faker_test.cpp
		${TESTSOURCES_CT_COMMON_CATH_COMMON_BOOST_ADDENDA}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_CLONE}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_CONTAINER}
		ct_common/cath/common/difference_test.cpp
		${TESTSOURCES_CT_COMMON_CATH_COMMON_EXCEPTION}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_FILE}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_GSL}
		ct_common/cath/common/invert_permutation_test.cpp
		ct_common/cath/common/less_than_helper_test.cpp
		${TESTSOURCES_CT_COMMON_CATH_COMMON_MATRIX}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_METAPROGRAMMING}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_OPTIONAL}
		ct_common/cath/common/program_exception_wrapper_test.cpp
		${TESTSOURCES_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA}
		${TESTSOURCES_CT_COMMON_CATH_COMMON_STRING}
		ct_common/cath/common/temp_check_offset_1_test.cpp
		${TESTSOURCES_CT_COMMON_CATH_COMMON_TUPLE}
		ct_common/cath/common/type_to_string_test.cpp
		${TESTSOURCES_CT_COMMON_CATH_COMMON_TYPE_TRAITS}
)

set(
	TESTSOURCES_CT_COMMON_CATH
		${TESTSOURCES_CT_COMMON_CATH_COMMON}
)

set(
	TESTSOURCES_CT_COMMON
		${TESTSOURCES_CT_COMMON_CATH}
)

set(
	TESTSOURCES_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR
		ct_display_colour/cath/display_colour/display_colour_gradient_test.cpp
		ct_display_colour/cath/display_colour/display_colour_list_test.cpp
		ct_display_colour/cath/display_colour/display_colour_test.cpp
)

set(
	TESTSOURCES_CT_DISPLAY_COLOUR_CATH
		${TESTSOURCES_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR}
)

set(
	TESTSOURCES_CT_DISPLAY_COLOUR
		${TESTSOURCES_CT_DISPLAY_COLOUR_CATH}
)

set(
	TESTSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE
		ct_options/cath/options/executable/env_var_option_name_handler_test.cpp
)

set(
	TESTSOURCES_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK
		ct_options/cath/options/options_block/check_pdb_options_block_test.cpp
		ct_options/cath/options/options_block/detail_help_options_block_test.cpp
		ct_options/cath/options/options_block/extract_pdb_options_block_test.cpp
		ct_options/cath/options/options_block/misc_help_version_options_block_test.cpp
		ct_options/cath/options/options_block/options_block_test.cpp
)

set(
	TESTSOURCES_CT_OPTIONS_CATH_OPTIONS
		${TESTSOURCES_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE}
		${TESTSOURCES_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK}
)

set(
	TESTSOURCES_CT_OPTIONS_CATH
		${TESTSOURCES_CT_OPTIONS_CATH_OPTIONS}
)

set(
	TESTSOURCES_CT_OPTIONS
		${TESTSOURCES_CT_OPTIONS_CATH}
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO
		ct_resolve_hits/cath/resolve_hits/algo/masked_bests_cache_test.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE_DETAIL
		ct_resolve_hits/cath/resolve_hits/file/detail/hmmer_parser_test.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE
		ct_resolve_hits/cath/resolve_hits/file/cath_id_score_category_test.cpp
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE_DETAIL}
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT
		ct_resolve_hits/cath/resolve_hits/html_output/resolve_hits_html_outputter_test.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_filter_spec_test.cpp
		ct_resolve_hits/cath/resolve_hits/options/spec/crh_single_output_spec_test.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC}
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor_list_test.cpp
		ct_resolve_hits/cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor_test.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS_HITS_PROCESSOR}
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_RESOLVE
		ct_resolve_hits/cath/resolve_hits/resolve/hit_resolver_test.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TEST
		ct_resolve_hits/cath/resolve_hits/test/resolve_hits_fixture.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TRIM
		ct_resolve_hits/cath/resolve_hits/trim/resolve_boundary_test.cpp
		ct_resolve_hits/cath/resolve_hits/trim/seq_seg_boundary_fns_test.cpp
		ct_resolve_hits/cath/resolve_hits/trim/trim_spec_test.cpp
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO}
		ct_resolve_hits/cath/resolve_hits/calc_hit_list_test.cpp
		ct_resolve_hits/cath/resolve_hits/cath_hit_resolver_test.cpp
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE}
		ct_resolve_hits/cath/resolve_hits/first_hit_is_better_test.cpp
		ct_resolve_hits/cath/resolve_hits/full_hit_list_test.cpp
		ct_resolve_hits/cath/resolve_hits/full_hit_test.cpp
		ct_resolve_hits/cath/resolve_hits/hit_extras_test.cpp
		ct_resolve_hits/cath/resolve_hits/hit_test.cpp
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT}
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS}
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_READ_AND_PROCESS_HITS}
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_RESOLVE}
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TEST}
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_TRIM}
)

set(
	TESTSOURCES_CT_RESOLVE_HITS_CATH
		${TESTSOURCES_CT_RESOLVE_HITS_CATH_RESOLVE_HITS}
)

set(
	TESTSOURCES_CT_RESOLVE_HITS
		${TESTSOURCES_CT_RESOLVE_HITS_CATH}
)

set(
	TESTSOURCES_CT_SEQ_CATH_SEQ
		ct_seq/cath/seq/seq_seg_run_test.cpp
		ct_seq/cath/seq/seq_seg_test.cpp
)

set(
	TESTSOURCES_CT_SEQ_CATH
		${TESTSOURCES_CT_SEQ_CATH_SEQ}
)

set(
	TESTSOURCES_CT_SEQ
		${TESTSOURCES_CT_SEQ_CATH}
)

set(
	TESTSOURCES_CT_TEST_CATH_TEST_PREDICATE
		ct_test/cath/test/predicate/files_equal_test.cpp
		ct_test/cath/test/predicate/istreams_equal_test.cpp
)

set(
	TESTSOURCES_CT_TEST_CATH_TEST
		${TESTSOURCES_CT_TEST_CATH_TEST_PREDICATE}
)

set(
	TESTSOURCES_CT_TEST_CATH
		${TESTSOURCES_CT_TEST_CATH_TEST}
)

set(
	TESTSOURCES_CT_TEST
		${TESTSOURCES_CT_TEST_CATH}
)

set(
	TESTSOURCES_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER
		ct_uni/cath/acquirer/alignment_acquirer/align_refining_test.cpp
		ct_uni/cath/acquirer/alignment_acquirer/do_the_ssaps_alignment_acquirer_test.cpp
		ct_uni/cath/acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER
		ct_uni/cath/acquirer/pdbs_acquirer/pdbs_acquirer_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ACQUIRER
		${TESTSOURCES_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER}
		${TESTSOURCES_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER}
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY
		ct_uni/cath/alignment/common_residue_selection_policy/common_residue_selection_policy_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER
		ct_uni/cath/alignment/dyn_prog_align/detail/string_aligner/string_aligner_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL
		ct_uni/cath/alignment/dyn_prog_align/detail/return_path_matrix_test.cpp
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER}
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE
		ct_uni/cath/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_TEST
		ct_uni/cath/alignment/dyn_prog_align/test/dyn_prog_score_source_fixture.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DETAIL}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_SCORE_SOURCE}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_TEST}
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_GAP
		ct_uni/cath/alignment/gap/alignment_gap_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_IO_OUTPUTTER
		ct_uni/cath/alignment/io/outputter/horiz_align_outputter_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_IO
		ct_uni/cath/alignment/io/align_scaffold_test.cpp
		ct_uni/cath/alignment/io/alignment_io_test.cpp
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_IO_OUTPUTTER}
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER_DETAIL
		ct_uni/cath/alignment/refiner/detail/alignment_split_mapping_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER_DETAIL}
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL
		ct_uni/cath/alignment/residue_name_align/detail/residue_name_align_map_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN_DETAIL}
		ct_uni/cath/alignment/residue_name_align/residue_name_aligner_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_TEST
		ct_uni/cath/alignment/test/alignment_fixture.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT_TOOLS
		ct_uni/cath/alignment/tools/alignment_breaks_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_ALIGNMENT
		ct_uni/cath/alignment/alignment_action_test.cpp
		ct_uni/cath/alignment/alignment_coord_extractor_test.cpp
		ct_uni/cath/alignment/alignment_test.cpp
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_GAP}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_IO}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_REFINER}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_RESIDUE_NAME_ALIGN}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_TEST}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT_TOOLS}
)

set(
	TESTSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC
		ct_uni/cath/display/display_colour_spec/broad_display_colour_spec_test.cpp
		ct_uni/cath/display/display_colour_spec/display_colour_spec_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_DISPLAY_OPTIONS
		ct_uni/cath/display/options/display_spec_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL
		ct_uni/cath/display/viewer/pymol/pymol_tools_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_DISPLAY_VIEWER
		${TESTSOURCES_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL}
		ct_uni/cath/display/viewer/pymol_viewer_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_DISPLAY
		${TESTSOURCES_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC}
		${TESTSOURCES_CT_UNI_CATH_DISPLAY_OPTIONS}
		${TESTSOURCES_CT_UNI_CATH_DISPLAY_VIEWER}
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE_DSSP_WOLF
		ct_uni/cath/file/dssp_wolf/dssp_file_test.cpp
		ct_uni/cath/file/dssp_wolf/tally_residue_ids_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE_HMMER_SCORES_FILE
		ct_uni/cath/file/hmmer_scores_file/hmmer_scores_entry_test.cpp
		ct_uni/cath/file/hmmer_scores_file/hmmer_scores_file_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE_NAME_SET
		ct_uni/cath/file/name_set/name_set_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE_OPTIONS
		ct_uni/cath/file/options/data_dirs_options_block_test.cpp
		ct_uni/cath/file/options/data_dirs_spec_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE_PDB
		ct_uni/cath/file/pdb/coarse_element_type_test.cpp
		ct_uni/cath/file/pdb/element_type_string_test.cpp
		ct_uni/cath/file/pdb/pdb_atom_test.cpp
		ct_uni/cath/file/pdb/pdb_test.cpp
		ct_uni/cath/file/pdb/proximity_calculator_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE_PRC_SCORES_FILE
		ct_uni/cath/file/prc_scores_file/prc_scores_entry_test.cpp
		ct_uni/cath/file/prc_scores_file/prc_scores_file_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE_SSAP_SCORES_FILE
		ct_uni/cath/file/ssap_scores_file/ssap_scores_entry_test.cpp
		ct_uni/cath/file/ssap_scores_file/ssap_scores_file_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_FILE
		${TESTSOURCES_CT_UNI_CATH_FILE_DSSP_WOLF}
		${TESTSOURCES_CT_UNI_CATH_FILE_HMMER_SCORES_FILE}
		${TESTSOURCES_CT_UNI_CATH_FILE_NAME_SET}
		${TESTSOURCES_CT_UNI_CATH_FILE_OPTIONS}
		${TESTSOURCES_CT_UNI_CATH_FILE_PDB}
		${TESTSOURCES_CT_UNI_CATH_FILE_PRC_SCORES_FILE}
		${TESTSOURCES_CT_UNI_CATH_FILE_SSAP_SCORES_FILE}
)

set(
	TESTSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_TEST
		ct_uni/cath/outputter/alignment_outputter/test/alignment_outputter_fixture.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER
		ct_uni/cath/outputter/alignment_outputter/ssap_ostream_alignment_outputter_test.cpp
		${TESTSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_TEST}
)

set(
	TESTSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS
		ct_uni/cath/outputter/superposition_output_options/superposition_output_options_block_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUTTER
		ct_uni/cath/outputter/superposition_outputter/json_file_superposition_outputter_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_OUTPUTTER
		${TESTSOURCES_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER}
		${TESTSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUT_OPTIONS}
		${TESTSOURCES_CT_UNI_CATH_OUTPUTTER_SUPERPOSITION_OUTPUTTER}
)

set(
	TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY
		ct_uni/cath/scan/detail/check_scan/test_only/check_scan_on_final_alignment_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN
		${TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN_TEST_ONLY}
)

set(
	TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL_STRIDE
		ct_uni/cath/scan/detail/stride/co_stride_test.cpp
		ct_uni/cath/scan/detail/stride/rep_strider_test.cpp
		ct_uni/cath/scan/detail/stride/roled_scan_stride_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL
		${TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL_CHECK_SCAN}
		${TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL_STRIDE}
)

set(
	TESTSOURCES_CT_UNI_CATH_SCAN_SCAN_TOOLS
		ct_uni/cath/scan/scan_tools/all_vs_all_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCAN_SPATIAL_INDEX
		ct_uni/cath/scan/spatial_index/spatial_index_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCAN
		${TESTSOURCES_CT_UNI_CATH_SCAN_DETAIL}
		ct_uni/cath/scan/quad_criteria_test.cpp
		ct_uni/cath/scan/scan_index_test.cpp
		ct_uni/cath/scan/scan_stride_test.cpp
		${TESTSOURCES_CT_UNI_CATH_SCAN_SCAN_TOOLS}
		${TESTSOURCES_CT_UNI_CATH_SCAN_SPATIAL_INDEX}
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_DETAIL
		ct_uni/cath/score/aligned_pair_score/detail/score_common_coord_handler_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE
		ct_uni/cath/score/aligned_pair_score/ssap_score/ssap_score_accuracy_test.cpp
		ct_uni/cath/score/aligned_pair_score/ssap_score/ssap_score_post_processing_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX
		ct_uni/cath/score/aligned_pair_score/substitution_matrix/substitution_matrix_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE
		ct_uni/cath/score/aligned_pair_score/aligned_pair_score_test.cpp
		${TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_DETAIL}
		ct_uni/cath/score/aligned_pair_score/rmsd_score_test.cpp
		${TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE}
		${TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SUBSTITUTION_MATRIX}
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST
		ct_uni/cath/score/aligned_pair_score_list/aligned_pair_score_list_factory_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS
		ct_uni/cath/score/homcheck_tools/ssap_and_prc_test.cpp
		ct_uni/cath/score/homcheck_tools/ssaps_and_prcs_of_query_test.cpp
		ct_uni/cath/score/homcheck_tools/superfamily_of_domain_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_LENGTH_GETTER
		ct_uni/cath/score/length_getter/length_getter_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE
		ct_uni/cath/score/score_classification/label_pair_is_positive/label_pair_is_positive_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION
		${TESTSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION_LABEL_PAIR_IS_POSITIVE}
		ct_uni/cath/score/score_classification/rbf_model_test.cpp
		ct_uni/cath/score/score_classification/score_classn_value_better_value_test.cpp
		ct_uni/cath/score/score_classification/score_classn_value_results_set_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER
		ct_uni/cath/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG
		ct_uni/cath/score/true_pos_false_neg/classn_rate_stat_test.cpp
		${TESTSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_CLASSN_STAT_PLOTTER}
		ct_uni/cath/score/true_pos_false_neg/classn_stat_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SCORE
		${TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE}
		${TESTSOURCES_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST}
		${TESTSOURCES_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS}
		${TESTSOURCES_CT_UNI_CATH_SCORE_LENGTH_GETTER}
		${TESTSOURCES_CT_UNI_CATH_SCORE_SCORE_CLASSIFICATION}
		${TESTSOURCES_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG}
)

set(
	TESTSOURCES_CT_UNI_CATH_SSAP
		ct_uni/cath/ssap/distance_score_formula_test.cpp
		ct_uni/cath/ssap/selected_pair_test.cpp
		ct_uni/cath/ssap/ssap_test.cpp
		ct_uni/cath/ssap/windowed_matrix_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_ACCESSIBILITY_CALC
		ct_uni/cath/structure/accessibility_calc/dssp_accessibility_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_ENTRY_QUERIER
		ct_uni/cath/structure/entry_querier/entry_querier_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_GEOMETRY
		ct_uni/cath/structure/geometry/angle_test.cpp
		ct_uni/cath/structure/geometry/coord_list_test.cpp
		ct_uni/cath/structure/geometry/coord_test.cpp
		ct_uni/cath/structure/geometry/orient_test.cpp
		ct_uni/cath/structure/geometry/orientation_covering_test.cpp
		ct_uni/cath/structure/geometry/pca_test.cpp
		ct_uni/cath/structure/geometry/quat_rot_test.cpp
		ct_uni/cath/structure/geometry/restrict_to_single_linkage_extension_test.cpp
		ct_uni/cath/structure/geometry/rotation_test.cpp
		ct_uni/cath/structure/geometry/superpose_fit_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET
		ct_uni/cath/structure/protein/protein_source_file_set/protein_source_file_set_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN
		ct_uni/cath/structure/protein/amino_acid_test.cpp
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET}
		ct_uni/cath/structure/protein/residue_test.cpp
		ct_uni/cath/structure/protein/sec_struc_test.cpp
		ct_uni/cath/structure/protein/sec_struc_type_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST
		ct_uni/cath/structure/sec_struc_calc/dssp/test/dssp_dupl_fixture.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP
		ct_uni/cath/structure/sec_struc_calc/dssp/bifur_hbond_list_test.cpp
		ct_uni/cath/structure/sec_struc_calc/dssp/dssp_hbond_calc_test.cpp
		ct_uni/cath/structure/sec_struc_calc/dssp/dssp_ss_calc_test.cpp
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST}
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_SEC
		ct_uni/cath/structure/sec_struc_calc/sec/sec_calc_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_SEC}
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_DETAIL_PLATE
		ct_uni/cath/structure/view_cache/detail/plate/rod_cache_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_DETAIL
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_DETAIL_PLATE}
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER
		ct_uni/cath/structure/view_cache/filter/filter_vs_full_score_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL
		ct_uni/cath/structure/view_cache/index/detail/dims/detail/view_cache_index_dim_linear_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL}
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS}
		ct_uni/cath/structure/view_cache/index/detail/vcie_match_criteria_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL}
		ct_uni/cath/structure/view_cache/index/view_cache_index_entry_test.cpp
		ct_uni/cath/structure/view_cache/index/view_cache_index_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_DETAIL}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX}
)

set(
	TESTSOURCES_CT_UNI_CATH_STRUCTURE
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_ACCESSIBILITY_CALC}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_ENTRY_QUERIER}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_GEOMETRY}
		ct_uni/cath/structure/get_residue_names_test.cpp
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_PROTEIN}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE_VIEW_CACHE}
)

set(
	TESTSOURCES_CT_UNI_CATH_SUPERPOSITION_IO
		ct_uni/cath/superposition/io/superposition_io_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH_SUPERPOSITION
		${TESTSOURCES_CT_UNI_CATH_SUPERPOSITION_IO}
		ct_uni/cath/superposition/superposition_context_test.cpp
		ct_uni/cath/superposition/superposition_test.cpp
)

set(
	TESTSOURCES_CT_UNI_CATH
		${TESTSOURCES_CT_UNI_CATH_ACQUIRER}
		${TESTSOURCES_CT_UNI_CATH_ALIGNMENT}
		${TESTSOURCES_CT_UNI_CATH_DISPLAY}
		${TESTSOURCES_CT_UNI_CATH_FILE}
		${TESTSOURCES_CT_UNI_CATH_OUTPUTTER}
		${TESTSOURCES_CT_UNI_CATH_SCAN}
		${TESTSOURCES_CT_UNI_CATH_SCORE}
		${TESTSOURCES_CT_UNI_CATH_SSAP}
		${TESTSOURCES_CT_UNI_CATH_STRUCTURE}
		${TESTSOURCES_CT_UNI_CATH_SUPERPOSITION}
)

set(
	TESTSOURCES_CT_UNI
		${TESTSOURCES_CT_UNI_CATH}
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
	TESTSOURCES
		${TESTSOURCES_CT_BIOCORE}
		${TESTSOURCES_CT_CATH_CLUSTER}
		${TESTSOURCES_CT_CATH_REFINE_ALIGN}
		${TESTSOURCES_CT_CATH_SUPERPOSE}
		${TESTSOURCES_CT_CHOPPING}
		${TESTSOURCES_CT_CLUSTAGGLOM}
		${TESTSOURCES_CT_CLUSTER}
		${TESTSOURCES_CT_COMMON}
		${TESTSOURCES_CT_DISPLAY_COLOUR}
		${TESTSOURCES_CT_OPTIONS}
		${TESTSOURCES_CT_RESOLVE_HITS}
		${TESTSOURCES_CT_SEQ}
		${TESTSOURCES_CT_TEST}
		${TESTSOURCES_CT_UNI}
		${TESTSOURCES_EXECUTABLES}
)

##### DON'T EDIT THIS FILE - IT'S AUTO-GENERATED #####
