/// \file
/// \brief The clustmap_options class definitions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "clustmap_options.hpp"

#include <boost/algorithm/string/join.hpp>

using namespace cath;
using namespace cath::clust;
using namespace cath::opts;
using namespace std::literals::string_literals;

using boost::algorithm::join;
using boost::none;
using boost::program_options::positional_options_description;
using boost::program_options::variables_map;
using std::string;

/// The name of the program that uses this executable_options
const string clustmap_options::PROGRAM_NAME("cath-map-clusters");

/// \brief Get the options for the "Detailed Help" block
str_str_str_pair_map clustmap_options::detail_help_spec() {
	return {
		{ "sorting-help", { "Show the criteria for sorting unmapped clusters", get_cmc_sorting_criteria_help_string() } },
	};
}

/// \brief Get the name of the program that uses this executable_options
string clustmap_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief Get the positional options, which in this case is the input block's PO_INPUT_FILE_OR_STDIN option
positional_options_description clustmap_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add( clustmap_input_options_block::PO_WORKING_CLUSTMEMB_FILE.c_str(), 1 );
	return positionals;
}

/// \brief Review all specified options and return a string containing any errors or a help string
///        (possibly using a description of all visible options)
///
/// This is a concrete definition of a virtual method that's pure in executable_options
///
/// This should only be called by executable_options, as the last step of the parse_options()
/// method, after all real parsing has completed.
///
/// \pre The options should have been parsed
///
/// \returns Any error/help string arising from the newly specified options
///          or an empty string if there aren't any
str_opt clustmap_options::do_get_error_or_help_string() const {
	// If detailed help was requested, then provide it
	if ( the_detail_help_ob.has_help_string() ) {
		return the_detail_help_ob.help_string();
	}

	const variables_map &local_vm = get_variables_map();

	if ( ! specifies_option( local_vm, clustmap_input_options_block::PO_WORKING_CLUSTMEMB_FILE ) ) {
		return "Must specify an input file"s;
	}

	const bool     read_batches_from_input = get_clustmap_input_spec().get_read_batches_from_input();
	const path_opt map_from_clustmemb_file = get_clustmap_input_spec().get_map_from_clustmemb_file();

	if ( get_clustmap_output_spec().get_append_batch_id() && read_batches_from_input ) {
		return "Cannot specify a batch ID for appending (--"
			+ clustmap_output_options_block::PO_APPEND_BATCH_ID
			+ ") when reading batches from input (--"
			+ clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT
			+ ")";
	}

	if ( specified_clust_thresh_options( local_vm ) && ! map_from_clustmemb_file && ! read_batches_from_input ) {
		return "Cannot specify mapping threshold options (--"
			+ join( clust_thresh_option_names(), ", --" )
			+ ") if not specifying map-from file (--"
			+ clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE
			+ ") and not using batch mode (--"
			+ clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT
			+ ")";
	}

	// If no error or help string, then return none
	return none;
}

/// \brief Get a string to prepend to the standard help
string clustmap_options::do_get_help_prefix_string() const {
	return "Usage: " + PROGRAM_NAME + R"( [options] <input_file>

)" + get_overview_string() + R"(

When <input_file> is -, the input is read from standard input.)";
}

/// \brief Get a string to append to the standard help
string clustmap_options::do_get_help_suffix_string() const {
	return R"(
The input cluster-membership data should contain lines like:

cluster_name domain_id

...where domain_id is a sequence/protein name,)"
	" optionally suffixed with the domain's segments in notation like '/100-199,350-399'."
	" The suffixes should be present for all of the domain IDs or for none of them. "
	R"(Domains shouldn't overlap with others in the same cluster-membership data.

Input data doesn't have to be grouped by cluster.

NOTE: When mapping (ie using --)" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE
	+ "), the mapping algorithm treats the two cluster membership files differently - see --"
	+ clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL
	+ " description.\n";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string clustmap_options::do_get_overview_string() const {
	return R"(Map names from previous clusters to new clusters based on (the overlaps between)
their members (which may be specified as regions within a parent sequence).
Renumber any clusters with no equivalents.)";
}

/// \brief Ctor, which initialises the detail_help_ob and adds the options_blocks to the parent executable_options
clustmap_options::clustmap_options() : the_detail_help_ob( detail_help_spec() ) {
	super::add_options_block( the_input_ob       );
	super::add_options_block( the_mapping_ob     );
	super::add_options_block( the_output_ob      );
	super::add_options_block( the_detail_help_ob );
}

/// \brief Getter for the spec of the cath-map-clusters input options_block
const clustmap_input_spec & clustmap_options::get_clustmap_input_spec() const {
	return the_input_ob.get_clustmap_input_spec();
}

/// \brief Getter for the spec of the cath-map-clusters mapping options_block
const clust_mapping_spec & clustmap_options::get_clust_mapping_spec() const {
	return the_mapping_ob.get_clust_mapping_spec();
}

/// \brief Getter for the spec of the cath-map-clusters output options_block
const clustmap_output_spec & clustmap_options::get_clustmap_output_spec() const {
	return the_output_ob.get_clustmap_output_spec();
}

/// \brief Return the string containing help on CMC sorting criteria
string cath::clust::get_cmc_sorting_criteria_help_string() {
	return R"(The sorting criteria for new, unmapped clusters are:

 * Descending on sum over domains of sqrt(total_dom_length) (ie clusters with more/longer sequences come earlier, with more emphasis on having more sequences)
 * Descending on number of sequences (ie clusters with more sequences come first)
 * Ascending on average mid-point index (ie clusters with domains earlier in their sequences come first)
 * Ascending on first domain ID
)";
}
