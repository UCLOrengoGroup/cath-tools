/// \file
/// \brief The cath_cluster_options class definitions

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

#include "cath_cluster_options.hpp"

using namespace ::cath;
using namespace ::cath::clust;

using ::boost::none;
using ::boost::program_options::positional_options_description;
using ::std::string;

/// The name of the program that uses this executable_options
const string cath_cluster_options::PROGRAM_NAME("cath-cluster");

/// \brief Get the name of the program that uses this executable_options
string cath_cluster_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief Get the positional options, which in this case is the input block's PO_LINKS_INFILE option
positional_options_description cath_cluster_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add( cath_cluster_input_options_block::PO_LINKS_INFILE.c_str(), 1 );
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
str_opt cath_cluster_options::do_get_error_or_help_string() const {
	const auto orslc = has_output_requiring_single_level_clustering( get_cath_cluster_output_spec() );
	if ( orslc && get_cath_cluster_clustering_spec().get_levels().size() > 1 ) {
		return "Cannot use --"
			+ *orslc
			+ " output when generating more than one level of clustering";
	}

	return none;
}

/// \brief Get a string to prepend to the standard help
string cath_cluster_options::do_get_help_prefix_string() const {
	return "Usage: "
		+ PROGRAM_NAME
		+ " --"
		+ cath_cluster_input_options_block::PO_LINK_DIRN
		+ R"( <dirn> --)"
		+ cath_cluster_clustering_options_block::PO_LEVELS
		+ R"( <levels> [options] <input_file>

)" + get_overview_string() + R"(

When <input_file> is -, the links are read from standard input.

The clustering is complete-linkage.)";
}

/// \brief Get a string to append to the standard help
string cath_cluster_options::do_get_help_suffix_string() const {
	return R"(
Links input format: `id1 id2 other columns afterwards`
...where --)" + cath_cluster_input_options_block::PO_COLUMN_IDX + R"( can be used to specify the column that contains the values

Names input format: `id score`
...where score is used to sort such that lower-scored entries appear earlier
)";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string cath_cluster_options::do_get_overview_string() const {
	return R"(Cluster items based on the links between them.)";
}

/// \brief Ctor, which initialises the detail_help_ob and adds the options_blocks to the parent executable_options
cath_cluster_options::cath_cluster_options() {
	super::add_options_block( the_input_ob  );
	super::add_options_block( clustering_ob  );
	super::add_options_block( the_output_ob );
}

/// \brief Getter for the cath-resolve-hits input options_block
const cath_cluster_input_spec & cath_cluster_options::get_cath_cluster_input_spec() const {
	return the_input_ob.get_cath_cluster_input_spec();
}

/// \brief Getter for the cath-resolve-hits clustering options_block
const cath_cluster_clustering_spec & cath_cluster_options::get_cath_cluster_clustering_spec() const {
	return clustering_ob.get_cath_cluster_clustering_spec();
}

/// \brief Getter for the cath-resolve-hits output options_block
const cath_cluster_output_spec & cath_cluster_options::get_cath_cluster_output_spec() const {
	return the_output_ob.get_cath_cluster_output_spec();
}
