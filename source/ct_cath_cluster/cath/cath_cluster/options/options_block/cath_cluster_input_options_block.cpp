/// \file
/// \brief The cath_cluster_input_options_block class definitions

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

#include "cath_cluster_input_options_block.hpp"

#include <filesystem>

#include <boost/algorithm/string/join.hpp>

#include "cath/common/boost_addenda/program_options/layout_values_with_descs.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/program_options/prog_opt_num_range.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::algorithm::join;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::filesystem::path;
using ::std::numeric_limits;
using ::std::string;
using ::std::unique_ptr;

/// \brief The option name for an optional file from which links should be read
const string cath_cluster_input_options_block::PO_LINKS_INFILE { "links-infile" };

/// \brief The option name for the direction of links in the input
const string cath_cluster_input_options_block::PO_LINK_DIRN    { "link-dirn"    };

/// \brief The option name for the index of the column from which the link values are to be read
const string cath_cluster_input_options_block::PO_COLUMN_IDX   { "column-idx"   };

/// \brief The option name for an optional file from which names should be read
const string cath_cluster_input_options_block::PO_NAMES_INFILE { "names-infile" };

/// \brief A standard do_clone method
unique_ptr<options_block> cath_cluster_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string cath_cluster_input_options_block::do_get_block_name() const {
	return "Input";
}

/// \brief Add this block's options to the provided options_description
void cath_cluster_input_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                             const size_t        &/*prm_line_length*/ ///< The line length to be used when inputting the description (not very clearly documented in Boost)
                                                                             ) {
	const auto &sep     = SUB_DESC_SEPARATOR;
	const auto &sub_sep = SUB_DESC_PAIR_SEPARATOR;

	const string link_dirn_varname  { "<dirn>"   };
	const string column_idx_varname { "<colnum>" };
	const string file_varname       { "<file>"   };

	const auto link_dirn_notifier    = [&] (const link_dirn &x) { the_spec.set_link_dirn   ( x     ); };
	const auto column_idx_notifier   = [&] (const size_t    &x) { the_spec.set_column_idx  ( x - 1 ); }; // Subtract 1 to move from offset-1 to offset-0
	const auto names_infile_notifier = [&] (const path      &x) { the_spec.set_names_infile( x     ); };

	const str_vec link_dirn_descs = layout_values_with_descs(
		all_link_dirns,
		[] (const link_dirn &x) { return to_string( x ); },
		&description_of_link_dirn,
		sub_sep
	);

	prm_desc.add_options()
		(
			PO_LINK_DIRN.c_str(),
			value<link_dirn>()
				->value_name   ( link_dirn_varname     )
				->notifier     ( link_dirn_notifier    )
				->required     (                       ),
			(   "Interpret each link value as "
			  + link_dirn_varname
			  + ", one of:"
			  + sep
			  + join( link_dirn_descs, sep ) ).c_str()
		)
		(
			PO_COLUMN_IDX.c_str(),
			value< prog_opt_num_range<size_t, 3, numeric_limits<uint32_t>::max(), size_t> >()
				->value_name   ( column_idx_varname    )
				->notifier     ( column_idx_notifier   )
				->default_value( 3                     ),
			(   "Parse the link values (distances/strengths) from column number "
			  + column_idx_varname
			  + "\nMust be ≥ 3 because columns 1 and 2 must contain the IDs" ).c_str()
		)
		(
			PO_NAMES_INFILE.c_str(),
			value<path>()
				->value_name   ( file_varname          )
				->notifier     ( names_infile_notifier ),
			(   "[RECOMMENDED] Read names and sorting scores from file "
			  + file_varname
			  + " (or '-' for stdin)" ).c_str()
		);
}

/// \brief Add any hidden options to the provided options_description
void cath_cluster_input_options_block::do_add_hidden_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                            const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                            ) {
	const string file_varname   { "<input_file>" };
	const auto links_infile_notifier = [&] (const path &x) { the_spec.set_links_infile( x ); };
	prm_desc.add_options()
		(
			PO_LINKS_INFILE.c_str(),
			value<path>()
				->value_name   ( file_varname          )
				->notifier     ( links_infile_notifier ),
			(   "Read the input links from file "
			  + file_varname
			  + "\nThis is the single positional option"
			  + "\nWhen "
			  + file_varname
			  + " is -, read from standard input" ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified cath_cluster_input_options_block invalid
///        or none otherwise
str_opt cath_cluster_input_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                             ) const {
	return get_invalid_description( the_spec );
}

/// \brief Return all options names for this block
str_vec cath_cluster_input_options_block::do_get_all_options_names() const {
	return {
		cath_cluster_input_options_block::PO_LINKS_INFILE,
		cath_cluster_input_options_block::PO_NAMES_INFILE,
	};
}

/// \brief Getter for the cath_cluster_input_spec that the cath_cluster_input_options_block configures
const cath_cluster_input_spec & cath_cluster_input_options_block::get_cath_cluster_input_spec() const {
	return the_spec;
}
