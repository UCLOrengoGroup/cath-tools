/// \file
/// \brief The clustmap_output_options_block class definitions

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

#include "clustmap_output_options_block.hpp"

#include "cath/common/clone/make_uptr_clone.hpp"

using namespace cath;
using namespace cath::clust;
using namespace cath::common;
using namespace cath::opts;

using boost::filesystem::path;
using boost::none;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for a batch ID to optionally append as a last column to the output
const string clustmap_output_options_block::PO_APPEND_BATCH_ID      { "append-batch-id"     };

/// \brief The option name for an optional file to which output should be redirected
const string clustmap_output_options_block::PO_OUTPUT_TO_FILE       { "output-to-file"      };

/// \brief The option name for an optional file to which a Markdown summary should be written
const string clustmap_output_options_block::PO_SUMMARISE_TO_FILE    { "summarise-to-file"   };

/// \brief The option name for whether to print out per-domain mapping info
const string clustmap_output_options_block::PO_PRINT_DOMAIN_MAPPING { "print-entry-results" };

/// \brief A standard do_clone method
unique_ptr<options_block> clustmap_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string clustmap_output_options_block::do_get_block_name() const {
	return "Output";
}

/// \brief Add this block's options to the provided options_description
void clustmap_output_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                          const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                          ) {
	const string id_varname     { "<id>"   };
	const string file_varname   { "<file>" };

	const auto append_batch_id_notifier      = [&] (const string &x) { the_spec.set_append_batch_id     ( x ); };
	const auto output_to_file_notifier       = [&] (const path   &x) { the_spec.set_output_to_file      ( x ); };
	const auto summarise_to_file_notifier    = [&] (const path   &x) { the_spec.set_summarise_to_file   ( x ); };
	const auto print_domain_mapping_notifier = [&] (const bool   &x) { the_spec.set_print_domain_mapping( x ); };

	prm_desc.add_options()
		(
			PO_APPEND_BATCH_ID.c_str(),
			value<string>()
				->value_name   ( id_varname                                         )
				->notifier     ( append_batch_id_notifier                           ),
			( "Append batch ID " + id_varname + " as an extra column in the results output "
				"(equivalent to the first column in a --multi-batch-file input file)" ).c_str()
		)
		(
			PO_OUTPUT_TO_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                                       )
				->notifier     ( output_to_file_notifier                            ),
			( "Write output to file " + file_varname + " (or, if unspecified, to stdout)" ).c_str()
		)
		(
			PO_SUMMARISE_TO_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                                       )
				->notifier     ( summarise_to_file_notifier                         ),
			( "Print a summary of the renumbering to file " + file_varname ).c_str()
		)
		(
			( PO_PRINT_DOMAIN_MAPPING ).c_str(),
			bool_switch()
				->notifier     ( print_domain_mapping_notifier                      )
				->default_value( clustmap_output_spec::DEFAULT_PRINT_DOMAIN_MAPPING ),
			"Output the entry (domain)-level mapping results"
		);

		static_assert( ! clustmap_output_spec::DEFAULT_PRINT_DOMAIN_MAPPING,
			"If clustmap_output_spec::DEFAULT_PRINT_DOMAIN_MAPPING isn't false, it might mess up the bool switch in here" );
}

/// \brief Generate a description of any problem that makes the specified clustmap_output_options_block invalid
///        or none otherwise
str_opt clustmap_output_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                         ) const {
	return none;
}

/// \brief Return all options names for this block
str_vec clustmap_output_options_block::do_get_all_options_names() const {
	return {
		clustmap_output_options_block::PO_APPEND_BATCH_ID,
		clustmap_output_options_block::PO_OUTPUT_TO_FILE,
		clustmap_output_options_block::PO_SUMMARISE_TO_FILE,
		clustmap_output_options_block::PO_PRINT_DOMAIN_MAPPING,
	};
}

/// \brief Getter for the clustmap_output_spec that the clustmap_output_options_block configures
const clustmap_output_spec & clustmap_output_options_block::get_clustmap_output_spec() const {
	return the_spec;
}
