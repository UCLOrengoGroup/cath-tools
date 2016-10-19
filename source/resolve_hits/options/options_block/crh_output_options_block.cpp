/// \file
/// \brief The crh_output_options_block class definitions

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

#include "crh_output_options_block.h"

#include "common/clone/make_uptr_clone.h"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;

using boost::filesystem::path;
using boost::none;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for the output file to which data should be written (or unspecified for stdout)
const string crh_output_options_block::PO_OUTPUT_FILE          { "output-file"         };

/// \brief The option name for whether to output the hits starts/stops *after* trimming
const string crh_output_options_block::PO_OUTPUT_TRIMMED_HITS  { "output-trimmed-hits" };

/// \brief The option name for whether to output HTML describing the hits and the results
const string crh_output_options_block::PO_GENERATE_HTML_OUTPUT { "html-output"         };

/// \brief A standard do_clone method
unique_ptr<options_block> crh_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string crh_output_options_block::do_get_block_name() const {
	return "Output";
}

/// \brief Add this block's options to the provided options_description
void crh_output_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                     ) {
	const string file_varname   { "<file>"   };

	const auto output_file_notifier          = [&] (const path &x) { the_spec.set_output_file         (           x ); };
	const auto output_trimmed_hits_notifier  = [&] (const bool &x) {          set_output_trimmed_hits ( the_spec, x ); };
	const auto generate_html_output_notifier = [&] (const bool &x) { the_spec.set_generate_html_output(           x ); };

	arg_desc.add_options()
		(
			PO_OUTPUT_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                                  )
				->notifier     ( output_file_notifier                          ),
			( "Write output to file " + file_varname + " (or, if unspecified, to stdout)" ).c_str()
		)
		(
			( PO_OUTPUT_TRIMMED_HITS ).c_str(),
			bool_switch()
				->notifier     ( output_trimmed_hits_notifier                              )
				->default_value(
					means_output_trimmed_hits( crh_output_spec::DEFAULT_BOUNDARY_OUTPUT )
				),
			"When writing out the final hits, output the hits' starts/stop as they are *after trimming*"
		)
		(
			( PO_GENERATE_HTML_OUTPUT ).c_str(),
			bool_switch()
				->notifier     ( generate_html_output_notifier                 )
				->default_value( crh_output_spec::DEFAULT_GENERATE_HTML_OUTPUT ),
			"Output the results as HTML"
		);

	static_assert( ! means_output_trimmed_hits( crh_output_spec::DEFAULT_BOUNDARY_OUTPUT ),
		"If crh_segment_spec::DEFAULT_OUTPUT_TRIMMED_HITS isn't false, it might mess up the bool switch in here" );
	static_assert( !                            crh_output_spec::DEFAULT_GENERATE_HTML_OUTPUT,
		"If crh_output_spec::DEFAULT_GENERATE_HTML_OUTPUT isn't false, it might mess up the bool switch in here" );
}

/// \brief Generate a description of any problem that makes the specified crh_output_options_block invalid
///        or none otherwise
str_opt crh_output_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                    ) const {
	return get_invalid_description( the_spec );
}

/// \brief Getter for the crh_output_spec that the crh_output_options_block configures
const crh_output_spec & crh_output_options_block::get_crh_output_spec() const {
	return the_spec;
}
