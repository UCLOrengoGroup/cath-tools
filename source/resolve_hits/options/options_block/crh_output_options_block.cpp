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

#include "crh_output_options_block.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>
#include <boost/range/join.hpp>

#include "common/algorithm/sort_uniq_build.hpp"
#include "common/clone/make_uptr_clone.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;

using boost::filesystem::path;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for an optional file to which the hits text should be output
const string crh_output_options_block::PO_HITS_TEXT_TO_FILE    { "hits-text-to-file"    };

/// \brief The option name for whether to suppress the default output of hits text to stdout
const string crh_output_options_block::PO_QUIET                { "quiet"                };

/// \brief The option name for whether to output the hits starts/stops *after* trimming
const string crh_output_options_block::PO_OUTPUT_TRIMMED_HITS  { "output-trimmed-hits"  };

/// \brief The option name for an optional file to which a summary of the input data should be output
const string crh_output_options_block::PO_SUMMARISE_TO_FILE    { "summarise-to-file"    };

/// \brief The option name for an optional file to which HTML should be output
const string crh_output_options_block::PO_HTML_OUTPUT_TO_FILE  { "html-output-to-file"  };

/// \brief The option name for an optional file to which JSON should be output
const string crh_output_options_block::PO_JSON_OUTPUT_TO_FILE  { "json-output-to-file"  };

/// \brief The option name for an optional file to which the CSS should be output
const string crh_output_options_block::PO_EXPORT_CSS_FILE      { "export-css-file"      };

/// \brief The option name for whether to output a summary of the HMMER alignment
const string crh_output_options_block::PO_OUTPUT_HMMER_ALN     { "output-hmmer-aln"     };

/// \brief A standard do_clone method
unique_ptr<options_block> crh_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string crh_output_options_block::do_get_block_name() const {
	return "Output ([...]-to-file options may be specified multiple times)";
}

/// \brief Add this block's options to the provided options_description
void crh_output_options_block::do_add_visible_options_to_description(options_description &arg_desc,           ///< The options_description to which the options are added
                                                                     const size_t        &/*arg_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                     ) {
	const string file_varname   { "<file>" };

	const auto hits_text_files_notifier     = [&] (const path_vec &x) { the_spec.set_hits_text_files     ( x           ); };
	const auto quiet_notifier               = [&] (const bool     &x) { the_spec.set_quiet               ( x           ); };
	const auto output_trimmed_hits_notifier = [&] (const bool     &x) {          set_output_trimmed_hits ( the_spec, x ); };
	const auto summarise_files_notifier     = [&] (const path_vec &x) { the_spec.set_summarise_files     ( x           ); };
	const auto html_output_files_notifier   = [&] (const path_vec &x) { the_spec.set_html_output_files   ( x           ); };
	const auto json_output_files_notifier   = [&] (const path_vec &x) { the_spec.set_json_output_files   ( x           ); };
	const auto export_css_file_notifier     = [&] (const path     &x) { the_spec.set_export_css_file     ( x           ); };

	arg_desc.add_options()
		(
			PO_HITS_TEXT_TO_FILE.c_str(),
			value<path_vec>()
				->value_name   ( file_varname                          )
				->notifier     ( hits_text_files_notifier              ),
			( "Write the resolved hits in plain text to file " + file_varname ).c_str()
		)
		(
			PO_QUIET.c_str(),
			bool_switch()
				->value_name   ( file_varname                          )
				->notifier     ( quiet_notifier                        )
				->default_value( crh_output_spec::DEFAULT_QUIET        ),
			"Suppress the default output of resolved hits in plain text to stdout"
		)
		(
			( PO_OUTPUT_TRIMMED_HITS ).c_str(),
			bool_switch()
				->notifier     ( output_trimmed_hits_notifier          )
				->default_value(
					means_output_trimmed_hits( crh_output_spec::DEFAULT_BOUNDARY_OUTPUT )
				),
			"When writing out the final hits, output the hits' starts/stop as they are *after trimming*"
		)
		(
			PO_SUMMARISE_TO_FILE.c_str(),
			value<path_vec>()
				->value_name   ( file_varname                          )
				->notifier     ( summarise_files_notifier              ),
			( "Write a brief text summary of the input data to file " + file_varname + " (or '-' for stdout)" ).c_str()
		)
		(
			PO_HTML_OUTPUT_TO_FILE.c_str(),
			value<path_vec>()
				->value_name   ( file_varname                          )
				->notifier     ( html_output_files_notifier            ),
			( "Write the results as HTML to file " + file_varname + " (or '-' for stdout)" ).c_str()
		)
		(
			PO_JSON_OUTPUT_TO_FILE.c_str(),
			value<path_vec>()
				->value_name   ( file_varname                          )
				->notifier     ( json_output_files_notifier            ),
			( "Write the results as JSON to file " + file_varname + " (or '-' for stdout)" ).c_str()
		)
		(
			PO_EXPORT_CSS_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                          )
				->notifier     ( export_css_file_notifier              ),
			( "Export the CSS used in the HTML output to " + file_varname + " (or '-' for stdout)").c_str()
		);

	static_assert( !                            crh_output_spec::DEFAULT_QUIET,
		"If crh_output_spec::DEFAULT_QUIET                isn't false, it might mess up the bool switch in here" );
	static_assert( ! means_output_trimmed_hits( crh_output_spec::DEFAULT_BOUNDARY_OUTPUT ),
		"If crh_segment_spec::DEFAULT_OUTPUT_TRIMMED_HITS isn't false, it might mess up the bool switch in here" );
}

/// \brief Add any hidden options to the provided options_description
void crh_output_options_block::do_add_hidden_options_to_description(options_description &arg_desc,       ///< The options_description to which the options are added
                                                                    const size_t        &arg_line_length ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                    ) {
	const auto output_hmmer_aln_notifier = [&] (const bool &x) { the_spec.set_output_hmmer_aln( x ); };

	deprecated_single_output_ob.add_hidden_options_to_description ( arg_desc, arg_line_length );
	deprecated_single_output_ob.add_visible_options_to_description( arg_desc, arg_line_length );

	arg_desc.add_options()
		(
			PO_OUTPUT_HMMER_ALN.c_str(),
			bool_switch()
				->notifier     ( output_hmmer_aln_notifier                 )
				->default_value( crh_output_spec::DEFAULT_OUTPUT_HMMER_ALN ),
			"[IN PUBLIC GENE3D COMMAND] Print a summary of the hmmer alignment in the output"
		);

	static_assert( ! crh_output_spec::DEFAULT_OUTPUT_HMMER_ALN,
		"If crh_output_spec::DEFAULT_OUTPUT_HMMER_ALN isn't false, it might mess up the bool switch in here" );
}

/// \brief Generate a description of any problem that makes the specified crh_output_options_block invalid
///        or none otherwise
str_opt crh_output_options_block::do_invalid_string(const variables_map &arg_variables_map ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                    ) const {
	const str_opt deprecated_invalid_str = deprecated_single_output_ob.invalid_string( arg_variables_map );
	if ( deprecated_invalid_str ) {
		return deprecated_invalid_str;
	}

	const str_vec specified_new  = specified_options           ( arg_variables_map, get_all_non_deprecated_option_names_that_clash_with_deprecated() );
	const str_vec specified_old  = specified_options_from_block( arg_variables_map, deprecated_single_output_ob                                      );

	if ( ! specified_old.empty() ) {
		if ( ! specified_new.empty() ) {
			return
				"Cannot mix old, deprecated options (--"
				+ boost::algorithm::join( specified_old, ", --" )
				+ ") with new, replacement options (--"
				+ boost::algorithm::join( specified_new, ", --" )
				+ "). Please use the new options only.";
		}

		BOOST_LOG_TRIVIAL( warning )
			<< "You're using deprecated output option(s) : --"
			<< boost::algorithm::join( specified_old, ", --" )
			<< ". Try replacing with : "
			<< get_deprecated_suggestion_str( deprecated_single_output_ob.get_crh_single_output_spec() );
	}

	return get_invalid_description( the_spec );
}

/// \brief Return all options names for this block
str_vec crh_output_options_block::do_get_all_options_names() const {
	// \TODO Remove the sort_uniq_build() (and replace with copy_build())
	// once the duplicated options have been removed from crh_single_output_options_block
	return sort_uniq_build<str_vec>(
		boost::range::join(
			get_all_non_deprecated_option_names(),
			deprecated_single_output_ob.get_all_options_names()
		)
	);
}

/// \brief Return all non-deprecated options names for this block that should clash with
///        deprecated options
///
/// \TODO Once the deprecated options have been removed, merge these functions' explicit options strings
///       into get_all_non_deprecated_option_names()
str_vec crh_output_options_block::get_all_non_deprecated_option_names_that_clash_with_deprecated() const {
	return {
		crh_output_options_block::PO_HITS_TEXT_TO_FILE,
		crh_output_options_block::PO_QUIET,
		crh_output_options_block::PO_SUMMARISE_TO_FILE,
		crh_output_options_block::PO_HTML_OUTPUT_TO_FILE,
		crh_output_options_block::PO_JSON_OUTPUT_TO_FILE,
	};
}
/// \brief Return all non-deprecated options names for this block that should clash with
///        deprecated options
///
/// \TODO Once the deprecated options have been removed, merge these functions' explicit options strings
///       into get_all_non_deprecated_option_names()
str_vec crh_output_options_block::get_all_non_deprecated_option_names_that_do_not_clash_with_deprecated() const {
	return {
		crh_output_options_block::PO_OUTPUT_TRIMMED_HITS,
		crh_output_options_block::PO_EXPORT_CSS_FILE,
		crh_output_options_block::PO_OUTPUT_HMMER_ALN,
	};
}

/// \brief Return all non-deprecated options names for this block
str_vec crh_output_options_block::get_all_non_deprecated_option_names() const {
	return sort_uniq_build<str_vec>(
		boost::range::join(
			get_all_non_deprecated_option_names_that_clash_with_deprecated(),
			get_all_non_deprecated_option_names_that_do_not_clash_with_deprecated()
		)
	);
}

/// \brief Getter for the crh_output_spec that the crh_output_options_block configures
const crh_output_spec & crh_output_options_block::get_crh_output_spec() const {
	return the_spec;
}

/// \brief Getter for the deprecated crh_single_output_options_block
const crh_single_output_options_block & crh_output_options_block::get_deprecated_single_output_options_block() const {
	return deprecated_single_output_ob;
}

/// \brief Get the deprecated crh_single_output_spec associated with the specified crh_output_options_block 
///
/// \relates crh_output_options_block
const crh_single_output_spec & cath::rslv::get_deprecated_single_output_spec(const crh_output_options_block &arg_crh_output_options_block ///< The crh_output_options_block to query
                                                                             ) {
	return arg_crh_output_options_block.get_deprecated_single_output_options_block().get_crh_single_output_spec();
}

/// \brief Return whether the specified crh_output_options_block implies any HTML output
///
/// \relates crh_output_options_block
bool cath::rslv::has_any_html_output(const crh_output_options_block &arg_crh_output_options_block ///< The crh_output_options_block to query
                                     ) {
	return (
		cath::rslv::has_html_output( arg_crh_output_options_block.get_crh_output_spec() )
		||
		get_out_format( arg_crh_output_options_block.get_deprecated_single_output_options_block() ) == crh_out_format::HTML
	);
}
