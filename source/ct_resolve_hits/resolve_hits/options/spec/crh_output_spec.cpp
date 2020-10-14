/// \file
/// \brief The crh_output_spec class definitions

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

#include "crh_output_spec.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/count.hpp>

#include "common/algorithm/append.hpp"
#include "common/algorithm/contains.hpp"
#include "common/algorithm/is_uniq.hpp"
#include "common/algorithm/sort_uniq_copy.hpp"
#include "resolve_hits/options/options_block/crh_output_options_block.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath;
using namespace std::literals::string_literals;

using boost::filesystem::path;
using boost::none;
using boost::numeric_cast;
using boost::range::count;
using std::string;

constexpr bool                crh_output_spec::DEFAULT_QUIET;
constexpr hit_boundary_output crh_output_spec::DEFAULT_BOUNDARY_OUTPUT;
constexpr bool                crh_output_spec::DEFAULT_OUTPUT_HMMER_ALN;

/// \brief Getter for any files to which hits text should be output
const path_vec & crh_output_spec::get_hits_text_files() const {
	return hits_text_files;
}

/// \brief Getter for whether to suppress the default output of hits text to stdout
const bool & crh_output_spec::get_quiet() const {
	return quiet;
}

/// \brief Getter for whether to output the hits starts/stops *after* trimming
const hit_boundary_output & crh_output_spec::get_boundary_output() const {
	return boundary_output;
}

/// \brief Getter for any files to which a summary of the input should be output
const path_vec & crh_output_spec::get_summarise_files() const {
	return summarise_files;
}

/// \brief Getter for any files to which HTML should be output
const path_vec & crh_output_spec::get_html_output_files() const {
	return html_output_files;
}

/// \brief Getter for any files to which JSON should be output
const path_vec & crh_output_spec::get_json_output_files() const {
	return json_output_files;
}

/// \brief Getter for any files to which the HTML's CSS should be output
const path_opt & crh_output_spec::get_export_css_file() const {
	return export_css_file;
}

/// \brief Getter for whether to output a summary of the HMMER alignment
const bool & crh_output_spec::get_output_hmmer_aln() const {
	return output_hmmer_aln;
}

/// \brief Setter for any files to which hits text should be output
crh_output_spec & crh_output_spec::set_hits_text_files(const path_vec &prm_hits_text_files ///< Any files to which hits text should be output
                                                       ) {
	hits_text_files = prm_hits_text_files;
	return *this;
}

/// \brief Setter for whether to suppress the default output of hits text to stdout
crh_output_spec & crh_output_spec::set_quiet(const bool &prm_quiet ///< Whether to suppress the default output of hits text to stdout
                                             ) {
	quiet = prm_quiet;
	return *this;
}

/// \brief Setter for whether to output the hits starts/stops *after* trimming
crh_output_spec & crh_output_spec::set_boundary_output(const hit_boundary_output &prm_boundary_output ///< Whether to output the hits starts/stops *after* trimming
                                                       ) {
	boundary_output = prm_boundary_output;
	return *this;
}

/// \brief Setter for any files to which a summary of the input should be output
crh_output_spec & crh_output_spec::set_summarise_files(const path_vec &prm_summarise_files ///< Any files to which a summary of the input should be output
                                                       ) {
	summarise_files = prm_summarise_files;
	return *this;
}

/// \brief Setter for any files to which HTML should be output
crh_output_spec & crh_output_spec::set_html_output_files(const path_vec &prm_html_output_files ///< Any files to which HTML should be output
                                                         ) {
	html_output_files = prm_html_output_files;
	return *this;
}

/// \brief Setter for any files to which JSON should be output
crh_output_spec & crh_output_spec::set_json_output_files(const path_vec &prm_json_output_files ///< Any files to which JSON should be output
                                                         ) {
	json_output_files = prm_json_output_files;
	return *this;
}

/// \brief Setter for any files to which the HTML's CSS should be output
crh_output_spec & crh_output_spec::set_export_css_file(const path_opt &prm_export_css_file ///< Any files to which the HTML's CSS should be output
                                                       ) {
	export_css_file = prm_export_css_file;
	return *this;
}

/// \brief Setter for whether to output a summary of the HMMER alignment
crh_output_spec & crh_output_spec::set_output_hmmer_aln(const bool &prm_output_hmmer_aln ///< Whether to output a summary of the HMMER alignment
                                                        ) {
	output_hmmer_aln = prm_output_hmmer_aln;
	return *this;
}

/// \brief Return whether the specified crh_output_spec implies any HTML output
///
/// \relates crh_output_spec
bool cath::rslv::has_html_output(const crh_output_spec &prm_output_spec ///< The crh_output_spec to query
                                 ) {
	return ! prm_output_spec.get_html_output_files().empty();
}

/// \brief Return whether the specified crh_output_spec implies any hits-text output
///
/// \relates crh_output_spec
bool cath::rslv::has_hits_text_output(const crh_output_spec &prm_output_spec ///< The crh_output_spec to query
                                      ) {
	return ! prm_output_spec.get_hits_text_files().empty();
}

/// \brief Return whether the specified crh_output_spec has any output files that match the specified file
///
/// \relates crh_output_spec
bool cath::rslv::has_any_out_files_matching(const crh_output_spec &prm_output_spec, ///< The crh_output_spec to query
                                            const path            &prm_query_path   ///< The file being searched for
                                            ) {
	return (
		contains( prm_output_spec.get_hits_text_files(),   prm_query_path )
		||
		contains( prm_output_spec.get_summarise_files(),   prm_query_path )
		||
		contains( prm_output_spec.get_html_output_files(), prm_query_path )
		||
		contains( prm_output_spec.get_json_output_files(), prm_query_path )
		||
		( prm_output_spec.get_export_css_file() == prm_query_path )
	);
}

/// \brief Get all output paths implied by the specified crh_output_spec
///
/// \relates crh_output_spec
path_vec cath::rslv::get_all_output_paths(const crh_output_spec &prm_output_spec ///< The crh_output_spec to query
                                          ) {
	path_vec the_paths;
	append( the_paths, prm_output_spec.get_hits_text_files()   );
	append( the_paths, prm_output_spec.get_summarise_files()   );
	append( the_paths, prm_output_spec.get_html_output_files() );
	append( the_paths, prm_output_spec.get_json_output_files() );
	if ( prm_output_spec.get_export_css_file() ) {
		the_paths.push_back( *prm_output_spec.get_export_css_file() );
	}
	return the_paths;
}

/// \brief Generate a description of any problem that makes the specified crh_output_spec invalid
///        or none otherwise
///
/// \relates crh_output_spec
str_opt cath::rslv::get_invalid_description(const crh_output_spec &prm_output_spec ///< The crh_output_spec to query
                                            ) {
	const auto all_sorted_paths = sort_copy( get_all_output_paths( prm_output_spec ) );

	const size_t num_stdouts = numeric_cast<size_t>( count( all_sorted_paths, path{ "-" } ) );
	if ( num_stdouts > 1 ) {
		return "Cannot send more than one type of output to stdout (which is specified as file \"-\")"s;
	}

	if ( num_stdouts > 0 && prm_output_spec.get_quiet() ) {
		return "Cannot send output to stdout when requesting --" + crh_output_options_block::PO_QUIET;
	}

	if ( ! is_uniq( all_sorted_paths ) ) {
		return "Cannot send more than one type of output to the same output file"s;
	}

	const bool hits_text_to_stdout = num_stdouts == 0 && ! prm_output_spec.get_quiet();
	if ( means_output_trimmed_hits( prm_output_spec.get_boundary_output() ) && ! ( hits_text_to_stdout || has_hits_text_output( prm_output_spec ) ) ) {
		return
			"Cannot specify trimmed boundaries if not outputting any hits text (with --"
			+ crh_output_options_block::PO_HITS_TEXT_TO_FILE
			+ ") but your output format may already contain the information you require";
	}

	return none;
}

/// \brief Convenience setter for whether to output the hits starts/stops *after* trimming
///
/// \relates crh_output_spec
crh_output_spec & cath::rslv::set_output_trimmed_hits(crh_output_spec &prm_output_spec,        ///< The crh_output_spec to modify
                                                      const bool      &prm_output_trimmed_hits ///< Whether to output the hits starts/stops *after* trimming
                                                      ) {
	return prm_output_spec.set_boundary_output( hit_boundary_output_of_output_trimmed_hits( prm_output_trimmed_hits ) );
}
