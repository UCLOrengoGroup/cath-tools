/// \file
/// \brief The crh_single_output_spec class definitions

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

#include "crh_single_output_spec.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/count_if.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::filesystem::path;
using boost::make_optional;
using boost::none;
using boost::range::count_if;
using std::string;
using std::vector;

constexpr hit_boundary_output crh_single_output_spec::DEFAULT_BOUNDARY_OUTPUT;
constexpr bool                crh_single_output_spec::DEFAULT_GENERATE_HTML_OUTPUT;
constexpr bool                crh_single_output_spec::DEFAULT_JSON_OUTPUT;
constexpr bool                crh_single_output_spec::DEFAULT_SUMMARISE;
constexpr bool                crh_single_output_spec::DEFAULT_OUTPUT_HMMSEARCH_ALN;

/// \brief Generate a string describing the specified crh_out_format
///
/// \relates crh_out_format
string cath::rslv::to_string(const crh_out_format &arg_out_format ///< The crh_out_format to describe
                             ) {
	switch ( arg_out_format ) {
		case ( crh_out_format::STANDARD ) : { return "standard"  ; }
		case ( crh_out_format::SUMMARY  ) : { return "a summary" ; }
		case ( crh_out_format::HTML     ) : { return "HTML"      ; }
		case ( crh_out_format::JSON     ) : { return "JSON"      ; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of crh_out_format not recognised whilst converting to_string()"));
}

/// \brief Getter for the output file to which data should be written
const path_opt & crh_single_output_spec::get_output_file() const {
	return output_file;
}


/// \brief Getter for whether to output the hits starts/stops *after* trimming
const hit_boundary_output & crh_single_output_spec::get_boundary_output() const {
	return boundary_output;
}

/// \brief Getter for whether to output HTML describing the hits and the results
const bool & crh_single_output_spec::get_generate_html_output() const {
	return generate_html_output;
}

/// \brief Getter for whether to output the results in JSON format
const bool & crh_single_output_spec::get_json_output() const {
	return json_output;
}

/// \brief Getter for whether to output a summary of the input data
const bool & crh_single_output_spec::get_summarise() const {
	return summarise;
}

/// \brief Getter for an optional file to which the cath-resolve-hits CSS should be dumped
const path_opt & crh_single_output_spec::get_export_css_file() const {
	return export_css_file;
}

/// \brief Getter for whether to output a summary of the hmmsearch output alignment
const bool & crh_single_output_spec::get_output_hmmsearch_aln() const {
	return output_hmmsearch_aln;
}

/// \brief Setter for the output file to which data should be written
crh_single_output_spec & crh_single_output_spec::set_output_file(const path &arg_output_file ///< The output file to which data should be written
                                                                 ) {
	output_file = arg_output_file;
	return *this;
}

/// \brief Setter for whether to output the hits starts/stops *after* trimming
crh_single_output_spec & crh_single_output_spec::set_boundary_output(const hit_boundary_output &arg_boundary_output ///< Whether to output the hits starts/stops *after* trimming
                                                                     ) {
	boundary_output = arg_boundary_output;
	return *this;
}

/// \brief Setter for whether to output HTML describing the hits and the results
crh_single_output_spec & crh_single_output_spec::set_generate_html_output(const bool &arg_generate_html_output ///< Whether to output HTML describing the hits and the results
                                                                          ) {
	generate_html_output = arg_generate_html_output;
	return *this;
}

/// \brief Setter for whether to output the results in JSON format
crh_single_output_spec & crh_single_output_spec::set_json_output(const bool &arg_json_output ///< Whether to output the results in JSON format
                                                                 ) {
	json_output = arg_json_output;
	return *this;
}

/// \brief Setter for whether to output a summary of the input data
crh_single_output_spec & crh_single_output_spec::set_summarise(const bool &arg_summarise ///< Whether to output a summary of the input data
                                                               ) {
	summarise = arg_summarise;
	return *this;
}

/// \brief Setter for an optional file to which the cath-resolve-hits CSS should be dumped
crh_single_output_spec & crh_single_output_spec::set_export_css_file(const path_opt &arg_export_css_file ///< An optional file to which the cath-resolve-hits CSS should be dumped
                                                                     ) {
	export_css_file = arg_export_css_file;
	return  *this;
}

/// \brief Setter for whether to output a summary of the hmmsearch output alignment
crh_single_output_spec & crh_single_output_spec::set_output_hmmsearch_aln(const bool &arg_output_hmmsearch_aln ///< Whether to output a summary of the hmmsearch output alignment
                                                                          ) {
	output_hmmsearch_aln = arg_output_hmmsearch_aln;
	return *this;
}

/// \brief Get the crh_out_format implied by the specified crh_single_output_spec
///
/// \pre There is a single, unambiguous output format
///
/// \relates crh_single_output_spec
crh_out_format cath::rslv::get_out_format(const crh_single_output_spec &arg_single_output_spec ///< The crh_single_output_spec to query
                                          ) {
	if ( arg_single_output_spec.get_summarise() ) {
		return crh_out_format::SUMMARY;
	}
	if ( arg_single_output_spec.get_generate_html_output() ) {
		return crh_out_format::HTML;
	}
	if ( arg_single_output_spec.get_json_output() ) {
		return crh_out_format::JSON;
	}
	return crh_out_format::STANDARD;
}

/// \brief Generate a description of any problem that makes the specified crh_single_output_spec invalid
///        or none otherwise
///
/// \relates crh_single_output_spec
str_opt cath::rslv::get_invalid_description(const crh_single_output_spec &arg_single_output_spec ///< The crh_single_output_spec to query
                                            ) {
	// Prepare a list of the mutually exclusive outputs that have been requested
	// Note: Don't include CSS export here because that takes a file argument and is independent of the output format
	const auto  mut_excl_output_opts = crh_out_format_opt_vec{
		make_optional( arg_single_output_spec.get_summarise(),            crh_out_format::SUMMARY ),
		make_optional( arg_single_output_spec.get_generate_html_output(), crh_out_format::HTML    ),
		make_optional( arg_single_output_spec.get_json_output(),          crh_out_format::JSON    )
	};
	const auto  mut_excl_outputs = transform_build<crh_out_format_vec>(
		mut_excl_output_opts | filtered( [] (const crh_out_format_opt &x) { return static_cast<bool>( x ); } ),
		[] (const crh_out_format_opt &x) { return *x; }
	);

	// Ensure that at most one output format has been specified
	if ( mut_excl_outputs.size() > 1 ) {
		return "Can only specify one output format but requested: " + join(
			mut_excl_outputs | transformed( to_string ),
			", "
		);
	}

	const auto out_format     = get_out_format( arg_single_output_spec );
	const auto out_format_str = to_string     ( out_format      );

	if ( means_output_trimmed_hits( arg_single_output_spec.get_boundary_output() ) && out_format != crh_out_format::STANDARD ) {
		return
			"Cannot specify trimmed boundaries if output format is "
			+ out_format_str
			+ " (but the output format may already contain the information you require)";
	}

	return none;
}

/// \brief Convenience setter for whether to output the hits starts/stops *after* trimming
///
/// \relates crh_single_output_spec
crh_single_output_spec & cath::rslv::set_output_trimmed_hits(crh_single_output_spec &arg_single_output_spec, ///< The crh_single_output_spec to modify
                                                             const bool             &arg_output_trimmed_hits ///< Whether to output the hits starts/stops *after* trimming
                                                             ) {
	return arg_single_output_spec.set_boundary_output( hit_boundary_output_of_output_trimmed_hits( arg_output_trimmed_hits ) );
}
