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

#include "crh_output_spec.h"

using namespace cath;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::filesystem::path;
using boost::make_optional;

constexpr hit_boundary_output crh_output_spec::DEFAULT_BOUNDARY_OUTPUT;
constexpr bool                crh_output_spec::DEFAULT_GENERATE_HTML_OUTPUT;
constexpr bool                crh_output_spec::DEFAULT_OUTPUT_HMMSEARCH_ALN;

/// \brief Getter for the output file to which data should be written
const path_opt & crh_output_spec::get_output_file() const {
	return output_file;
}


/// \brief Getter for whether to output the hits starts/stops *after* trimming
const hit_boundary_output & crh_output_spec::get_boundary_output() const {
	return boundary_output;
}

/// \brief Getter for whether to output HTML describing the hits and the results
const bool & crh_output_spec::get_generate_html_output() const {
	return generate_html_output;
}

/// \brief Getter for whether to output a summary of the hmmsearch output alignment
const bool & crh_output_spec::get_output_hmmsearch_aln() const {
	return output_hmmsearch_aln;
}

/// \brief Setter for the output file to which data should be written
crh_output_spec & crh_output_spec::set_output_file(const path &arg_output_file ///< The output file to which data should be written
                                                   ) {
	output_file = arg_output_file;
	return *this;
}

/// \brief Setter for whether to output the hits starts/stops *after* trimming
crh_output_spec & crh_output_spec::set_boundary_output(const hit_boundary_output &arg_boundary_output ///< Whether to output the hits starts/stops *after* trimming
                                                       ) {
	boundary_output = arg_boundary_output;
	return *this;
}

/// \brief Setter for whether to output HTML describing the hits and the results
crh_output_spec & crh_output_spec::set_generate_html_output(const bool &arg_generate_html_output ///< Whether to output HTML describing the hits and the results
                                                            ) {
	generate_html_output = arg_generate_html_output;
	return *this;
}

/// \brief Setter for whether to output a summary of the hmmsearch output alignment
crh_output_spec & crh_output_spec::set_output_hmmsearch_aln(const bool &arg_output_hmmsearch_aln ///< Whether to output a summary of the hmmsearch output alignment
                                                            ) {
	output_hmmsearch_aln = arg_output_hmmsearch_aln;
	return *this;
}

/// \brief Convenience setter for whether to output the hits starts/stops *after* trimming
///
/// \relates crh_output_spec
crh_output_spec & cath::rslv::set_output_trimmed_hits(crh_output_spec &arg_output_spec,        ///< The crh_output_spec to modify
                                                      const bool      &arg_output_trimmed_hits ///< Whether to output the hits starts/stops *after* trimming
                                                      ) {
	return arg_output_spec.set_boundary_output( hit_boundary_output_of_output_trimmed_hits( arg_output_trimmed_hits ) );
}

/// \brief Generate a description of any problem that makes the specified crh_output_spec invalid
///        or none otherwise
///
/// \relates crh_output_spec
str_opt cath::rslv::get_invalid_description(const crh_output_spec &arg_output_spec ///< The crh_output_spec to query
                                            ) {
	return make_optional(
		means_output_trimmed_hits( arg_output_spec.get_boundary_output() ) && arg_output_spec.get_generate_html_output(),
		"Cannot output trimmed boundaries in HTML output"s
	);
}