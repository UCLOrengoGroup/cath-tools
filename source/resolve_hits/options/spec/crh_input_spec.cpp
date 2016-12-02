/// \file
/// \brief The crh_input_spec class definitions

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

#include "crh_input_spec.hpp"

using namespace cath;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::filesystem::path;
using boost::none;

constexpr bool                  crh_input_spec::DEFAULT_READ_FROM_STDIN;
constexpr hits_input_format_tag crh_input_spec::DEFAULT_INPUT_FORMAT;
constexpr residx_t              crh_input_spec::DEFAULT_MIN_GAP_LENGTH;
constexpr bool                  crh_input_spec::DEFAULT_INPUT_HITS_ARE_GROUPED;

/// \brief Getter for the input file from which data should be read
const path_opt & crh_input_spec::get_input_file() const {
	return input_file;
}

/// \brief Getter for whether to read the input data from stdin
const bool & crh_input_spec::get_read_from_stdin() const {
	return read_from_stdin;
}

/// \brief Getter for the format of the input data
const hits_input_format_tag & crh_input_spec::get_input_format() const {
	return input_format;
}

/// \brief Getter for the minimum gap length to consider when parsing an alignment
const residx_t & crh_input_spec::get_min_gap_length() const {
	return min_gap_length;
}

/// \brief Getter for whether the code can assume that the input data is pre-grouped by query_id
const bool & crh_input_spec::get_input_hits_are_grouped() const {
	return input_hits_are_grouped;
}

/// \brief Setter for the input file from which data should be read
crh_input_spec & crh_input_spec::set_input_file(const path &arg_input_file ///< The input file from which data should be read
                                                ) {
	input_file = arg_input_file;
	return *this;
}

/// \brief Setter for whether to read the input data from stdin
crh_input_spec & crh_input_spec::set_read_from_stdin(const bool &arg_read_from_stdin ///< Whether to read the input data from stdin
                                                     ) {
	read_from_stdin = arg_read_from_stdin;
	return *this;
}

/// \brief Setter for the format of the input data
crh_input_spec & crh_input_spec::set_input_format(const hits_input_format_tag &arg_input_format ///< The format of the input data
                                                  ) {
	input_format = arg_input_format;
	return *this;
}

/// \brief Setter for the minimum gap length to consider when parsing an alignment
crh_input_spec & crh_input_spec::set_min_gap_length(const residx_t &arg_min_gap_length ///< The minimum gap length to consider when parsing an alignment
                                                    ) {
	min_gap_length = arg_min_gap_length;
	return *this;
}

/// \brief Setter for whether the code can assume that the input data is pre-grouped by query_id
crh_input_spec & crh_input_spec::set_input_hits_are_grouped(const bool &arg_input_hits_are_grouped ///< Whether the code can assume that the input data is pre-grouped by query_id
                                                            ) {
	input_hits_are_grouped = arg_input_hits_are_grouped;
	return *this;
}

/// \brief Generate a description of any problem that makes the specified crh_input_spec invalid
///        or none otherwise
///
/// \relates crh_input_spec
str_opt cath::rslv::get_invalid_description(const crh_input_spec &arg_spec ///< The crh_input_spec to query
                                            ) {
	if ( arg_spec.get_input_file() && arg_spec.get_read_from_stdin() ) {
		return "Cannot read from both a file and stdin"s;
	}

	return none;
}
