/// \file
/// \brief The resolve_hits_input_spec class definitions

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

#include "resolve_hits_input_spec.h"

using namespace cath::opts;
using namespace cath;
using namespace std::literals::string_literals;

using boost::none;

constexpr bool resolve_hits_input_spec::DEFAULT_READ_FROM_STDIN;
constexpr bool resolve_hits_input_spec::DEFAULT_INPUT_HITS_ARE_GROUPED;
constexpr bool resolve_hits_input_spec::DEFAULT_DOMTBLOUT;
constexpr bool resolve_hits_input_spec::DEFAULT_HMMERALN;
constexpr bool resolve_hits_input_spec::DEFAULT_RAW_SCORE_IS_EVALUE;
constexpr bool resolve_hits_input_spec::DEFAULT_APPLY_CATH_RULES;

/// \brief Getter for TODOCUMENT
const opt_path & resolve_hits_input_spec::get_input_file() const {
	return input_file;
}

/// \brief Getter for TODOCUMENT
const bool & resolve_hits_input_spec::get_read_from_stdin() const {
	return read_from_stdin;
}

/// \brief Getter for TODOCUMENT
const bool & resolve_hits_input_spec::get_input_hits_are_grouped() const {
	return input_hits_are_grouped;
}

/// \brief Getter for TODOCUMENT
const bool & resolve_hits_input_spec::get_domtblout() const {
	return domtblout;
}

/// \brief Getter for TODOCUMENT
const bool & resolve_hits_input_spec::get_hmmeraln() const {
	return hmmeraln;
}

/// \brief Getter for TODOCUMENT
const bool & resolve_hits_input_spec::get_raw_score_is_evalue() const {
	return raw_score_is_evalue;
}

/// \brief Getter for TODOCUMENT
const bool & resolve_hits_input_spec::get_apply_cath_rules() const {
	return apply_cath_rules;
}

/// \brief Setter for TODOCUMENT
resolve_hits_input_spec & resolve_hits_input_spec::set_input_file(const boost::filesystem::path &arg_input_file ///< TODOCUMENT
                                                                  ) {
	input_file = arg_input_file;
	return *this;
}

/// \brief Setter for TODOCUMENT
resolve_hits_input_spec & resolve_hits_input_spec::set_read_from_stdin(const bool &arg_read_from_stdin ///< TODOCUMENT
                                                                       ) {
	read_from_stdin = arg_read_from_stdin;
	return *this;
}

/// \brief Setter for TODOCUMENT
resolve_hits_input_spec & resolve_hits_input_spec::set_input_hits_are_grouped(const bool &arg_input_hits_are_grouped ///< TODOCUMENT
                                                                             ) {
	input_hits_are_grouped = arg_input_hits_are_grouped;
	return *this;
}

/// \brief Setter for TODOCUMENT
resolve_hits_input_spec & resolve_hits_input_spec::set_domtblout(const bool &arg_domtblout ///< TODOCUMENT
                                                                 ) {
	domtblout = arg_domtblout;
	return *this;
}

/// \brief Setter for TODOCUMENT
resolve_hits_input_spec & resolve_hits_input_spec::set_hmmeraln(const bool &arg_hmmeraln ///< TODOCUMENT
                                                                ) {
	hmmeraln = arg_hmmeraln;
	return *this;
}

/// \brief Setter for TODOCUMENT
resolve_hits_input_spec & resolve_hits_input_spec::set_raw_score_is_evalue(const bool &arg_raw_score_is_evalue ///< TODOCUMENT
                                                                           ) {
	raw_score_is_evalue = arg_raw_score_is_evalue;
	return *this;
}

/// \brief Setter for TODOCUMENT
resolve_hits_input_spec & resolve_hits_input_spec::set_apply_cath_rules(const bool &arg_apply_cath_rules ///< TODOCUMENT
                                                                        ) {
	apply_cath_rules = arg_apply_cath_rules;
	return *this;
}

/// \brief TODOCUMENT
///
/// \relates resolve_hits_input_spec
opt_str cath::opts::get_invalid_description(const resolve_hits_input_spec &arg_spec ///< TODOCUMENT
                                            ) {
	if ( arg_spec.get_domtblout() && arg_spec.get_hmmeraln() ) {
		return "Format cannot be both domtblout and hmmeraln"s;
	}
	if ( arg_spec.get_input_file() && arg_spec.get_read_from_stdin() ) {
		return "Cannot read from both a file and stdin"s;
	}
	// if ( ! arg_spec.get_input_file() && ! arg_spec.get_read_from_stdin() ) {
	// 	return "Must read from either a file or stdin"s;
	// }
	if ( arg_spec.get_raw_score_is_evalue() && ( arg_spec.get_domtblout() || arg_spec.get_hmmeraln() ) ) {
		return "Cannot treat a raw score as an evalue for domtblout or hmmeraln input"s;
	}

	return none;
}
