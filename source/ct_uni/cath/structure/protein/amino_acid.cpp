/// \file
/// \brief The amino_acid class definitions

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

#include "amino_acid.hpp"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"

#include <set>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::std;

/// \brief Generate a string describing the specified amino_acid_type
///
/// \relates amino_acid_type
string cath::to_string(const amino_acid_type &prm_amino_acid_type ///< The amino_acid_type to describe
                       ) {
	switch ( prm_amino_acid_type ) {
		case ( amino_acid_type::AA      ) : { return "AA"      ; };
		case ( amino_acid_type::HETATOM ) : { return "HETATOM" ; };
		case ( amino_acid_type::DNA     ) : { return "DNA"     ; };
	}
	BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of amino_acid_type not recognised whilst converting to_string()"));
}

/// \brief Make a vector of amino acids from a vector of amino acid chars
///
/// \relates amino_acid
amino_acid_vec cath::make_amino_acids_of_chars(const char_vec &prm_amino_acid_chars ///< A vector of chars for amino acids
                                               ) {
	return transform_build<amino_acid_vec>(
		prm_amino_acid_chars,
		[] (const char &x) { return amino_acid{ x }; }
	);
}

/// \brief Get the three-letter-code char_3_arr associated with the specified one letter
///
/// eg 'A' -> "ALA"
///
/// \relates amino_acid
char_3_arr cath::get_code_of_amino_acid_letter(const char &prm_one_letter_aa ///< The single-letter amino acid (eg 'A' for alanine)
                                               ) {
	return amino_acid( prm_one_letter_aa ).get_code();
}

/// \brief Get the three-letter-code string associated with the specified one letter
///
/// eg 'A' -> "ALA"
///
/// \relates amino_acid
string cath::get_code_str_of_amino_acid_letter(const char &prm_one_letter_aa ///< The single-letter amino acid (eg 'A' for alanine)
                                               ) {
	return string_of_char_arr( get_code_of_amino_acid_letter( prm_one_letter_aa ) );
}

/// \brief TODOCUMENT
///
/// \relates amino_acid
char cath::get_letter_of_amino_acid_code(const string_view &prm_three_letter_aa ///< TODOCUMENT
                                         ) {
	return *( amino_acid( prm_three_letter_aa ).get_letter_if_amino_acid() );
}

/// \brief Insert a description of the specified amino_acid into the specified ostream
///
/// \relates amino_acid
ostream & cath::operator<<(ostream          &prm_os,        ///< The ostream into which the description should be inserted
                           const amino_acid &prm_amino_acid ///< The amino_acid to describe
                           ) {
	prm_os << get_code_string( prm_amino_acid );
	return prm_os;
}

/// \brief TODOCUMENT
///
/// \relates amino_acid
istream & cath::operator>>(istream    &is,            ///< The istream from which to parse the amino_acid
                           amino_acid &prm_amino_acid ///< The amino acid to populate
                           ) {
	string input_string;
	is >> input_string;

	try {
		prm_amino_acid = amino_acid(input_string);
	}
	catch (const invalid_argument_exception &) {
		BOOST_THROW_EXCEPTION(invalid_argument("invalid_evaluator_argument"));
		exit(1);
	}
	return is;
}
