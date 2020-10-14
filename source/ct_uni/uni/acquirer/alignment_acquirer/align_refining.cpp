/// \file
/// \brief The align_refining class definitions

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

#include "align_refining.hpp"

#include <boost/algorithm/string/case_conv.hpp>

#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "common/program_options/validator.hpp"

using namespace cath::common;

using boost::any;
using boost::to_upper;
using std::istream;
using std::ostream;
using std::string;

/// \brief Generate a string describing the specified align_refining
///
/// \relates align_refining
string cath::align::to_string(const align_refining &prm_align_refining ///< The align_refining to describe
                              ) {
	switch ( prm_align_refining ) {
		case ( align_refining::NO    ) : { return "NO"    ; }
		case ( align_refining::LIGHT ) : { return "LIGHT" ; }
		case ( align_refining::HEAVY ) : { return "HEAVY" ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("align_refining value not recognised in to_string()"));
}

/// \brief Insert a description of the specified align_refining into the specified ostream
///
/// \relates align_refining
ostream & cath::align::operator<<(ostream              &prm_os,            ///< The ostream into which the description should be inserted
                                  const align_refining &prm_align_refining ///< The align_refining to describe
                                  ) {
	prm_os << to_string( prm_align_refining );
	return prm_os;
}

/// \brief Simple extraction operator for align_refining
///
/// \relates align_refining
istream & cath::align::operator>>(istream        &prm_is,            ///< The istream from which to extract the align_refining
                                  align_refining &prm_align_refining ///< The align_refining to populate
                                  ) {
	string input_string;
	prm_is >> input_string;
	to_upper( input_string );
	if (     input_string == "NO"     ) {
		prm_align_refining = align_refining::NO;
	}
	else if ( input_string == "LIGHT" ) {
		prm_align_refining = align_refining::LIGHT;
	}
	else if ( input_string == "HEAVY" ) {
		prm_align_refining = align_refining::HEAVY;
	}
	else {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse align_refining from string \""+ input_string + "\"" ));
	}
	return prm_is;
}

/// \brief Generate a string containing a description of the specified align_refining
///
/// \relates align_refining
string cath::align::description_of_align_refining(const align_refining &prm_align_refining ///< The align_refining to describe
                                                  ) {
	switch ( prm_align_refining ) {
		case ( align_refining::NO    ) : { return "Don't refine the alignment" ; }
		case ( align_refining::LIGHT ) : { return "Refine any alignments with few entries; glue alignments one more entry at a time" ; }
		case ( align_refining::HEAVY ) : { return "Perform heavy (slow) refining on the alignment, including when gluing alignments together" ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("align_refining value not recognised in description_of_align_refining()"));
}

/// \brief Provide Boost program_options validation for align_refining
///
/// \relates align_refining
void cath::align::validate(any           &prm_value,         ///< The value to populate
                           const str_vec &prm_value_strings, ///< The string values to validate
                           align_refining *, int) {
	prm_value = lex_castable_validator<align_refining>::perform_validate( prm_value, prm_value_strings );
}
