/// \file
/// \brief The residue_name class definitions

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

#include <cstddef>
#include <istream>
#include <ostream>
#include <string>

#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>

#include "cath/biocore/residue_name.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath;
using namespace ::cath::common;

using ::boost::algorithm::is_digit;
using ::boost::lexical_cast;
using ::std::istream;
using ::std::ostream;
using ::std::string;

/// \brief Simple to_string() overload for residue_name
///
/// \relates residue_name
string cath::to_string(const residue_name &prm_residue_name ///< The residue_name to be output as a string
                       ) {
	if ( prm_residue_name.is_null() ) {
		return "null_res";
	}
	const string insert_string = has_insert( prm_residue_name )
	                             ? string{ insert( prm_residue_name ) }
	                             : string{};
	return ::std::to_string( prm_residue_name.residue_number() ) + insert_string;
}

/// \brief Simple insertion operator for residue_name
///
/// \relates residue_name
ostream & cath::operator<<(ostream            &prm_os,          ///< The ostream to which the residue_name should be output
                           const residue_name &prm_residue_name ///< The residue_name to output
                           ) {
	prm_os << to_string( prm_residue_name );
	return prm_os;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
istream & cath::operator>>(istream      &prm_istream,     ///< TODOCUMENT
                           residue_name &prm_residue_name ///< TODOCUMENT
                           ) {
	string input_string;
	prm_istream >> input_string;
	prm_residue_name = make_residue_name( input_string );
	return prm_istream;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
string cath::insert_string(const residue_name &prm_residue_name ///< The residue_name to query
                           ) {
	return has_insert( prm_residue_name ) ? string{ insert( prm_residue_name ) } : "";
}

/// \brief TODOCUMENT
///
/// \relates residue_name
string cath::make_residue_name_string_with_insert_or_space(const residue_name &prm_residue_name ///< TODOCUMENT
                                                           ) {
	const string insert_or_space = has_insert( prm_residue_name ) ? string{ insert( prm_residue_name ) } : " ";
	return lexical_cast<string>( prm_residue_name.residue_number() ) + insert_or_space;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
residue_name cath::make_residue_name(const string &prm_residue_name ///< TODOCUMENT
                                     ) {
	// Sanity check the input
	if ( prm_residue_name.empty() ) {
		return {};
	}
	if ( prm_residue_name.length() == 1 && ! is_digit()( prm_residue_name.at( 0 ) ) && prm_residue_name.at( 0 ) != '-' ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot split residue name into number and insert if it only contains an insert code"));
	}

	// Grab the last character in the string and check whether it's a numeric digit or an insert code
	const size_t residue_name_length = prm_residue_name.length();
	const char   res_name_last_char  = prm_residue_name.at( residue_name_length - 1 );
	const bool   res_name_has_insert = ! is_digit()( res_name_last_char );

	// If the last character's an insert code then make a residue_namethere isn't any insert code then
	if ( res_name_has_insert ) {
		const int residue_number = lexical_cast<int>( prm_residue_name.substr(0, ( residue_name_length - 1 ) ) );
		return residue_name( residue_number, res_name_last_char );
	}
	// Otherwise, no insert code so just return a residue_name of the string converted to an int
	return residue_name( lexical_cast<int>( prm_residue_name ) );
}
