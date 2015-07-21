/// \file
/// \brief The residue_name class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "residue_name.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "exception/invalid_argument_exception.h"

using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;
using boost::optional;

/// \brief Throw if insert code is invalid
void residue_name::sanity_check() const {
	if ( insert_code ) {
		if ( ! is_alnum()( *insert_code ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Residue name's insert code '"
				+ string(1, *insert_code)
				+ "' is not a valid alphanumeric character"
			));
		}
	}
}

/// \brief TODOCUMENT
void residue_name::sanity_check_is_not_null_residue() const {
	if ( get_is_null_residue_name() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot access number or insert code of null residue"));
	}
}

/// \brief Ctor for residue_name
residue_name::residue_name() {
	sanity_check();
}

/// \brief Ctor for residue_name
residue_name::residue_name(const int &arg_residue_number ///< TODOCUMENT
                           ) : is_null_residue_name( false              ),
                               residue_number      ( arg_residue_number ) {
	sanity_check();
}

/// \brief Ctor for residue_name
residue_name::residue_name(const int  &arg_residue_number, ///< TODOCUMENT
                           const char &arg_insert_code     ///< TODOCUMENT
                           ) : is_null_residue_name( false              ),
                               residue_number      ( arg_residue_number ),
                               insert_code         ( arg_insert_code    ) {
	sanity_check();
}

/// \brief TODOCUMENT
const bool & residue_name::get_is_null_residue_name() const {
	return is_null_residue_name;
}

/// \brief TODOCUMENT
const int & residue_name::get_residue_number() const {
	sanity_check_is_not_null_residue();
	return residue_number;
}

/// \brief TODOCUMENT
const optional<char> & residue_name::get_opt_insert_code() const {
	sanity_check_is_not_null_residue();
	return insert_code;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
bool cath::operator==(const residue_name &arg_residue_name_a, ///< TODOCUMENT
                      const residue_name &arg_residue_name_b  ///< TODOCUMENT
					  ) {
	const bool a_is_null = arg_residue_name_a.get_is_null_residue_name();
	const bool b_is_null = arg_residue_name_b.get_is_null_residue_name();
	if ( a_is_null                                != b_is_null                                ) {
		return false;
	}
	if ( ! a_is_null ) {
		if ( arg_residue_name_a.get_residue_number()  != arg_residue_name_b.get_residue_number()  ) {
			return false;
		}
		if ( arg_residue_name_a.get_opt_insert_code() != arg_residue_name_b.get_opt_insert_code() ) {
			return false;
		}
	}
	return true;
}

/// \brief Simple insertion operator for residue_name
///
/// \relates residue_name
ostream & cath::operator<<(ostream            &arg_os,          ///< The ostream to which the residue_name should be output
                           const residue_name &arg_residue_name ///< The residue_name to output
                           ) {
	if ( arg_residue_name.get_is_null_residue_name() ) {
		arg_os << "null_res";
	}
	else {
		arg_os << arg_residue_name.get_residue_number();
		if ( has_insert_code( arg_residue_name ) ) {
			arg_os << get_insert_code( arg_residue_name );
		}
	}
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
istream & cath::operator>>(istream      &arg_istream,     ///< TODOCUMENT
                           residue_name &arg_residue_name ///< TODOCUMENT
						   ) {
	string input_string;
	arg_istream >> input_string;
	arg_residue_name = make_residue_name( input_string );
	return arg_istream;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
bool cath::has_insert_code(const residue_name &arg_residue_name ///< TODOCUMENT
                           ) {
	return static_cast<bool>( arg_residue_name.get_opt_insert_code() );
}

/// \brief TODOCUMENT
///
/// \relates residue_name
char cath::get_insert_code(const residue_name &arg_residue_name ///< TODOCUMENT
                           ) {
	return *arg_residue_name.get_opt_insert_code();
}

/// \brief TODOCUMENT
///
/// \relates residue_name
string cath::get_insert_code_string(const residue_name &arg_residue_name ///< TODOCUMENT
                                    ) {
	return has_insert_code( arg_residue_name ) ? string( 1, get_insert_code( arg_residue_name ) ) : "";
}

/// \brief TODOCUMENT
///
/// \relates residue_name
string cath::make_residue_name_string_with_insert_or_space(const residue_name &arg_residue_name ///< TODOCUMENT
                                                           ) {
	const string insert_or_space = has_insert_code( arg_residue_name ) ? string( 1, get_insert_code( arg_residue_name ) ) : " ";
	return lexical_cast<string>( arg_residue_name.get_residue_number() ) + insert_or_space;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
residue_name cath::make_residue_name(const string &arg_residue_name ///< TODOCUMENT
                                     ) {
	// Sanity check the input
	if ( arg_residue_name.empty() ) {
		return residue_name();
	}
	if ( arg_residue_name.length() == 1 && ! is_digit()( arg_residue_name.at( 0 ) ) && arg_residue_name.at( 0 ) != '-' ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot split residue name into number and insert if it only contains an insert code"));
	}

	// Grab the last character in the string and check whether it's a numeric digit or an insert code
	const size_t residue_name_length = arg_residue_name.length();
	const char   res_name_last_char  = arg_residue_name.at( residue_name_length - 1 );
	const bool   res_name_has_insert = ! is_digit()( res_name_last_char );

	// If the last character's an insert code then make a residue_namethere isn't any insert code then
	if ( res_name_has_insert ) {
		const int residue_number = lexical_cast<int>( arg_residue_name.substr(0, ( residue_name_length - 1 ) ) );
		return residue_name( residue_number, res_name_last_char );
	}
	// Otherwise, no insert code so just return a residue_name of the string converted to an int
	else {
		return residue_name( lexical_cast<int>( arg_residue_name ) );
	}
}

/// \brief TODOCUMENT
///
/// \relates residue_name
residue_name cath::make_residue_name_with_non_insert_char(const int  &arg_residue_number,       ///< TODOCUMENT
                                                          const char &arg_possible_insert_code, ///< TODOCUMENT
                                                          const char &arg_non_insert_char       ///< TODOCUMENT
                                                          ) {
	return ( arg_possible_insert_code == arg_non_insert_char ) ? residue_name( arg_residue_number                           )
	                                                           : residue_name( arg_residue_number, arg_possible_insert_code );
}
