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

#include "residue_name.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;

using ::boost::algorithm::is_alnum;
using ::boost::algorithm::is_digit;
using ::boost::lexical_cast;
using ::boost::optional;
using ::std::istream;
using ::std::ostream;
using ::std::string;

/// \brief Throw if insert code is invalid
void residue_name::sanity_check() const {
	if ( insert ) {
		if ( ! is_alnum()( *insert ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Residue name's insert code '"
				+ string{ *insert }
				+ "' is not a valid alphanumeric character"
			));
		}
	}
}

/// \brief TODOCUMENT
void residue_name::sanity_check_is_not_null_residue() const {
	if ( is_null() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot access number or insert code of null residue"));
	}
}

/// \brief Ctor for residue_name
residue_name::residue_name() {
	sanity_check();
}

/// \brief Ctor for residue_name
residue_name::residue_name(const int &prm_residue_number ///< TODOCUMENT
                           ) : res_num             ( prm_residue_number ),
                               is_null_residue_name( false              ) {
	sanity_check();
}

/// \brief Ctor for residue_name
residue_name::residue_name(const int  &prm_residue_number, ///< TODOCUMENT
                           const char &prm_insert     ///< TODOCUMENT
                           ) : res_num             ( prm_residue_number ),
                               insert              ( prm_insert         ),
                               is_null_residue_name( false              ) {
	sanity_check();
}

/// \brief TODOCUMENT
const bool & residue_name::is_null() const {
	return is_null_residue_name;
}

/// \brief TODOCUMENT
const int & residue_name::residue_number() const {
	sanity_check_is_not_null_residue();
	return res_num;
}

/// \brief TODOCUMENT
const optional<char> & residue_name::opt_insert() const {
	sanity_check_is_not_null_residue();
	return insert;
}

/// \brief TODOCUMENT
///
/// \relates residue_name
bool cath::operator==(const residue_name &prm_residue_name_a, ///< TODOCUMENT
                      const residue_name &prm_residue_name_b  ///< TODOCUMENT
                      ) {
	return (
		( prm_residue_name_a.is_null() == prm_residue_name_b.is_null() )
		&&
		(
			prm_residue_name_a.is_null()
			||
			(
				prm_residue_name_a.residue_number()  == prm_residue_name_b.residue_number()
				&&
				prm_residue_name_a.opt_insert() == prm_residue_name_b.opt_insert()
			)
		)
	);
}


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

/// \brief Get the specified residue_name's number or the specified value if the residue_name is null
///
/// \relates residue_name
int cath::residue_number_or_value_if_null(const residue_name &prm_residue_name, ///< The residue_name to query
                                          const int          &prm_value         ///< The value to use if the residue_name is null
                                          ) {
	return prm_residue_name.is_null() ? prm_value : prm_residue_name.residue_number();
}

/// \brief Get the specified residue_name's insert code or the specified value if the residue_name is null
///
/// \relates residue_name
optional<char> cath::opt_insert_or_value_if_null(const residue_name   &prm_residue_name, ///< The residue_name to query
                                                 const optional<char> &prm_value         ///< The value to use if the residue_name is null
                                                 ) {
	return prm_residue_name.is_null() ? prm_value : prm_residue_name.opt_insert();
}

/// \brief Get whether the specified residue_name has an insert code
///
/// \relates residue_name
bool cath::has_insert(const residue_name &prm_residue_name ///< The residue_name to query
                      ) {
	return static_cast<bool>( prm_residue_name.opt_insert() );
}

/// \brief Get whether the specified residue_name has an insert code or the specified value if the residue_name is null
///
/// \relates residue_name
bool cath::has_insert_or_value_if_null(const residue_name &prm_residue_name, ///< The residue_name to query
                                       const bool         &prm_value         ///< The value to use if the residue_name is null
                                       ) {
	return prm_residue_name.is_null() ? prm_value : has_insert( prm_residue_name );
}


/// \brief Get the specified residue_name's insert code
///
/// \relates residue_name
const char & cath::insert(const residue_name &prm_residue_name ///< The residue_name to query
                          ) {
	return *prm_residue_name.opt_insert();
}

/// \brief Get the specified residue_name's insert code or the specified value if the residue_name is null
///
/// \relates residue_name
char cath::insert_or_value_if_null(const residue_name &prm_residue_name, ///< The residue_name to query
                                   const char         &prm_value         ///< The value to use if the residue_name is null
                                   ) {
	return prm_residue_name.is_null() ? prm_value : insert( prm_residue_name );
}

/// \brief Get the specified residue_name's insert code or the specified value if the residue_name is null or if it has no insert code
///
/// \relates residue_name
char cath::insert_or_value_if_null_or_absent(const residue_name &prm_residue_name, ///< The residue_name to query
                                             const char         &prm_value         ///< The value to use if the residue_name is null or the insert is absent
                                             ) {
	return ( prm_residue_name.is_null() || ! has_insert( prm_residue_name ) )
		? prm_value
		: insert( prm_residue_name );
}

/// \brief Return whether the specified residue_name has a strictly negative residue number
///
/// \relates residue_name
bool cath::has_strictly_negative_residue_number(const residue_name &prm_residue_name ///< The residue_name to query
                                                ) {
	return ( ! prm_residue_name.is_null() && prm_residue_name.residue_number() < 0 );
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
		return residue_name();
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

/// \brief TODOCUMENT
///
/// \relates residue_name
residue_name cath::make_residue_name_with_non_insert_char(const int  &prm_residue_number,  ///< TODOCUMENT
                                                          const char &prm_possible_insert, ///< TODOCUMENT
                                                          const char &prm_non_insert_char  ///< TODOCUMENT
                                                          ) {
	return ( prm_possible_insert == prm_non_insert_char ) ? residue_name( prm_residue_number                      )
	                                                      : residue_name( prm_residue_number, prm_possible_insert );
}
