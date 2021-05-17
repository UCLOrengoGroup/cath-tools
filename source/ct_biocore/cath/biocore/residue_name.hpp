/// \file
/// \brief The residue_name class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_NAME_HPP
#define _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_NAME_HPP

#include <optional>

#include <boost/operators.hpp>

#include <fmt/core.h>

#include "cath/common/exception/invalid_argument_exception.hpp"

namespace cath {

	/// \brief Represent a PDB residue name (eg 324A)
	class residue_name final : private boost::equality_comparable<residue_name> {
	  private:
		/// \brief The residue number
		int res_num = 0;

		/// \brief The (optional insert code)
		::std::optional<char> insert;

		/// \brief Whether this is a null residue
		bool is_null_residue_name = true;

		constexpr void sanity_check() const;
		constexpr void sanity_check_is_not_null_residue() const;

	  public:
		constexpr residue_name();
		constexpr explicit residue_name( const int & );
		constexpr residue_name( const int &, const char & );

		[[nodiscard]] constexpr const bool &is_null() const;

		[[nodiscard]] constexpr const int &residue_number() const;

		[[nodiscard]] constexpr const ::std::optional<char> &opt_insert() const;

		/// \brief TODOCUMENT
		///
		/// \relates residue_name
		///
		/// \param prm_lhs TODOCUMENT
		/// \param prm_rhs TODOCUMENT
		friend constexpr bool operator==( const residue_name &prm_lhs, const residue_name &prm_rhs ) {
			return ( ( prm_lhs.is_null() == prm_rhs.is_null() )
			         && ( prm_lhs.is_null()
			              || ( prm_lhs.residue_number() == prm_rhs.residue_number()
			                   && prm_lhs.opt_insert() == prm_rhs.opt_insert() ) ) );
		}
	};

	constexpr bool constexpr_is_alnum( const char &x ) {
		return ( ( x >= '0' && x <= '9' ) || ( x >= 'A' && x <= 'Z' ) || ( x >= 'a' && x <= 'z' ) );
	}

	/// \brief Throw if insert code is invalid
	constexpr void residue_name::sanity_check() const {
		if ( insert ) {
			if ( !constexpr_is_alnum( *insert ) ) {
				BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
				  ::fmt::format( "Residue name's insert code '{}' is not a valid alphanumeric character", *insert ) ) );
			}
		}
	}

	/// \brief TODOCUMENT
	constexpr void residue_name::sanity_check_is_not_null_residue() const {
		if ( is_null() ) {
			BOOST_THROW_EXCEPTION(
			  common::invalid_argument_exception( "Cannot access number or insert code of null residue" ) );
		}
	}

	/// \brief Ctor for residue_name
	constexpr residue_name::residue_name() {
		sanity_check();
	}

	/// \brief Ctor for residue_name
	///
	/// \param prm_residue_number TODOCUMENT
	constexpr residue_name::residue_name( const int &prm_residue_number ) :
	        res_num( prm_residue_number ), is_null_residue_name( false ) {
		sanity_check();
	}

	/// \brief Ctor for residue_name
	///
	/// \param prm_residue_number TODOCUMENT
	/// \param prm_insert         TODOCUMENT
	constexpr residue_name::residue_name( const int &prm_residue_number, const char &prm_insert ) :
	        res_num( prm_residue_number ), insert( prm_insert ), is_null_residue_name( false ) {
		sanity_check();
	}

	/// \brief TODOCUMENT
	constexpr const bool &residue_name::is_null() const {
		return is_null_residue_name;
	}

	/// \brief TODOCUMENT
	constexpr const int &residue_name::residue_number() const {
		sanity_check_is_not_null_residue();
		return res_num;
	}

	/// \brief TODOCUMENT
	constexpr const ::std::optional<char> &residue_name::opt_insert() const {
		sanity_check_is_not_null_residue();
		return insert;
	}

	bool operator<( const residue_name &, const residue_name & )  = delete; // Better to avoid residue_name comparisons (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues
	bool operator<=( const residue_name &, const residue_name & ) = delete; // Better to avoid residue_name comparisons (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues
	bool operator>( const residue_name &, const residue_name & )  = delete; // Better to avoid residue_name comparisons (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues
	bool operator>=( const residue_name &, const residue_name & ) = delete; // Better to avoid residue_name comparisons (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues

	/// \brief Get the specified residue_name's number or the specified value if the residue_name is null
	///
	/// \relates residue_name
	/// \param prm_residue_name The residue_name to query
	/// \param prm_value        The value to use if the residue_name is null
	constexpr int residue_number_or_value_if_null( const residue_name &prm_residue_name, const int &prm_value ) {
		return prm_residue_name.is_null() ? prm_value : prm_residue_name.residue_number();
	}

	/// \brief Get the specified residue_name's insert code or the specified value if the residue_name is null
	///
	/// \relates residue_name
	///
	/// \param prm_residue_name The residue_name to query
	/// \param prm_value        The value to use if the residue_name is null
	constexpr ::std::optional<char> opt_insert_or_value_if_null( const residue_name &         prm_residue_name,
	                                                             const ::std::optional<char> &prm_value ) {
		return prm_residue_name.is_null() ? prm_value : prm_residue_name.opt_insert();
	}

	/// \brief Get whether the specified residue_name has an insert code
	///
	/// \relates residue_name
	///
	/// \param prm_residue_name The residue_name to query
	constexpr bool has_insert( const residue_name &prm_residue_name ) {
		return static_cast<bool>( prm_residue_name.opt_insert() );
	}

	/// \brief Get whether the specified residue_name has an insert code or the specified value if the residue_name is null
	///
	/// \relates residue_name
	///
	/// \param prm_residue_name The residue_name to query
	/// \param prm_value        The value to use if the residue_name is null
	constexpr bool has_insert_or_value_if_null( const residue_name &prm_residue_name, const bool &prm_value ) {
		return prm_residue_name.is_null() ? prm_value : has_insert( prm_residue_name );
	}

	/// \brief Get the specified residue_name's insert code
	///
	/// \relates residue_name
	///
	/// \param prm_residue_name The residue_name to query
	constexpr const char &insert( const residue_name &prm_residue_name ) {
		return *prm_residue_name.opt_insert();
	}

	/// \brief Get the specified residue_name's insert code or the specified value if the residue_name is null
	///
	/// \relates residue_name
	///
	/// \param prm_residue_name The residue_name to query
	/// \param prm_value        The value to use if the residue_name is null
	constexpr char insert_or_value_if_null( const residue_name &prm_residue_name, const char &prm_value ) {
		return prm_residue_name.is_null() ? prm_value : insert( prm_residue_name );
	}

	/// \brief Get the specified residue_name's insert code or the specified value if the residue_name is null or if it has no insert code
	///
	/// \relates residue_name
	///
	/// \param prm_residue_name The residue_name to query
	/// \param prm_value        The value to use if the residue_name is null or the insert is absent
	constexpr char insert_or_value_if_null_or_absent( const residue_name &prm_residue_name, const char &prm_value ) {
		return ( prm_residue_name.is_null() || !has_insert( prm_residue_name ) ) ? prm_value : insert( prm_residue_name );
	}

	/// \brief Return whether the specified residue_name has a strictly negative residue number
	///
	/// \relates residue_name
	///
	/// \param prm_residue_name The residue_name to query
	constexpr bool has_strictly_negative_residue_number( const residue_name &prm_residue_name ) {
		return ( !prm_residue_name.is_null() && prm_residue_name.residue_number() < 0 );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue_name
	///
	/// \param prm_residue_number  TODOCUMENT
	/// \param prm_possible_insert TODOCUMENT
	/// \param prm_non_insert_char TODOCUMENT
	constexpr residue_name make_residue_name_with_non_insert_char( const int & prm_residue_number,
	                                                               const char &prm_possible_insert,
	                                                               const char &prm_non_insert_char ) {
		return ( prm_possible_insert == prm_non_insert_char ) ? residue_name( prm_residue_number )
		                                                      : residue_name( prm_residue_number, prm_possible_insert );
	}

	std::string   to_string( const residue_name & );
	std::ostream &operator<<( std::ostream &, const residue_name & );
	std::istream &operator>>( std::istream &, residue_name & );
	std::string   insert_string( const residue_name & );
	std::string   make_residue_name_string_with_insert_or_space( const residue_name & );
	residue_name  make_residue_name( const std::string & );

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_NAME_HPP
