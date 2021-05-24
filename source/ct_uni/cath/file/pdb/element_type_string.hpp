/// \file
/// \brief The element_type_string class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_ELEMENT_TYPE_STRING_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_ELEMENT_TYPE_STRING_HPP

#include <boost/utility/string_view.hpp>

#include "cath/common/char_arr_type_aliases.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/string/string_parse_tools.hpp"
#include "cath/common/string/string_view_of_char_arr.hpp"
#include "cath/file/pdb/coarse_element_type.hpp"

#include <string>

namespace cath::file {

	/// \brief Hold an element type string along with a string_view to the non-whitespace
	///        part of the string.
	///
	/// The point of this class is to bundle the string_view and the string to which it refers together
	/// to ensure that the copy/move construction/assignment is done correctly.
	/// In all four cases, that requires the string_view to be reset to point to the new string.
	///
	/// \todo If there are any other uses for bundling a string_view with its string, then
	///       consider generalising this class.
	///
	/// \todo Come C++17 replace ::std::string_view with std::string_view
	class element_type_string final {
	private:
		/// \brief The untrimmed string describing the element type of this atom.
		///
		/// This is whitespace-untrimmed so that the correct whitespace can be returned whilst writing
		char_4_arr element_type_untrimmed;

		/// \brief The trimmed string describing the element type of this atom.
		///
		/// \todo Come C++17, replace ::std::string_view with std::string_view
		std::pair<char, char> trim_offsets;

	public:
		explicit constexpr element_type_string(char_4_arr);

		[[nodiscard]] constexpr const char_4_arr & get_element_type_untrimmed() const;
		[[nodiscard]] constexpr ::std::string_view get_element_type() const;
	};

	/// \brief Constructor from lvalue string
	constexpr element_type_string::element_type_string(char_4_arr prm_string ///< The source string
	                                                   ) : element_type_untrimmed( std::move( prm_string ) ), //< Don't change these brackets to braces - it breaks the build on the older Clang on Travis-CI
	                                                       trim_offsets {
	                                                       	common::dumb_trim_string_view_to_offsets<char>(
	                                                       		common::string_view_of_char_arr( element_type_untrimmed )
	                                                       	)
	                                                       } {
		if ( trim_offsets.first > trim_offsets.second ) {
			BOOST_THROW_EXCEPTION(common::out_of_range_exception("Error trimming of element_type_string"));
		}
	}

	/// \brief Getter for the original trimmed string
	constexpr const char_4_arr & element_type_string::get_element_type_untrimmed() const {
		return element_type_untrimmed;
	}

	/// \brief Getter for the trimmed element type string (as a string_view)
	constexpr ::std::string_view element_type_string::get_element_type() const {
		return {
			std::next( ::std::cbegin( element_type_untrimmed ), trim_offsets.first ),
			static_cast<size_t>( trim_offsets.second ) - static_cast<size_t>( trim_offsets.first )
		};
	}

	/// \brief Get the untrimmed string_view for the specified element_type_string
	constexpr ::std::string_view get_element_type_untrimmed_str_ref(const element_type_string &prm_element_type_string ///< The element_type_string to query
	                                                               ) {
		return common::string_view_of_char_arr( prm_element_type_string.get_element_type_untrimmed() );
	}

	/// \brief Get the coarse_element_type corresponding to the specified trimmed element string
	constexpr coarse_element_type get_coarse_element_type(const ::std::string_view &prm_trimmed_element_str ///< The trimmed element string (as it appears in PDB ATOM records)
	                                                      ) {
		if ( prm_trimmed_element_str.empty() || prm_trimmed_element_str.length() > 2 ) {
			return coarse_element_type::NON_CORE;
		}
		if ( prm_trimmed_element_str.front() == 'C' ) {
			if ( prm_trimmed_element_str == "CA" ) {
				return coarse_element_type::CARBON_ALPHA;
			}
			if ( prm_trimmed_element_str == "C"  ) {
				return coarse_element_type::CARBON;
			}
			if ( prm_trimmed_element_str == "CB" ) {
				return coarse_element_type::CARBON_BETA;
			}
			return coarse_element_type::NON_CORE;
		}
		if ( prm_trimmed_element_str == "N" ) {
			return coarse_element_type::NITROGEN;
		}
		if ( prm_trimmed_element_str == "O" ) {
			return coarse_element_type::OXYGEN;
		}
		return coarse_element_type::NON_CORE;
	}

	/// \brief Get the coarse_element_type corresponding to the specified element_type_string
	constexpr coarse_element_type get_coarse_element_type(const element_type_string &prm_element_type_string ///< The element_type_string to query
	                                                      ) {
		return get_coarse_element_type( prm_element_type_string.get_element_type() );
	}

} // namespace cath::file

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_ELEMENT_TYPE_STRING_HPP
