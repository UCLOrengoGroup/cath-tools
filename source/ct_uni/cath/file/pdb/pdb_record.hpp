/// \file
/// \brief The pdb_record header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_RECORD_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_RECORD_HPP

#include <boost/range/sub_range.hpp>

#include "cath/common/char_arr_type_aliases.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

#include <array>
#include <iosfwd>

namespace cath {
	namespace file {

		/// \brief Represent the type of record of a PDB file line
		enum class pdb_record : bool {
			ATOM,  ///< An ATOM record
			HETATM ///< A HETATM record
		};

		/// \brief Return the pdb_record corresponding to the specified substring
		///
		/// \pre prm_substring must refer to a string "ATOM  ", "ATOM" or "HETATM" else
		///      an invalid_argument_exception will be thrown
		///
		/// \todo Come C++17, convert this to use string_view rather than sub_range<const string>
		inline pdb_record pdb_rec_of_substring(const boost::sub_range<const std::string> &prm_substring ///< The substring to examine
		                                       ) {
			constexpr auto atom_6_str   = make_char_arr( "ATOM  " );
			constexpr auto atom_4_str   = make_char_arr( "ATOM"   );
			constexpr auto hetatm_6_str = make_char_arr( "HETATM" );
			if ( boost::range::equal( prm_substring, atom_6_str ) || boost::range::equal( prm_substring, atom_4_str ) ) {
				return pdb_record::ATOM;
			}
			if ( boost::range::equal( prm_substring, hetatm_6_str ) ) {
				return pdb_record::HETATM;
			}
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
				"Unable to recognise pdb_record type "
				+ std::string{
					common::cbegin( prm_substring ),
					common::cend  ( prm_substring )
				}
			));
		}

		/// \brief Return the pdb_record corresponding to the six chars after the specified start in the
		///        specified string
		///
		/// \pre The six chars must be "ATOM  " or "HETATM" else
		///      an invalid_argument_exception will be thrown
		inline pdb_record pdb_rec_of_six_chars_in_string(const std::string &prm_string,   ///< The string in which the substring of six chars appears
		                                                 const size_t      &prm_start = 0 ///< The index of the first char of the substring within the string
		                                                 ) {
			constexpr size_t NUM_CHARS = 6;
			if ( prm_start + NUM_CHARS > prm_string.length() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Cannot determine pdb_rec_of_six_chars_in_string() of a string that isn't long enough for six characters after the specified start"
				));
			}
			const auto begin_itr = std::next( prm_string.begin(), cath::debug_numeric_cast<ptrdiff_t>( prm_start ) );
			const auto end_itr   = std::next( begin_itr,          NUM_CHARS                                        );
			return pdb_rec_of_substring( boost::sub_range<const std::string>{ begin_itr, end_itr } );
		}

		/// \brief Return the pdb_record corresponding to the six chars in the specified string
		///
		/// \pre The string must contain exactly six chars and must be "ATOM  " or "HETATM" else
		///      an invalid_argument_exception will be thrown
		inline pdb_record pdb_rec_of_str(const std::string &prm_string ///< TODOCUMENT
		                                 ) {
			constexpr size_t NUM_CHARS = 6;
			if ( prm_string.length() != NUM_CHARS ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot convert a string to a pdb_record if it doesn't have 6 characters"));
			}
			return pdb_rec_of_six_chars_in_string( prm_string );
		}

		std::istream & operator>>(std::istream &,
		                          pdb_record &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_record &);

	} // namespace file
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_RECORD_HPP
