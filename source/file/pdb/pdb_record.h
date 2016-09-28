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

#ifndef PDB_RECORD_H_INCLUDED
#define PDB_RECORD_H_INCLUDED

#include <boost/range/sub_range.hpp>

#include "common/cpp14/cbegin_cend.h"
#include "common/debug_numeric_cast.h"
#include "exception/invalid_argument_exception.h"

#include <array>
#include <iosfwd>

namespace cath {
	namespace file {

		/// \brief Represent the type of record of a PDB file line
		enum class pdb_record {
			ATOM,  ///< An ATOM record
			HETATM ///< A HETATM record
		};

		/// \brief Return the pdb_record corresponding to the specified substring
		///
		/// \pre arg_substring must refer to a string "ATOM  ", "ATOM" or "HETATM" else
		///      an invalid_argument_exception will be thrown
		///
		/// \todo Come C++17, convert this to use string_view rather than sub_range<const string>
		inline pdb_record pdb_rec_of_substring(const boost::sub_range<const std::string> &arg_substring ///< The substring to examine
		                                       ) {
			constexpr std::array<char, 6> atom_6_str   = { { 'A','T','O','M',' ',' ', } };
			constexpr std::array<char, 4> atom_4_str   = { { 'A','T','O','M',         } };
			constexpr std::array<char, 6> hetatm_6_str = { { 'H','E','T','A','T','M', } };
			if ( boost::range::equal( arg_substring, atom_6_str ) || boost::range::equal( arg_substring, atom_4_str ) ) {
				return pdb_record::ATOM;
			}
			else if ( boost::range::equal( arg_substring, hetatm_6_str ) ) {
				return pdb_record::HETATM;
			}
			else {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Unable to recognise pdb_record type "
					+ std::string{
						common::cbegin( arg_substring ),
						common::cend  ( arg_substring )
					}
				));
				return pdb_record::ATOM; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
			}
		}

		/// \brief Return the pdb_record correspoding to the six chars after the specified start in the
		///        specified string
		///
		/// \pre The six chars must be "ATOM  " or "HETATM" else
		///      an invalid_argument_exception will be thrown
		inline pdb_record pdb_rec_of_six_chars_in_string(const std::string &arg_string,   ///< The string in which the substring of six chars appears
		                                                 const size_t      &arg_start = 0 ///< The index of the first char of the substring within the string
		                                                 ) {
			constexpr size_t NUM_CHARS = 6;
			if ( arg_start + NUM_CHARS > arg_string.length() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Cannot determine pdb_rec_of_six_chars_in_string() of a string that isn't long enough for six characters after the specified start"
				));
			}
			const auto begin_itr = std::next( arg_string.begin(), cath::debug_numeric_cast<ptrdiff_t>( arg_start ) );
			const auto end_itr   = std::next( begin_itr,          NUM_CHARS                                        );
			return pdb_rec_of_substring( boost::sub_range<const std::string>{ begin_itr, end_itr } );
		}

		/// \brief Return the pdb_record corresponding to the six chars in the specified string
		///
		/// \pre The string must contain exactly six chars and must be "ATOM  " or "HETATM" else
		///      an invalid_argument_exception will be thrown
		inline pdb_record pdb_rec_of_str(const std::string &arg_string ///< TODOCUMENT
		                                 ) {
			constexpr size_t NUM_CHARS = 6;
			if ( arg_string.length() != NUM_CHARS ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot convert a string to a pdb_record if it doesn't have 6 characters"));
			}
			return pdb_rec_of_six_chars_in_string( arg_string );
		}

		std::istream & operator>>(std::istream &,
		                          pdb_record &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_record &);

	}
}

#endif
