/// \file
/// \brief The prc_scores_line_parser class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PRC_SCORES_FILE_DETAIL_PRC_SCORES_LINE_PARSER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PRC_SCORES_FILE_DETAIL_PRC_SCORES_LINE_PARSER_HPP

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/file/prc_scores_file/prc_scores_entry.hpp"

#include <functional>

namespace cath {
	namespace file {

		/// \brief Helper for parsing PRC scores lines that, when reused over multiple lines,
		///        reuses a vector<string> and hence avoids wasting time on needless reallocations
		class prc_scores_line_parser final {
		private:
			/// \brief Reused vector<string> of for storing parts parsed from lines
			str_vec string_parts;

		public:
			prc_scores_entry parse_line(const std::string &);
		};

		/// \brief Parse the specified string into a prc_scores_entry
		///
		/// This has been optimised based on profiling because it can be used quite intensively
		/// if used on many comprehensive PRC searches. Optimisations:
		///  * Uses explicit char comparisons for spaces rather than locale-based boost::algorithm::is_space / std::isspace()
		///  * Searches along the string rather than using split()
		inline prc_scores_entry prc_scores_line_parser::parse_line(const std::string &prm_prc_line ///< The PRC scores file line to be parsed
		                                                           ) {
			const auto is_space_in_prc     = [] (const char &x) { return ( x == ' ' || x == '\t' || x == '\n' ); };
			const auto is_not_space_in_prc = [] (const char &x) { return ( x != ' ' && x != '\t' && x != '\n' ); };
			auto       start_itr = ::std::cbegin( prm_prc_line );
			const auto end_itr   = ::std::cend  ( prm_prc_line );
			auto       next_itr  = std::find_if( start_itr, end_itr, is_space_in_prc );
			size_t     ctr       = 0;
			while ( start_itr != end_itr ) {
				if ( ctr >= string_parts.size() ) {
					string_parts.emplace_back( start_itr, next_itr );
				}
				string_parts[ ctr ].assign( start_itr, next_itr );
				start_itr = std::find_if( next_itr,  end_itr, is_not_space_in_prc );
				next_itr  = std::find_if( start_itr, end_itr, is_space_in_prc     );
				++ctr;
			}

			if ( ctr != 12 ) {
				BOOST_THROW_EXCEPTION(common::runtime_error_exception("Unable to parse prc_scores_entry from line that doesn't contain 12 parts"));
			}

			return prc_scores_entry{
				       string_parts[  0 ],
				stoul( string_parts[  1 ] ),
				stoul( string_parts[  2 ] ),
				stoul( string_parts[  3 ] ),
				stoul( string_parts[  4 ] ),
				       string_parts[  5 ],
				stoul( string_parts[  6 ] ),
				stoul( string_parts[  7 ] ),
				stoul( string_parts[  8 ] ),
				stod ( string_parts[  9 ] ),
				stod ( string_parts[ 10 ] ),
				stod ( string_parts[ 11 ] )
				};
		}

	} // namespace file
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PRC_SCORES_FILE_DETAIL_PRC_SCORES_LINE_PARSER_HPP
