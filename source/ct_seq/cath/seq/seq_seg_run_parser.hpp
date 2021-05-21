/// \file
/// \brief The seq_seg_run class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_SEG_RUN_PARSER_HPP
#define _CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_SEG_RUN_PARSER_HPP

#include <boost/spirit/include/qi.hpp>

#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/seq/seq_seg_run.hpp"
#include "cath/seq/seq_type_aliases.hpp"

namespace cath {
	namespace seq {

		/// \brief Parse seq_seg_run from strings quickly and reusing the same intermediate vector
		///        to avoid reallocating every time
		class seq_seg_run_parser {
		private:
			/// \brief The vector in which to store the bounds
			residx_vec bounds;

		public:
			template <typename BegItr,
			          typename EndItr>
			seq_seg_run parse(const BegItr &,
			                  const EndItr &);
		};

		/// \brief Parse a seq_seg_run from the specified range of chars
		///
		/// \todo Change read_hit_list_from_istream() in calc_hit_list.cpp to use this
		template <typename BegItr,
		          typename EndItr>
		inline seq_seg_run seq_seg_run_parser::parse(const BegItr &prm_begin_itr, ///< The begin of the char range from which to parse the segments
		                                             const EndItr &prm_end_itr    ///< The end   of the char range from which to parse the segments
		                                             ) {
			// Clear any previous entries in bounds
			bounds.clear();

			// Prepare a lambda that pushes the specified value onto bounds
			const auto bounds_pusher = [&] (const residx_t &x) { bounds.push_back( x ); };

			// Parse the text
			auto parse_itr = prm_begin_itr;
			const bool ok = boost::spirit::qi::parse(
				parse_itr,
				prm_end_itr,
				   boost::spirit::uint_   [ bounds_pusher ]
				>> boost::spirit::qi::omit[ '-'           ]
				>> boost::spirit::uint_   [ bounds_pusher ]
				>> *(
					   boost::spirit::qi::char_( ",_" )
					>> boost::spirit::uint_   [ bounds_pusher ]
					>> boost::spirit::qi::omit[ "-"           ]
					>> boost::spirit::uint_   [ bounds_pusher ]
				)
			);

			// Check for problems
			if ( ! ok || parse_itr != prm_end_itr ) {
				BOOST_THROW_EXCEPTION(common::runtime_error_exception( "Error on attempt to parse line : " + std::string{ prm_begin_itr, prm_end_itr } ));
			}
			if ( bounds.empty() ) {
				BOOST_THROW_EXCEPTION(common::runtime_error_exception( "No bounds" ));
			}
			if ( bounds.size() % 2 != 0 ) {
				BOOST_THROW_EXCEPTION(common::runtime_error_exception( "Odd number of bounds" ));
			}

			// Return the result of building a seq_seg_run from bounds
			return seq_seg_run{ segments_from_bounds( bounds ) };
		}

	} // namespace seq
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_SEG_RUN_PARSER_HPP
