/// \file
/// \brief The split_build header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_STRING_ALGORITHM_SPLIT_BUILD_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_STRING_ALGORITHM_SPLIT_BUILD_HPP

#include <boost/algorithm/string/split.hpp>

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		template <typename OutputContainer,
		          typename RangeT,
		          typename PredicateT>
		inline OutputContainer split_build(const RangeT                               &prm_input,                                          ///< TODOCUMENT
		                                   PredicateT                                  prm_pred,                                           ///< TODOCUMENT
		                                   boost::algorithm::token_compress_mode_type  prm_compress = boost::algorithm::token_compress_off ///< TODOCUMENT
		                                   ) {
			OutputContainer output_container;
			boost::algorithm::split(
				output_container,
				prm_input,
				prm_pred,
				prm_compress
			);
			return output_container;
		}
	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_STRING_ALGORITHM_SPLIT_BUILD_HPP
