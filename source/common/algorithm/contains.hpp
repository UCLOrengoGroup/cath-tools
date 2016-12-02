/// \file
/// \brief The contains class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_CONTAINS_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_CONTAINS_H

#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include "common/cpp14/cbegin_cend.hpp"

#include <map>
#include <set>

/// \todo Consider adding predicate versions of all contains functions
namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		template <typename I, typename T>
		bool contains(const I &arg_begin, ///< TODOCUMENT
		              const I &arg_end,   ///< TODOCUMENT
		              const T &arg_value  ///< TODOCUMENT
		              ) {
			return ( arg_end != std::find( arg_begin, arg_end, arg_value ) );
		}

		/// \brief TODOCUMENT
		template <typename R, typename T>
		bool contains(const R &arg_range, ///< TODOCUMENT
		              const T &arg_value  ///< TODOCUMENT
		              ) {
			return ( common::cend( arg_range ) != boost::range::find( arg_range, arg_value ) );
		}



		/// \brief TODOCUMENT
		template <typename K, typename V, typename T>
		bool contains(const std::map<K, V> &arg_map,   ///< TODOCUMENT
		              const T              &arg_value  ///< TODOCUMENT
		              ) {
			return ( arg_map.count( arg_value ) > 0 );
		}

		/// \brief TODOCUMENT
		template <typename K, typename T>
		bool contains(const std::set<K> &arg_set,   ///< TODOCUMENT
		              const T           &arg_value  ///< TODOCUMENT
		              ) {
			return ( arg_set.count( arg_value ) > 0 );
		}

		/// \brief TODOCUMENT
		template <typename K, typename V, typename T>
		bool contains(const boost::ptr_map<K, V> &arg_map,   ///< TODOCUMENT
		              const T                    &arg_value  ///< TODOCUMENT
		              ) {
			return ( arg_map.count( arg_value ) > 0 );
		}



		/// \brief TODOCUMENT
		template <class I, class P>
		bool contains_if(const I &arg_begin,   ///< TODOCUMENT
		                 const I &arg_end,     ///< TODOCUMENT
		                 P        arg_uni_pred ///< TODOCUMENT
		                 ) {
			return ( arg_end != std::find_if( arg_begin, arg_end, arg_uni_pred ) );
		}

		/// \brief TODOCUMENT
		template <class R, class P>
		bool contains_if(const R &arg_range,   ///< TODOCUMENT
		                 P        arg_uni_pred ///< TODOCUMENT
		                 ) {
			return ( common::cend( arg_range ) != boost::range::find_if( arg_range, arg_uni_pred ) );
		}



		/// \brief TODOCUMENT
		template <typename I>
		bool contains_adjacent_match(const I &arg_begin, ///< TODOCUMENT
		                             const I &arg_end    ///< TODOCUMENT
		                             ) {
			return ( arg_end != std::adjacent_find( arg_begin, arg_end ) );
		}

		/// \brief TODOCUMENT
		template <typename I, typename P>
		bool contains_adjacent_match(const I &arg_begin,   ///< TODOCUMENT
		                             const I &arg_end,     ///< TODOCUMENT
		                             P        arg_bin_pred ///< TODOCUMENT
		                             ) {
			return ( arg_end != std::adjacent_find( arg_begin, arg_end, arg_bin_pred ) );
		}

		/// \brief TODOCUMENT
		template <typename R>
		bool contains_adjacent_match(const R &arg_range ///< TODOCUMENT
		                             ) {
			return ( common::cend( arg_range ) != boost::range::adjacent_find( arg_range ) );
		}

		/// \brief TODOCUMENT
		template <typename R, typename P>
		bool contains_adjacent_match(const R &arg_range,   ///< TODOCUMENT
		                             P        arg_bin_pred ///< TODOCUMENT
		                             ) {
			return ( common::cend( arg_range ) != boost::range::adjacent_find( arg_range, arg_bin_pred ) );
		}
	} // namespace common
} // namespace cath

#endif
