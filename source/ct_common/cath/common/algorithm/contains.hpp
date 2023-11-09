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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONTAINS_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONTAINS_HPP

#include <map>
#include <set>
#include <unordered_map>

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/find_if.hpp>

/// \todo Consider adding predicate versions of all contains functions
namespace cath::common {

	/// \brief TODOCUMENT
	template <typename I, typename T>
	bool contains(const I &prm_begin, ///< TODOCUMENT
	              const I &prm_end,   ///< TODOCUMENT
	              const T &prm_value  ///< TODOCUMENT
	              ) {
		return ( prm_end != std::find( prm_begin, prm_end, prm_value ) );
	}

	/// \brief TODOCUMENT
	template <typename R, typename T>
	bool contains(const R &prm_range, ///< TODOCUMENT
	              const T &prm_value  ///< TODOCUMENT
	              ) {
		return ( ::std::cend( prm_range ) != boost::range::find( prm_range, prm_value ) );
	}



	/// \brief TODOCUMENT
	template <typename K, typename V, typename T>
	bool contains(const std::map<K, V> &prm_map,   ///< TODOCUMENT
	              const T              &prm_value  ///< TODOCUMENT
	              ) {
		return ( prm_map.count( prm_value ) > 0 );
	}

	/// \brief TODOCUMENT
	template <typename K, typename T>
	bool contains(const std::set<K> &prm_set,   ///< TODOCUMENT
	              const T           &prm_value  ///< TODOCUMENT
	              ) {
		return ( prm_set.count( prm_value ) > 0 );
	}

	/// \brief TODOCUMENT
	template <typename K, typename V, typename T>
	bool contains(const std::unordered_map<K, V> &prm_map,   ///< TODOCUMENT
	              const T                        &prm_value  ///< TODOCUMENT
	              ) {
		return ( prm_map.count( prm_value ) > 0 );
	}

	/// \brief TODOCUMENT
	template <typename K, typename V, typename T>
	bool contains(const boost::ptr_map<K, V> &prm_map,   ///< TODOCUMENT
	              const T                    &prm_value  ///< TODOCUMENT
	              ) {
		return ( prm_map.count( prm_value ) > 0 );
	}



	/// \brief TODOCUMENT
	template <class I, class P>
	bool contains_if(const I &prm_begin,   ///< TODOCUMENT
	                 const I &prm_end,     ///< TODOCUMENT
	                 P        prm_uni_pred ///< TODOCUMENT
	                 ) {
		return ( prm_end != std::find_if( prm_begin, prm_end, prm_uni_pred ) );
	}

	/// \brief TODOCUMENT
	template <class R, class P>
	bool contains_if(const R &prm_range,   ///< TODOCUMENT
	                 P        prm_uni_pred ///< TODOCUMENT
	                 ) {
		return ( ::std::cend( prm_range ) != boost::range::find_if( prm_range, prm_uni_pred ) );
	}



	/// \brief TODOCUMENT
	template <typename I>
	bool contains_adjacent_match(const I &prm_begin, ///< TODOCUMENT
	                             const I &prm_end    ///< TODOCUMENT
	                             ) {
		return ( prm_end != std::adjacent_find( prm_begin, prm_end ) );
	}

	/// \brief TODOCUMENT
	template <typename I, typename P>
	bool contains_adjacent_match(const I &prm_begin,   ///< TODOCUMENT
	                             const I &prm_end,     ///< TODOCUMENT
	                             P        prm_bin_pred ///< TODOCUMENT
	                             ) {
		return ( prm_end != std::adjacent_find( prm_begin, prm_end, prm_bin_pred ) );
	}

	/// \brief TODOCUMENT
	template <typename R>
	bool contains_adjacent_match(const R &prm_range ///< TODOCUMENT
	                             ) {
		return ( ::std::cend( prm_range ) != boost::range::adjacent_find( prm_range ) );
	}

	/// \brief TODOCUMENT
	template <typename R, typename P>
	bool contains_adjacent_match(const R &prm_range,   ///< TODOCUMENT
	                             P        prm_bin_pred ///< TODOCUMENT
	                             ) {
		return ( ::std::cend( prm_range ) != boost::range::adjacent_find( prm_range, prm_bin_pred ) );
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONTAINS_HPP
