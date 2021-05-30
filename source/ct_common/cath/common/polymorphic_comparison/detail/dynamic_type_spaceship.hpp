/// \file
/// \brief The dynamic_type_spaceship class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_POLYMORPHIC_COMPARISON_DETAIL_DYNAMIC_TYPE_SPACESHIP_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_POLYMORPHIC_COMPARISON_DETAIL_DYNAMIC_TYPE_SPACESHIP_HPP

#include <boost/logic/tribool.hpp>

#include <cassert>
#include <typeinfo>

#include "cath/common/config.hpp"

namespace cath::common::detail {

	/// \brief Class containing public static function for performing a spaceship-style less-than
	///        comparison on the dynamic type of two objects of matching static type
	class dynamic_type_spaceship final {
	public:
		dynamic_type_spaceship() = delete;

		template <typename T>
		static boost::logic::tribool compare_lt(const T &,
		                                        const T &);
	};

	/// \brief Perform a spaceship-style less-than comparison on the dynamic type of two objects of matching static type
	///
	/// This returns:
	///  * true          if the first object's dynamic type compares strictly less    than the second object's,
	///  * false         if the first object's dynamic type compares strictly greater than the second object's or
	///  * indeterminate otherwise.
	template <typename T>
	boost::logic::tribool dynamic_type_spaceship::compare_lt(const T &prm_object1, ///< The first  object whose dynamic type should be compared
	                                                         const T &prm_object2  ///< The second object whose dynamic type should be compared
	                                                         ) {
		if ( typeid( prm_object1 ).before( typeid( prm_object2 ) ) ) {
			return true;
		}
		if ( typeid( prm_object2 ).before( typeid( prm_object1 ) ) ) {
			return false;
		}
		if constexpr ( IS_IN_DEBUG_MODE ) {
			// Neither dynamic type compares before the other but nor do they compare equal
			assert( typeid( prm_object1 ) == typeid( prm_object2 ) );
		}
		return boost::logic::indeterminate;
	}

} // namespace cath::common::detail

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_POLYMORPHIC_COMPARISON_DETAIL_DYNAMIC_TYPE_SPACESHIP_HPP
