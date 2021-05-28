/// \file
/// \brief The is_less_than_with_same_dynamic_type_comparable concept template class

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_POLYMORPHIC_COMPARISON_DETAIL_IS_LESS_THAN_WITH_SAME_DYNAMIC_TYPE_COMPARABLE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_POLYMORPHIC_COMPARISON_DETAIL_IS_LESS_THAN_WITH_SAME_DYNAMIC_TYPE_COMPARABLE_HPP

#include <boost/concept_check.hpp>

namespace cath::common::detail {

	/// \brief Check type provides a const less_than_with_same_dynamic_type() method
	///        whose return type is convertible to bool.
	///
	/// This uses the Boost Concept tools.
	/// See the "Creating Concept Checking Classes" page for the Boost Concept Check Library pages
	/// (currently http://www.boost.org/libs/concept_check/creating_concepts.htm)
	template <typename X>
	struct is_less_than_with_same_dynamic_type_comparable {
	public:
		/// \brief The bit that does the actual concept-checking work
		BOOST_CONCEPT_USAGE( is_less_than_with_same_dynamic_type_comparable ) {
			[[maybe_unused]] bool answer( value_1.less_than_with_same_dynamic_type( value_2 ) );
		}

		/// \brief Ctor just to prevent compiler complaining that value references won't ever get initialised
		is_less_than_with_same_dynamic_type_comparable(const X &prm_value_1,
		                                               const X &prm_value_2
		                                               ) : value_1( prm_value_1 ),
		                                                   value_2( prm_value_2 ) {
		}

	private:
		/// \brief Value stored as const reference to permit X to be abstract
		///        (and to ensure less_than_with_same_dynamic_type() doesn't require non-const)
		const X &value_1;

		/// \brief Value stored as const reference to permit X to be abstract
		///        (and to ensure less_than_with_same_dynamic_type() doesn't require non-const)
		const X &value_2;
	};

} // namespace cath::common::detail

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_POLYMORPHIC_COMPARISON_DETAIL_IS_LESS_THAN_WITH_SAME_DYNAMIC_TYPE_COMPARABLE_HPP
