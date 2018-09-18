/// \file
/// \brief The polymorphic_equality_comparable class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_POLYMORPHIC_COMPARISON_POLYMORPHIC_EQUALITY_COMPARABLE_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_POLYMORPHIC_COMPARISON_POLYMORPHIC_EQUALITY_COMPARABLE_HPP

#include <boost/concept/assert.hpp>
#include <boost/serialization/nvp.hpp>

#include "common/polymorphic_comparison/detail/dynamic_type_spaceship.hpp"
#include "common/polymorphic_comparison/detail/is_equal_with_same_dynamic_type_comparable.hpp"

#include <cassert>

namespace cath {
	namespace common {

		/// \brief Do Boost.Operators-like supplying of sensible operator==() for comparing objects within a class hierarchy
		///
		/// This is very similar to polymorphic_less_than_comparable and so the rest of the documentation is copied from that...
		///
		/// /// \copydetails polymorphic_less_than_comparable
		template <typename T>
		class polymorphic_equality_comparable {
		private:
			friend class boost::serialization::access;

			/// \brief Friend, non-member operator<() that uses dynamic type first, and then T's equal_with_same_dynamic_type()
			friend bool operator==(const T &prm_object1, ///< The first T to compare
			                       const T &prm_object2  ///< The first T to compare
			                       ) {
				// Make the compiler output sensible errors on any attempt to instantiate polymorphic_equality_comparable<>
				// on a type that doesn't provide a const `equal_with_same_dynamic_type(const T &)` method returning a bool-convertible type.
				BOOST_CONCEPT_ASSERT(( detail::is_equal_with_same_dynamic_type_comparable<T> ));

				// Perform a spaceship-style less-than comparison on the dynamic types of the two objects
				// which returns a tribool (true, false or indeterminate)
				const boost::logic::tribool dyn_type_cmp = detail::dynamic_type_spaceship::compare_lt( prm_object1, prm_object2 );

				// * If the result is either true (first object's dynamic type compares strictly less-than) or false (strictly greater-than), then
				//    * return false;
				// * otherwise:
				//    * return the result of calling equal_with_same_dynamic_type().
				return ! boost::logic::indeterminate( dyn_type_cmp  )
					? false
					: prm_object1.equal_with_same_dynamic_type( prm_object2 );
			}

			template<class archive> void serialize(archive &/*ar*/,
			                                       const unsigned int /*version*/) {
			}

		protected:
			polymorphic_equality_comparable() noexcept = default;
			~polymorphic_equality_comparable() noexcept = default;
			polymorphic_equality_comparable(const polymorphic_equality_comparable &) noexcept = default;
			polymorphic_equality_comparable(polymorphic_equality_comparable &&) noexcept = default;
			polymorphic_equality_comparable & operator=(const polymorphic_equality_comparable &) noexcept = default;
			polymorphic_equality_comparable & operator=(polymorphic_equality_comparable &&) noexcept = default;
		};
	
	} // namespace common
} // namespace cath

#endif
