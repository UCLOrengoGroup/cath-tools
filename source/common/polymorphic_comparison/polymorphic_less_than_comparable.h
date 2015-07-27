/// \file
/// \brief The polymorphic_less_than_comparable class header

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

#ifndef POLYMORPHIC_LESS_THAN_COMPARABLE_H_INCLUDED
#define POLYMORPHIC_LESS_THAN_COMPARABLE_H_INCLUDED

#include <boost/concept/assert.hpp>
#include <boost/serialization/nvp.hpp>

#include "common/polymorphic_comparison/detail/dynamic_type_spaceship.h"
#include "common/polymorphic_comparison/detail/is_less_than_with_same_dynamic_type_comparable.h"

#include <cassert>

namespace cath {
	namespace common {

		/// \brief Do Boost.Operators-like supplying of sensible operator<() for comparing objects within a class hierarchy
		///
		/// In the style of Boost.Operators, this involves two roles:
		///  * requiring client classes to implement less_than_with_same_dynamic_type()
		///  * providing an operator<() that first sorts on dynamic type and then
		///    falls back on do_less_than_with_same_dynamic_type() for pairs with matching dynamic type
		///
		/// This template class should be used via the Curiously Recurring Template Pattern (CRTP),
		/// ie to use it, make your class T privately inherit from polymorphic_less_than_comparable<T>.
		/// This will then require that T implements less_than_with_same_dynamic_type(const T &).
		///
		/// If the class hierarchy is only required to sort based on dynamic type, then the base class
		/// only need implement its less_than_with_same_dynamic_type() to always return false.
		/// If further sorting is required, less_than_with_same_dynamic_type() should be made into an NVI
		/// pass-through to a pure-virtual do_less_than_with_same_dynamic_type(). Eg:
		///
		/// \code
		/// class my_class : private common::polymorphic_less_than_comparable< my_class > {
		/// private:
		///   virtual bool do_less_than_with_same_dynamic_type(const my_class &) const = 0;
		///  // ...
		///
		/// public:
		///   bool less_than_with_same_dynamic_type(const my_class &) const;
		///  // ...
		///
		/// }
		///
		/// bool my_class::less_than_with_same_dynamic_type(const my_class &arg_other_object
		///                                                 ) const {
		///   return do_less_than_with_same_dynamic_type(arg_other_object);
		/// }
		/// \endcode
		///
		/// Questions
		/// ---------
		///
		/// * Why is this implemented in this way? *
		///
		/// Because this is basically doing the same sort of thing as Boost.Operators and this
		/// is how Boost.Operators does it.
		///
		/// * In particular, why does this need to use the CRTP? *
		///
		/// So that the compiler rejects erroneous attempts to compare types from
		/// two separate inheritance hierarchies that both happen to use this class.
		///
		/// * In particular, what's the story with the friend function defined in the class? *
		///
		/// This is how Boost.Operators does it. The library contains the comment:
		/// > Note that friend functions defined in a class are implicitly inline.
		/// > See the C++ std, 11.4 [class.friend] paragraph 5
		///
		/// ---
		///
		/// General operator<() policy
		/// -------------------------
		///
		/// ???
		///
		/// (note that there is no claim that this policy is yet comprehensively applied at the time of writing)
		///
		/// * Use less_than_helper
		template <typename T>
		class polymorphic_less_than_comparable {
		private:
			friend class boost::serialization::access;

			/// \brief Friend, non-member operator<() that uses dynamic type first, and then T's less_than_with_same_dynamic_type()
			friend bool operator<(const T &arg_object1, ///< The first T to compare
			                      const T &arg_object2  ///< The first T to compare
			                      ) {
				// Make the compiler output sensible errors on any attempt to instantiate polymorphic_less_than_comparable<>
				// on a type that doesn't provide a const `less_than_with_same_dynamic_type(const T &)` method returning a bool-convertible type.
				BOOST_CONCEPT_ASSERT(( detail::is_less_than_with_same_dynamic_type_comparable<T> ));

				// Perform a spaceship-style less-than comparison on the dynamic types of the two objects
				// which returns a tribool (true, false or indeterminate)
				const boost::logic::tribool dyn_type_cmp = detail::dynamic_type_spaceship::compare_lt( arg_object1, arg_object2 );



				// If the result is true (first object's dynamic type compares strictly less-than) or false (strictly greater-than),
				// then just (the bool equivalent of) that value; otherwise, return the result of calling less_than_with_same_dynamic_type().
				return ( dyn_type_cmp || ! dyn_type_cmp ) ? static_cast<bool>(dyn_type_cmp)
				                                          : arg_object1.less_than_with_same_dynamic_type( arg_object2 );
			}

			template<class archive> void serialize(archive &/*ar*/,
			                                       const unsigned int /*version*/) {
			}

		protected:
			~polymorphic_less_than_comparable() noexcept = default;
		};
	
	}
}

#endif
