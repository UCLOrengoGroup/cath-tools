/// \file
/// \brief The test tools header

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
///
/// \todo Extend all of the below to also provide checking of less-than functionality and consistency
///       (using common implementations as far as possible)

#ifndef TEST_TOOLS_H_INCLUDED
#define TEST_TOOLS_H_INCLUDED

#include <boost/test/unit_test.hpp>

#include "common/cpp14/cbegin_cend.h"

#include <iostream> // ***** TEMPORARY *****

namespace cath {
	namespace common {
		namespace test {

			namespace detail {

				/// \brief Helper to check that equality/inequality operators return opposite values that are both correct
				///
				/// This requires that T and U are equality/inequality comparable (with T on the LHS and U on the RHS).
				///
				/// Note. This does not check the symmetric operations (U == T, U != T) so it can be used if there is some
				/// situation where the reverse ==, != operations are forbidden. This is probably rare and the code below
				/// that uses this function always uses it in both directions
				///
				/// \todo This currently also requires that T and U provide ostream insertion operator overloads so that Boost Test
				///       can output comprehensive errors on failure. That's normally quite useful but it may be worth
				///       providing an option to only perform the comparisons with bool values for use on types for
				///       which it isn't worth providing an insertion operator overload. If motivated, this could probably be done
				///       with a compile-time check for whether T and U both have ostream insertion operator overloads.
				template <typename T, typename U>
				void check_equality_and_inequality_are_consistent(const T    &arg_value_1,        ///< The first argument to compare
				                                                  const U    &arg_value_2,        ///< The second argument to compare
				                                                  const bool &arg_should_be_equal ///< Whether the two arguments should compare equal
				                                                  ) {
					const bool are_equal     = ( arg_value_1 == arg_value_2 );
					const bool are_not_equal = ( arg_value_1 != arg_value_2 );
					BOOST_REQUIRE_NE( are_equal, are_not_equal );
					if ( arg_should_be_equal ) {
						BOOST_CHECK_EQUAL( arg_value_1, arg_value_2 );
					}
					else {
						BOOST_CHECK_NE( arg_value_1, arg_value_2 );
					}
					BOOST_CHECK_EQUAL( are_equal, arg_should_be_equal );
				}

				/// \brief Boost Test that the equality/inequality operators are correct for comparing arg_value against itself
				///        and, if copy-constructible, against a copy-constructed copy of itself (in both directions)
				///
				/// Requirements on T:
				///  * T is equality/inequality comparable with itself
				///  * (Currently) T has an ostream insertion operator overload
				template <bool T_is_copy_constructible>
				struct check_equality_operators_on_value_impl final {
					template <typename T>
					static void check(const T &arg_value ///< The value to be compared with itself and copy-constructed copies of itself
					                  ) {
						detail::check_equality_and_inequality_are_consistent(    arg_value,      arg_value,   true );
					}
				};

				/// \brief Boost Test that the equality/inequality operators are correct for comparing arg_value against itself
				///        and, if copy-constructible, against a copy-constructed copy of itself (in both directions)
				///
				/// Requirements on T:
				///  * T is equality/inequality comparable with itself
				///  * (Currently) T has an ostream insertion operator overload
				template <>
				struct check_equality_operators_on_value_impl<true> final {
					template <typename T>
					static void check(const T &arg_value ///< The value to be compared with itself and copy-constructed copies of itself
					                  ) {
						detail::check_equality_and_inequality_are_consistent( T( arg_value ),    arg_value,   true );
						detail::check_equality_and_inequality_are_consistent(    arg_value,   T( arg_value ), true );
						detail::check_equality_and_inequality_are_consistent(    arg_value,      arg_value,   true );
					}
				};
			}

			/// \brief Boost Test that the equality/inequality operators are correct for comparing arg_value against itself
			///        and, if copy-constructible, against a copy-constructed copy of itself (in both directions)
			///
			/// Requirements on T:
			///  * T is equality/inequality comparable with itself
			///  * (Currently) T has an ostream insertion operator overload
			template <typename T>
			void check_equality_operators_on_value(const T &arg_value ///< The value to be compared with itself and copy-constructed copies of itself
			                                       ) {
				detail::check_equality_operators_on_value_impl<std::is_copy_constructible<T>::value>::check( arg_value );
			}

			/// \brief Boost Test that the equality/inequality operators give correct values for two different objects
			///
			/// Requirements on T, U:
			///  * T, U are equality/inequality comparable with themselves and each other (in both directions)
			///  * (Currently) T, U are both copy constructible
			///  * (Currently) T, U both have ostream insertion operator overloads
			///
			/// Tests:
			///  - The equality/inequality operators both give correct results for arg1 versus arg2
			///  - Equality/inequality are also both correct if arguments are swapped (ie arg2 == arg2, arg2 != arg1)
			///  - The equality/inequality operators are both correct for comparing each argument to itself
			///    or to a copy-constructed copy of itself (in either direction)
			template <bool EQUAL, typename T, typename U>
			void check_equality_operators_on_vals(const T &arg1, ///< The first value to compare
			                                      const U &arg2  ///< The second value to compare
												  ) {
				// arg1 == arg1, arg2 == arg2
				check_equality_operators_on_value( arg1 );
				check_equality_operators_on_value( arg2 );

				// arg1 == arg2, arg2 == arg1
				detail::check_equality_and_inequality_are_consistent( arg1, arg2, EQUAL );
				detail::check_equality_and_inequality_are_consistent( arg2, arg1, EQUAL );
			}

			/// \brief Boost Test that the equality/inequality operators give correct values for two unequal objects
			///
			/// This is just a convenience interface that wraps check_equality_operators_on_vals
			template <typename T, typename U>
			void check_equality_operators_on_diff_vals(const T &arg1, ///< The first  value to compare
			                                           const U &arg2  ///< The second value to compare
			                                           ) {
				check_equality_operators_on_vals<false>( arg1, arg2 );
			}

			/// \brief Boost Test that the equality/inequality operators give correct values for two equal objects
			///
			/// This is just a convenience interface that wraps check_equality_operators_on_vals
			template <typename T, typename U>
			void check_equality_operators_on_equal_vals(const T &arg1, ///< The first value to compare
			                                            const U &arg2  ///< The second value to compare
			                                            ) {
				check_equality_operators_on_vals<true>( arg1, arg2 );
			}

			/// \brief
			template <typename R>
			void check_equality_operators_on_diff_vals_range(const R &arg_range ///< The range of values to compare
			                                                 ) {
				for (auto itr_a = common::cbegin( arg_range ); itr_a != common::cend( arg_range ); ++itr_a) {
					for (auto itr_b = common::cbegin( arg_range ); itr_b != common::cend( arg_range ); ++itr_b) {
						if ( itr_a != itr_b ) {
							check_equality_operators_on_diff_vals( *itr_a, *itr_b );
						}
					}
				}
			}
		}
	}
}

#endif
