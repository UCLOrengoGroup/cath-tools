/// \file
/// \brief The less_than_helper class header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Binaries project and then tweaked, eg namespaced in cath)
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

#ifndef LESS_THAN_HELPER_H_INCLUDED
#define LESS_THAN_HELPER_H_INCLUDED

#include <boost/logic/tribool.hpp>

#include <functional>

namespace cath {
	namespace common {

		/// \brief Facilitate implementations of less-than and equivalence functions for arbitrary types
		///
		/// It's quite fiddly to implement operator<() (and to a lesser extent operator==() ),
		/// because it requires stepping through the class's fields and for each:
		///  * returning true  if the first  value has a lower value
		///  * returning false if the second value has a lower value
		///  * continuing otherwise
		///
		/// This class helps to remove some of the boiler-plate from that task so it can be reduced to:
		///
		/// ~~~~~.cpp
		/// auto the_helper = make_less_than_helper( my_object_a, my_object_b );
		/// the_helper.register_comparison_field( &my_class::getter_1                            ); // Can register pointer-to-member-functions...
		/// the_helper.register_comparison_field( [] (const my_class &x) { return x.getter_2() } ); // ...or lambdas
		/// //...
		/// the_helper.register_comparison_field( [] (const my_class &x) { return x.getter_n() } );
		/// return final_less_than_result( the_helper );
		/// ~~~~~
		///
		/// The same can be done for equivalence by replacing the final_less_than_result() with
		/// final_equivalent_result() at the end.
		///
		/// This might be slightly slower than the full, hand-typed version because it doesn't
		/// immediately return when it finds a field with a clear true or false result.
		/// That said, it should be pretty close because the remaining calls to register_comparison_field()
		/// will return after checking `result`, without having to call the getters (ie comparison fields)
		/// or perform their comparisons.
		template <typename T>
		class less_than_helper final {
			/// \brief A reference to the first object
			///
			/// reference_wrapper is used over a reference to avoid needlessly making the class noncopyable
			std::reference_wrapper<const T> val_a;

			/// \brief A reference to the second object
			///
			/// reference_wrapper is used over a reference to avoid needlessly making the class noncopyable
			std::reference_wrapper<const T> val_b;

			/// \brief The result based on the fields seen *so far*
			///
			/// This is a spaceship less-than value (like Perl's spaceship <=> ). Vales:
			///  * true          : the first  object is less than the second (and no more testing is necessary)
			///  * false         : the second object is less than the first  (and no more testing is necessary)
			///  * indeterminate : it is not yet known whether the first object evaluates less than the second object
			boost::logic::tribool result = boost::logic::indeterminate;

		public:
			less_than_helper(const T &,
			                 const T &);

			template <typename U, typename V>
			void register_comparison_values(const U &,
			                                const V &);

			template <typename F>
			void register_comparison_field(F);

			boost::logic::tribool final_less_than_spaceship_result() const;
		};

		/// \brief TODOCUMENT
		template <typename T>
		less_than_helper<T>::less_than_helper(const T &arg_val_a, ///< TODOCUMENT
		                                      const T &arg_val_b  ///< TODOCUMENT
		                                      ) : val_a( arg_val_a ),
		                                          val_b( arg_val_b ) {
		}

		/// \brief TODOCUMENT
		template <typename T>
		template <typename U, typename V>
		void less_than_helper<T>::register_comparison_values(const U &arg_value_a, ///< TODOCUMENT
		                                                     const V &arg_value_b  ///< TODOCUMENT
		                                                     ) {
			if ( boost::logic::indeterminate( result ) ) {
				if ( arg_value_a < arg_value_b ) {
					result = true;
				}
				else if ( arg_value_b < arg_value_a ) {
					result = false;
				}
			}
		}

		/// \brief TODOCUMENT
		template <typename T>
		template <typename F>
		void less_than_helper<T>::register_comparison_field(F arg_getter ///< TODOCUMENT
		                                                    ) {
			if ( boost::logic::indeterminate( result ) ) {
				register_comparison_values(
					std::bind( arg_getter, std::cref( val_a.get() ) )(),
					std::bind( arg_getter, std::cref( val_b.get() ) )()
				);
			}
		}

		/// \brief TODOCUMENT
		template <typename T>
		boost::logic::tribool less_than_helper<T>::final_less_than_spaceship_result() const {
			return result;
		}

		/// \brief TODOCUMENT
		template <typename T>
		bool final_less_than_result(const less_than_helper<T> &arg_less_than_helper ///< TODOCUMENT
		                            ) {
			return static_cast<bool>( arg_less_than_helper.final_less_than_spaceship_result() );
		}

		/// \brief TODOCUMENT
		template <typename T>
		bool final_equivalent_result(const less_than_helper<T> &arg_less_than_helper ///< TODOCUMENT
		                             ) {
			return boost::logic::indeterminate( arg_less_than_helper.final_less_than_spaceship_result() );
		}

		/// \brief TODOCUMENT
		template <typename T>
		less_than_helper<T> make_less_than_helper(const T &arg_val_a, ///< TODOCUMENT
		                                          const T &arg_val_b  ///< TODOCUMENT
		                                          ) {
			return less_than_helper<T>( arg_val_a, arg_val_b );
		}

	}
}

#endif
