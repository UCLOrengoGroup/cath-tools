/// \file
/// \brief The vector_of_vector header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_CONTAINER_VECTOR_OF_VECTOR_H
#define _CATH_TOOLS_SOURCE_COMMON_CONTAINER_VECTOR_OF_VECTOR_H

#include <boost/range/combine.hpp>
#include <boost/range/irange.hpp>

#include "common/size_t_literal.hpp"
#include "exception/invalid_argument_exception.hpp"

#include <iostream>
#include <vector>

namespace cath {
	namespace common {

		/// \brief A rectangular vector of vector of T that stores in a single vector
		///
		/// This is a first step towards a windowed-matrix class and is being added
		/// as a first step towards making the main SSAP matrix into local rather than
		/// global variables
		template <typename T>
		class vector_of_vector final {
		private:
			/// \brief The length along the first dimension
			size_t length_a = 0;

			/// \brief The length along the second dimension
			size_t length_b = 0;

			/// \brief The single vector that stores all the elements
			std::vector<T> data;

			/// \brief Type alias for std::vector<T>
			using vec_t = std::vector<T>;

		public:
			/// \brief Type alias for the vector<T>'s size_type
			using size_type = typename vec_t::size_type;

			/// \brief Type alias for the vector<T>'s reference
			using reference = typename vec_t::reference;

			/// \brief Type alias for the vector<T>'s const_reference
			using const_reference = typename vec_t::const_reference;

			vector_of_vector() = default;
			explicit vector_of_vector(const std::initializer_list<std::initializer_list<T>> &);
			vector_of_vector(const size_t &,
			                 const size_t &,
			                 const T & = T{});

			const size_t & get_length_a() const;
			const size_t & get_length_b() const;

			vector_of_vector & assign(const size_t &,
			                          const size_t &,
			                          const T &);

			vector_of_vector & resize(const size_t &,
			                          const size_t &,
			                          const T & = T{});

			reference get(const size_t &,
			              const size_t &);
			const_reference get(const size_t &,
			                    const size_t &) const;
			vector_of_vector & set(const size_t &,
			                       const size_t &,
			                       const T &);
		};

		/// \brief Ctor from a nested initializer_list
		///
		/// This may be quite inefficient and is currently just here for convenience
		/// re test data
		///
		/// \pre arg_init_list must be rectangular (ie all inner intializer_list<T>s should have the same size)
		///      else an invalid_argument_exception will be thrown
		template <typename T>
		inline vector_of_vector<T>::vector_of_vector(const std::initializer_list<std::initializer_list<T>> &arg_init_list ///< The nested initializer_list containing the data with which this vector_of_vector should be initialised
		                                             ) : length_a( arg_init_list.size()                                 ),
		                                                 length_b( ( length_a > 0 ) ? arg_init_list.begin()->size() : 0 ),
		                                                 data    ( length_a * length_b                                  ) {
			for (const auto &x_idx_and_rng : boost::combine( boost::irange( 0_z, arg_init_list.size() ), arg_init_list ) ) {
				const size_t                   &x_idx = boost::get<0>( x_idx_and_rng );
				const std::initializer_list<T> &rng   = boost::get<1>( x_idx_and_rng );

				if ( length_b != rng.size() ) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception(
						"The nested initializer_list for vector_of_vector ctor must be rectangular (ie all inner intializer_list<T>s should have the same size)"
					));
				}
				for (const auto &y_idx_and_val : boost::combine( boost::irange( 0_z, rng.size() ), rng ) ) {
					const size_t &y_idx = boost::get<0>( y_idx_and_val );
					const T      &value = boost::get<1>( y_idx_and_val );
					this->set( x_idx, y_idx, value );
				}
			}
		}

		/// \brief Construct from the two specified lengths and the (optional) value with which to populate
		template <typename T>
		inline vector_of_vector<T>::vector_of_vector(const size_t &arg_length_a, ///< The length of the first  dimension
		                                             const size_t &arg_length_b, ///< The length of the second dimension
		                                             const T      &arg_value     ///< The (optional) value with which to populate
		                                             ) : length_a ( arg_length_a                           ),
		                                                 length_b ( arg_length_b                           ),
		                                                 data     ( arg_length_a * arg_length_b, arg_value ) {
		}

		/// \brief Getter for the length of the first dimension
		template <typename T>
		const size_t & vector_of_vector<T>::get_length_a() const {
			return length_a;
		}

		/// \brief Getter for the length of the second dimension
		template <typename T>
		const size_t & vector_of_vector<T>::get_length_b() const {
			return length_b;
		}

		/// \brief Assign from the two specified lengths and the value with which to populate
		///
		/// This is analogous to std::vector::assign
		template <typename T>
		vector_of_vector<T> & vector_of_vector<T>::assign(const size_t &arg_length_a, ///< The new length of the first  dimension
		                                                  const size_t &arg_length_b, ///< The new length of the second dimension
		                                                  const T      &arg_value     ///< The value with which to populate
		                                                  ) {
			data.assign( arg_length_a * arg_length_b, arg_value );
			length_a = arg_length_a;
			length_b = arg_length_b;
			return *this;
		}

		/// \brief Resize (if necessary) to the specified lengths, using the specified value
		///        if any extra values are to be added
		///
		/// NOTE: this is different to std::vector::resize; unlike that method, this will
		///       overwrite with arg_value if either of the specified lengths differs from
		///       the current values. If neither length is to be changed, this will have no
		///       effect.
		template <typename T>
		vector_of_vector<T> & vector_of_vector<T>::resize(const size_t &arg_length_a, ///< The new length of the first  dimension
		                                                  const size_t &arg_length_b, ///< The new length of the second dimension
		                                                  const T      &arg_value     ///< The value with which to populate any new elements that are to be created
		                                                  ) {
			if ( length_a != arg_length_a || length_b != arg_length_b ) {
				assign( arg_length_a, arg_length_b, arg_value );
			}
			return *this;
		}

		/// \brief Get the element at the specified indices (non-const overload)
		template <typename T>
		inline auto vector_of_vector<T>::get(const size_t &arg_index_a, ///< The index in the first  dimension
		                                     const size_t &arg_index_b  ///< The index in the second dimension
		                                     ) -> reference {
			return data[ ( arg_index_a * length_b ) + arg_index_b ];
		}

		/// \brief Get the element at the specified indices (const overload)
		template <typename T>
		inline auto vector_of_vector<T>::get(const size_t &arg_index_a, ///< The index in the first  dimension
		                                     const size_t &arg_index_b  ///< The index in the second dimension
		                                     ) const -> const_reference {
			return data[ ( arg_index_a * length_b ) + arg_index_b ];
		}

		/// \brief Set the element at the specified indices to the specified value
		template <typename T>
		inline vector_of_vector<T> & vector_of_vector<T>::set(const size_t &arg_index_a, ///< The index in the first  dimension
		                                                      const size_t &arg_index_b, ///< The index in the second dimension
		                                                      const T      &arg_value    ///< The value to set
		                                                      ) {
			data[ ( arg_index_a * length_b ) + arg_index_b ] = arg_value;
			return *this;
		}

	} // namespace common
} // namespace cath

#endif
