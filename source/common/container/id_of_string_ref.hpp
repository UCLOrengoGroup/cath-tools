/// \file
/// \brief The id_of_string_ref header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_CONTAINER_ID_OF_STRING_REF_H
#define _CATH_TOOLS_SOURCE_COMMON_CONTAINER_ID_OF_STRING_REF_H

#include "common/cpp14/cbegin_cend.hpp"

#include <string>
#include <unordered_map>

namespace cath {
	namespace common {

		namespace detail {

			/// \brief Wrap a reference_wrapper<T> for use in an unordered_map with hashing on the T value
			///
			/// Equality is also performed on the T value.
			///
			/// \todo If this continues to be used, move it to its own header
			template <typename T>
			class ref_wrap_uom_wrap final {
			private:
				/// \brief The reference_wrapper of the value
				std::reference_wrapper<T> value;

			public:
				/// \brief Ctor from a reference
				ref_wrap_uom_wrap(T &arg_value ///< The value from which to construct
				                  ) : value{ arg_value } {
				}

				/// \brief Prevent construction from an rvalue
				ref_wrap_uom_wrap(T&& x ) = delete;

				/// \brief Ctor from a reference_wrapper<T>
				ref_wrap_uom_wrap(const std::reference_wrapper<T> &arg_value ///< The reference_wrapper<T> from which this should be constructed
				                  ) : value{ arg_value } {
				}
				/// \brief Copy ctor
				ref_wrap_uom_wrap(const ref_wrap_uom_wrap<T> &arg_value ///< The ref_wrap_uom_wrap from which this should be constructed
				                  ) : value{ arg_value.get_ref_wrap() } {
				}

				/// \brief Assignment from a reference_wrapper<T>
				ref_wrap_uom_wrap & operator=(const std::reference_wrapper<T> &arg_rhs ///< The reference_wrapper<T> to assign
				                              ) {
					value = arg_rhs;
					return *this;
				}
				/// \brief Copy assignment operator
				ref_wrap_uom_wrap & operator=(const ref_wrap_uom_wrap<T> &arg_rhs ///< The ref_wrap_uom_wrap to assign
				                              ) {
					value = arg_rhs.get_ref_wrap();
					return *this;
				}

				/// \brief Get the reference_wrapper
				const std::reference_wrapper<T> & get_ref_wrap() {
					return value;
				}

				/// \brief Conversion to a reference to T
				operator T & () const {
					return value.get();
				}

				/// \brief Get a reference to T
				T & get() const {
					return value.get();
				}

				/// \brief Invoke the wrapped value on the specified arguments
				template  <typename... Ts>
				std::result_of_t< T & (Ts &&...) > operator()(Ts &&... args ///< The values to pass to the value
				                                              ) const {
					return value( std::forward<Ts>( args )... );
				}
			};

			/// \brief Equality operator for ref_wrap_uom_wrap
			template <typename T>
			inline bool operator==(const ref_wrap_uom_wrap<T> &arg_lhs, ///< The first  ref_wrap_uom_wrap to compare
			                       const ref_wrap_uom_wrap<T> &arg_rhs  ///< The second ref_wrap_uom_wrap to compare
			                       ) {
				return ( arg_lhs.get() == arg_rhs.get() );
			}

			/// \brief Hasher for a ref_wrap_uom_wrap<T>
			template <typename T>
			struct ref_wrap_hasher final {
				/// \brief The function operator that performs the hash on the T value
				///        using std::hash<decay_t<T>>
				size_t operator()(const ref_wrap_uom_wrap<T> &arg_value
				                  ) const {
					return std::hash<std::decay_t<T>>{}( arg_value.get() );
				}
			};
		}

		/// \brief A type alias for a reference_wrapper of const string
		///
		/// \todo If this continues to be used, move it to common/type_aliases.hpp
		using string_cref = std::reference_wrapper<const std::string>;

		/// \brief Map reference_wrapper<const string> to numeric IDs that count incrementally from 0
		///
		/// Note: this is an unordered_map of reference_wrapper<string> not of boost::string_ref / std::string_view
		///
		/// This makes lookups easier
		///
		/// This class isn't thread-safe
		///
		/// \TODO Consider adding (or changing this into) id_of_string_ref_ref / id_of_string_ref_view.
		///       The Boost string_ref currently has a std::hash specialisation
		///       but it's #if-ed out (#if 0).
		///
		/// \TODO Alternatively consider using Boost.MultiIndex as discussed in
		///       "Why You Should Use Boost MultiIndex (Part II)"
		class id_of_string_ref final {
		private:
			/// \brief Type alias for the type of the IDs
			using id_type = size_t;

			/// \brief Type alias for the map type
			using map_type = std::unordered_map<detail::ref_wrap_uom_wrap<const std::string>, id_type, detail::ref_wrap_hasher<const std::string>>;

			/// \brief The unordered_map that stores the string_ref-to-id lookup map
			map_type the_map;

		public:
			/// \brief A type alias for the const_iterator type as part of making this a range
			using const_iterator = map_type::const_iterator;

			id_of_string_ref() = default;

			/// \brief Insert a new string and return its new ID
			///
			/// Can be used if the name already exists
			inline const std::pair<const detail::ref_wrap_uom_wrap<const std::string>, id_type> & emplace(const string_cref &arg_string ///< The string to insert
			                                                                                              ) {
				return *( the_map.emplace( arg_string, the_map.size() ).first );
			}

			/// \brief Get the ID corresponding to the specified string
			inline id_type operator[](const string_cref &arg_string ///< The string to lookup
			                          ) const {
				return the_map.find( arg_string )->second;
			}

			/// \brief Return whether this contains the specified string_cref
			inline bool contains(const string_cref &arg_string ///< The string to lookup
			                     ) const {
				return ( the_map.find( arg_string ) != common::cend( the_map ) );
			}

			/// \brief Return whether this id_of_string_ref is empty
			inline bool empty() const {
				return the_map.empty();
			}

			/// \brief Return the number of strings that are stored
			inline size_t size() const {
				return the_map.size();
			}

			/// \brief Reserve space for the specified number of strings
			inline void reserve(const size_t &arg_count ///< The number of strings for which space should be reserved
			                    ) {
				the_map.reserve( arg_count );
			}

			/// \brief Clear the id_of_string_ref of all strings
			inline void clear() {
				the_map.clear();
			}

			/// \brief Standard const begin() operator to make id_of_string_ref into a range
			inline const_iterator begin() const {
				return common::cbegin( the_map );
			}

			/// \brief Standard const end() operator to make id_of_string_ref into a range
			inline const_iterator end() const {
				return common::cend( the_map );
			}
		};

	} // namespace common
} // namespace cath

#endif
