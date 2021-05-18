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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_REF_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_REF_HPP

#include "cath/common/container/detail/ref_wrap_hasher.hpp"
#include "cath/common/container/detail/ref_wrap_uom_wrap.hpp"
#include "cath/common/type_aliases.hpp"

#include <string>
#include <unordered_map>

namespace cath {
	namespace common {

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
			inline const std::pair<const detail::ref_wrap_uom_wrap<const std::string>, id_type> & emplace(const string_cref &prm_string ///< The string to insert
			                                                                                              ) {
				return *( the_map.emplace( prm_string, the_map.size() ).first );
			}

			/// \brief Get the ID corresponding to the specified string
			inline id_type operator[](const string_cref &prm_string ///< The string to lookup
			                          ) const {
				return the_map.find( prm_string )->second;
			}

			/// \brief Return whether this contains the specified string_cref
			inline bool contains(const string_cref &prm_string ///< The string to lookup
			                     ) const {
				return ( the_map.find( prm_string ) != ::std::cend( the_map ) );
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
			inline void reserve(const size_t &prm_count ///< The number of strings for which space should be reserved
			                    ) {
				the_map.reserve( prm_count );
			}

			/// \brief Clear the id_of_string_ref of all strings
			inline void clear() {
				the_map.clear();
			}

			/// \brief Standard const begin() operator to make id_of_string_ref into a range
			inline const_iterator begin() const {
				return ::std::cbegin( the_map );
			}

			/// \brief Standard const end() operator to make id_of_string_ref into a range
			inline const_iterator end() const {
				return ::std::cend( the_map );
			}
		};

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_REF_HPP
