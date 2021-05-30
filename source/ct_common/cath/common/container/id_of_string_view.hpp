/// \file
/// \brief The id_of_string_view header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_VIEW_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_VIEW_HPP

#include <optional>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "cath/common/optional/make_optional_if.hpp"

namespace cath::common {
	namespace detail {

		/// \brief Hasher for ::std::string_view
		struct string_view_hasher final {
			/// \brief The function operator that hashes the character range referred to by the string_view
			size_t operator()(const ::std::string_view &prm_value ///< The string_view value to hash
			                  ) const {
				return boost::hash_range( prm_value.begin(), prm_value.end() );
			}
		};

	} // namespace detail

	/// \brief Map string_view to numeric IDs that count incrementally from 0
	///
	/// This class isn't thread-safe
	///
	/// \TODO Consider adding (or changing this into) id_of_string_view_ref / id_of_string_view_view.
	///       The Boost string_view currently has a std::hash specialisation
	///       but it's #if-ed out (#if 0).
	///
	/// \TODO Alternatively consider using Boost.MultiIndex as discussed in
	///       "Why You Should Use Boost MultiIndex (Part II)"
	class id_of_string_view final {
	private:
		/// \brief Type alias for the type of the IDs
		using id_type = size_t;

		/// \brief Type alias for the map type
		using map_type = std::unordered_map<::std::string_view, id_type, detail::string_view_hasher>;

		/// \brief The unordered_map that stores the string_view-to-id lookup map
		map_type the_map;

	public:
		/// \brief A type alias for the const_iterator type as part of making this a range
		using const_iterator = map_type::const_iterator;

		id_of_string_view() = default;

		/// \brief Insert a new string and return its new ID
		///
		/// Can be used if the name already exists
		inline const std::pair<const ::std::string_view, id_type> & emplace(const ::std::string_view &prm_string ///< The string to insert
		                                                                    ) {
			return *( the_map.emplace( prm_string, the_map.size() ).first );
		}

		/// \brief Get the ID corresponding to the specified string
		inline ::std::optional<id_type> operator[](const ::std::string_view &prm_string ///< The string to lookup
		                                           ) const {
			const auto find_itr = the_map.find( prm_string );
			return if_then_optional(
				find_itr != ::std::cend( the_map ),
				find_itr->second
			);
		}

		/// \brief Return whether this contains the specified ::std::string_view
		inline bool contains(const ::std::string_view &prm_string ///< The string to lookup
		                     ) const {
			return ( the_map.find( prm_string ) != ::std::cend( the_map ) );
		}

		/// \brief Return whether this id_of_string_view is empty
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

		/// \brief Clear the id_of_string_view of all strings
		inline void clear() {
			the_map.clear();
		}

		/// \brief Standard const begin() operator to make id_of_string_view into a range
		inline const_iterator begin() const {
			return ::std::cbegin( the_map );
		}

		/// \brief Standard const end() operator to make id_of_string_view into a range
		inline const_iterator end() const {
			return ::std::cend( the_map );
		}
	};

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_VIEW_HPP
