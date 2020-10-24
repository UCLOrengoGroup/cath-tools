/// \file
/// \brief The id_of_string header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_HPP

#include <boost/range/algorithm/find_if.hpp>

#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

#include <string>
#include <unordered_map>

namespace cath {
	namespace common {

		/// \brief Map strings to numeric IDs that count incrementally from 0
		///
		/// This makes lookups easier
		///
		/// This class isn't thread-safe
		///
		/// Not currently able to avoid constructing a new string if passed a char *
		///
		/// \TODO Consider adding (or changing this into) id_of_string_ref / id_of_string_view.
		///       The Boost string_ref currently has a std::hash specialisation
		///       but it's #if-ed out (#if 0).
		///
		/// \TODO Alternatively consider using Boost.MultiIndex as discussed in
		///       "Why You Should Use Boost MultiIndex (Part II)"
		class id_of_string final {
		public:
			/// \brief Type alias for the type of the IDs
			using id_type = size_t;

		private:
			/// \brief The unordered_map that stores the string-to-id lookup map
			std::unordered_map<std::string, id_type> the_map;

		public:
			/// \brief A const_iterator type alias as part of making this a range over pair<const string, id_type>
			using const_iterator = std::unordered_map<std::string, id_type>::const_iterator;

			/// \brief An iterator type alias that just duplicates const_iterator to appease some Boost code (Range?)
			using iterator = const_iterator;

			id_of_string() = default;

			/// \brief Insert a new string and return its new ID
			///
			/// Can be used if the name already exists
			inline const std::pair<const std::string, id_type> & emplace(std::string prm_string ///< The string to insert
			                                                             ) {
				return *( the_map.emplace( std::move( prm_string ), the_map.size() ).first );
			}

			/// \brief Get the ID corresponding to the specified string
			inline id_type operator[](const std::string &prm_string ///< The string to lookup
			                          ) const {
				return the_map.find( prm_string )->second;
			}

			/// \brief Return whether this id_of_string is empty
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

			/// \brief Clear the id_of_string of all strings
			inline void clear() {
				the_map.clear();
			}

			/// \brief Standard const begin() method, as part of making this into a range over pair<const string, id_type>
			inline const_iterator begin() const {
				return common::cbegin( the_map );
			}

			/// \brief Standard const end() method, as part of making this into a range over pair<const string, id_type>
			inline const_iterator end() const {
				return common::cend( the_map );
			}

		};

		/// \brief Find the string associated with the specified ID in the specified id_of_string
		///
		/// Note: This is inefficient! Don't use in performance-critical code
		inline const std::string & string_of_id(const id_of_string          &prm_id_of_string, ///< The id_of_string to query
		                                        const id_of_string::id_type &prm_id            ///< The ID of interest
		                                        ) {
			return boost::find_if(
				prm_id_of_string,
				[&] (const std::pair<const std::string, id_of_string::id_type> &x) {
					return ( x.second == prm_id );
				}
			)->first;
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_ID_OF_STRING_HPP
