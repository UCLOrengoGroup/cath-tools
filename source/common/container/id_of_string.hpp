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

#ifndef _CATH_TOOLS_SOURCE_COMMON_CONTAINER_ID_OF_STRING_H
#define _CATH_TOOLS_SOURCE_COMMON_CONTAINER_ID_OF_STRING_H

#include "common/cpp14/cbegin_cend.hpp"
#include "exception/invalid_argument_exception.hpp"

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
		class id_of_string final {
		private:
			/// \brief Type alias for the type of the IDs
			using id_type = size_t;

			/// \brief The unordered_map that stores the string-to-id lookup map
			std::unordered_map<std::string, id_type> the_map;

		public:
			id_of_string() = default;

			/// \brief Insert a new string and return its new ID
			inline id_type emplace(std::string arg_string ///< The string to insert
			                       ) {
#ifndef NDEBUG
				if ( the_map.count( arg_string ) > 0 ) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use id_of_string::emplace() to insert a string that's already present"));
				}
#endif
				const auto prev_size = the_map.size();
				the_map.emplace( std::move( arg_string ), prev_size );
				return prev_size;
			}

			/// \brief Insert a new string if it isn't already present, and return its ID either way
			inline id_type emplace_if_not_present(std::string arg_string ///< The string to insert
			                                      ) {
				const auto find_itr = the_map.find( arg_string );
				return ( find_itr == common::cend( the_map ) )
					? emplace( std::move( arg_string ) )
					: find_itr->second;
			}

			/// \brief Get the ID corresponding to the specified string
			inline id_type get(const std::string &arg_string ///< The string to lookup
			                   ) const {
				return the_map.at( arg_string );
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
			inline void reserve(const size_t &arg_count ///< The number of strings for which space should be reserved
			                    ) {
				the_map.reserve( arg_count );
			}

			/// \brief Clear the id_of_string of all strings
			inline void clear() {
				the_map.clear();
			}
		};

	} // namespace common
} // namespace cath

#endif
