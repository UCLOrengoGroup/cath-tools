/// \file
/// \brief The rapidjson_writer_list header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_RAPIDJSON_ADDENDA_RAPIDJSON_WRITER_LIST_H
#define _CATH_TOOLS_SOURCE_COMMON_RAPIDJSON_ADDENDA_RAPIDJSON_WRITER_LIST_H

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/range/adaptor/indirected.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/rapidjson_addenda/rapidjson_writer.hpp"
#include "common/type_aliases.hpp"

#include <string>
// #include <vector>

namespace cath {
	namespace common {


		/// \brief A list of rapidjson_writers
		template <json_style Style = json_style::PRETTY,
		          typename   OStrm = rapidjson::StringBuffer>
		class rapidjson_writer_list final {
		private:
			/// \brief Type alias for the type of the rapidjson_writer to store
			using rapidjson_writer_t = rapidjson_writer<Style, OStrm>;

			/// \brief Type alias for a vector of unique_ptr of rapidjson_writer_t
			using rapidjson_writer_t_uptr_vec = uptr_vec<rapidjson_writer_t>;

			/// \brief Type alias for the type of buffer to which the JSON should be written
			uptr_vec<rapidjson_writer_t> rapidjson_writers;

			/// \brief An iterator type alias as part of making this a non-const range
			///
			/// Note that this pipes through boost::indirected_range so the range is
			/// over elements of type `rapidjson_writer_t &` rather than `unique_ptr<rapidjson_writer_t>`
			using iterator = common::range_iterator_t< boost::indirected_range<rapidjson_writer_t_uptr_vec> >;

			/// \brief Standard non-const begin() method, as part of making this a range over rapidjson_writers
			///
			/// Note that this pipes through boost::indirected_range so the range is
			/// over elements of type `rapidjson_writer_t &` rather than `unique_ptr<rapidjson_writer_t>`
			auto begin() -> iterator {
				return std::begin( rapidjson_writers | boost::adaptors::indirected );
			}

			/// \brief Standard non-const end() method, as part of making this a range over rapidjson_writers
			///
			/// Note that this pipes through boost::indirected_range so the range is
			/// over elements of type `rapidjson_writer_t &` rather than `unique_ptr<rapidjson_writer_t>`
			auto end() -> iterator {
				return std::end( rapidjson_writers | boost::adaptors::indirected );
			}


		public:
			rapidjson_writer_list() = default;

			/// \brief Reserve space for the specified number of rapidjson_writers
			inline void reserve(const size_t &arg_size ///< The number of rapidjson_writers for which space should be reserved
			                    ) {
				return rapidjson_writers.reserve( arg_size );
			}

			/// \brief Get the number of rapidjson_writers currently stored
			inline size_t size() const {
				return rapidjson_writers.size();
			}
			/// \brief Get whether the list of rapidjson_writers is currently empty
			inline bool empty() const {
				return rapidjson_writers.empty();
			}

			/// \brief Emplace back a new rapidjson_writer with the specified arguments (which are perfect-forwarded to the ctor)
			template <typename... Ts>
			inline void emplace_back(Ts &&...args ///< The arguments
			                         ) {
				rapidjson_writers.emplace_back( std::make_unique<rapidjson_writer_t>( std::forward<Ts>( args )... ) );
			}

			/// \brief A const_iterator type alias as part of making this a const range
			///
			/// Note that this pipes through boost::indirected_range so the range is
			/// over elements of type `const rapidjson_writer_t &` rather than `unique_ptr<rapidjson_writer_t>`
			using const_iterator = common::range_const_iterator_t< boost::indirected_range<const rapidjson_writer_t_uptr_vec> >;

			/// \brief Standard const begin() method, as part of making this a range over rapidjson_writers
			///
			/// Note that this pipes through boost::indirected_range so the range is
			/// over elements of type `rapidjson_writer_t &` rather than `unique_ptr<rapidjson_writer_t>`
			auto begin() const -> const_iterator {
				return common::cbegin( rapidjson_writers | boost::adaptors::indirected );
			}

			/// \brief Standard const end() method, as part of making this a range over rapidjson_writers
			///
			/// Note that this pipes through boost::indirected_range so the range is
			/// over elements of type `rapidjson_writer_t &` rather than `unique_ptr<rapidjson_writer_t>`
			auto end() const -> const_iterator {
				return common::cend( rapidjson_writers | boost::adaptors::indirected );
			}

			/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
			rapidjson_writer_list            (const rapidjson_writer_list &) = delete;
			/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
			rapidjson_writer_list            (rapidjson_writer_list &&     ) = delete;
			/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
			rapidjson_writer_list & operator=(const rapidjson_writer_list &) = delete;
			/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
			rapidjson_writer_list & operator=(rapidjson_writer_list &&     ) = delete;

			/// \brief Explicit ctor for when 1 or more arguments are specified, perfect-forwards all parameters to the outstream ctor
			template <typename T>
			explicit rapidjson_writer_list(const std::vector<T> &arg_values ///< The values with which to build a list of corresponding rapidjson_writers
			                               ) {
				rapidjson_writers.reserve( arg_values.size() );
				for (const T &value : arg_values) {
					emplace_back( value );
				}
			}

			/// \brief Write a raw string in the JSON
			template <typename... Ts>
			rapidjson_writer_list & write_raw_string(const std::string &arg_string ///< The string to write to the rapidjson_writers
			                                         ) {
				for (auto &writer : *this) {
					writer.write_raw_string( arg_string );
				}
				return *this;
			}

			/// \brief Write a key to the JSON
			rapidjson_writer_list & write_key(const std::string &arg_key ///< The key string
			                                  ) {
				for (auto &writer : *this) {
					writer.write_key( arg_key );
				}
				return *this;
			}

			/// \brief Write a key to the JSON
			rapidjson_writer_list & write_key(const char * const arg_key ///< The key string (can be std::string or char *)
			                                  ) {
				for (auto &writer : *this) {
					writer.write_key( arg_key );
				}
				return *this;
			}

			/// \brief Write a string value to the JSON
			rapidjson_writer_list & write_value(const char * const arg_value ///< The string value (can be std::string or char *)
			                                    ) {
				for (auto &writer : *this) {
					writer.write_value( arg_value );
				}
				return *this;
			}

			/// \brief Write a bool value to the JSON
			template <typename T>
			rapidjson_writer_list & write_value(const T &arg_value ///< The bool value to write to the JSON
			                                    ) {
				for (auto &writer : *this) {
					writer.write_value( arg_value );
				}
				return *this;
			}

			/// \brief Write a null value to the JSON
			rapidjson_writer_list & write_null() {
				for (auto &writer : *this) {
					writer.write_null();
				}
				return *this;
			}

			/// \brief Write the start of an array to the JSON
			rapidjson_writer_list & start_array() {
				for (auto &writer : *this) {
					writer.start_array();
				}
				return *this;
			}

			/// \brief Write the end of an array to the JSON
			rapidjson_writer_list & end_array() {
				for (auto &writer : *this) {
					writer.end_array();
				}
				return *this;
			}

			/// \brief Write the start of an object to the JSON
			rapidjson_writer_list & start_object() {
				for (auto &writer : *this) {
					writer.start_object();
				}
				return *this;
			}

			/// \brief Write the end of an object to the JSON
			rapidjson_writer_list & end_object() {
				for (auto &writer : *this) {
					writer.end_object();
				}
				return *this;
			}

			/// \brief Get a C-style string of the JSON
			const char * get_c_string(const size_t &arg_index ///< The index of the writer from which to grab the string
			                          ) const {
				return rapidjson_writers[ arg_index ]->get_c_string();
			}

			/// \brief Get a std::string of the JSON
			std::string get_cpp_string(const size_t &arg_index ///< The index of the writer from which to grab the string
			                           ) const {
				return rapidjson_writers[ arg_index ]->get_cpp_string();
			}

			/// \brief Get whether the JSON is complete (ie finished the initial array/object/value)
			///        after which no more data can be written
			bool is_complete() const {
				return boost::algorithm::all_of(
					*this,
					[] (const rapidjson_writer<Style, OStrm> &x) { return x.is_complete(); }
				);
			}

		};

		/// \brief Write a key and a value
		template <typename T, json_style Style, typename OStrm>
		rapidjson_writer_list<Style, OStrm> & write_key_value(rapidjson_writer_list<Style, OStrm> &arg_writer_list, ///< The rapidjson_writer_list to which they key and value should be written
		                                                      const std::string                   &arg_key,         ///< The key to write
		                                                      const T                             &arg_value        ///< The value to write
		                                                      ) {
			return arg_writer_list.write_key( arg_key ).write_value( arg_value );
		}

	} // namespace common
} // namespace cath

#endif
