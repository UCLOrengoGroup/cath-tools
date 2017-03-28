/// \file
/// \brief The rapidjson_writer header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_RAPIDJSON_ADDENDA_RAPIDJSON_WRITER_H
#define _CATH_TOOLS_SOURCE_COMMON_RAPIDJSON_ADDENDA_RAPIDJSON_WRITER_H

#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include "common/json_style.hpp"

#include <string>

namespace cath {
	namespace common {

		/// \brief Wrap rapidjson's Writer/PrettyWriter simple writing functionality
		template <json_style Style = json_style::PRETTY>
		class rapidjson_writer final {
		private:
			/// \brief Type alias for the type of buffer to which the JSON should be written
			using buffer_t = rapidjson::StringBuffer;

			/// \brief Type alias for the type of rapidjson writer
			using writer_t = std::conditional_t< Style == json_style::PRETTY,
			                                     rapidjson::PrettyWriter< buffer_t >,
			                                     rapidjson::Writer      < buffer_t > >;

			/// \brief The buffer to which the JSON should be written
			buffer_t buffer;

			/// \brief The rapidjson writer
			writer_t writer{ buffer };

		public:
			rapidjson_writer() = default;

			/// \brief Write a raw string in the JSON
			template <typename... Ts>
			rapidjson_writer & write_raw_string(const std::string &arg_string
			                                    ) {
				writer.RawValue( arg_string.c_str(), arg_string.length(), rapidjson::kNullType );
				return *this;
			}

			/// \brief Write a key to the JSON
			template <typename... Ts>
			rapidjson_writer & write_key(Ts &&...arg_values ///< The key string (can be std::string or char *)
			                             ) {
				writer.Key( std::forward<Ts>( arg_values )... );
				return *this;
			}

			/// \brief Write a string value to the JSON
			template <typename... Ts>
			rapidjson_writer & write_string(Ts &&...arg_values ///< The string value (can be std::string or char *)
			                                ) {
				writer.String( std::forward<Ts>( arg_values )... );
				return *this;
			}

			/// \brief Write a bool value to the JSON
			rapidjson_writer & write_bool(const bool &arg_value ///< The bool value to write to the JSON
			                              ) {
				writer.Bool( arg_value );
				return *this;
			}

			/// \brief Write a double value to the JSON
			rapidjson_writer & write_double(const double &arg_value ///< The double value to write to the JSON
			                                ) {
				writer.Double( arg_value );
				return *this;
			}

			/// \brief Write a int value to the JSON
			rapidjson_writer & write_int(const int &arg_value ///< The int value to write to the JSON
			                             ) {
				writer.Int( arg_value );
				return *this;
			}

			/// \brief Write a int64_t value to the JSON
			rapidjson_writer & write_int64(const int64_t &arg_value ///< The int64_t value to write to the JSON
			                               ) {
				writer.Int64( arg_value );
				return *this;
			}

			/// \brief Write a uint value to the JSON
			rapidjson_writer & write_uint(const uint &arg_value ///< The uint value to write to the JSON
			                              ) {
				writer.Uint( arg_value );
				return *this;
			}

			/// \brief Write a uint64_t value to the JSON
			rapidjson_writer & write_uint64(const uint64_t &arg_value ///< The uint64_t value to write to the JSON
			                                ) {
				writer.Uint64( arg_value );
				return *this;
			}

			/// \brief Write a null value to the JSON
			rapidjson_writer & write_null() {
				writer.Null();
				return *this;
			}

			/// \brief Write the start of an array to the JSON
			rapidjson_writer & start_array() {
				writer.StartArray();
				return *this;
			}

			/// \brief Write the end of an array to the JSON
			rapidjson_writer & end_array() {
				writer.EndArray();
				return *this;
			}

			/// \brief Write the start of an object to the JSON
			rapidjson_writer & start_object() {
				writer.StartObject();
				return *this;
			}

			/// \brief Write the end of an object to the JSON
			rapidjson_writer & end_object() {
				writer.EndObject();
				return *this;
			}

			/// \brief Get a C-style string of the JSON
			const char * get_c_string() const {
				return buffer.GetString();
			}

			/// \brief Get a std::string of the JSON
			std::string get_cpp_string() const {
				return { get_c_string() };
			}

			/// \brief Get whether the JSON is complete (ie finished the initial array/object/value)
			///        after which no more data can be written
			bool is_complete() const {
				return writer.IsComplete();
			}
		};

	} // namespace common
} // namespace cath

#endif
