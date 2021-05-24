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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_RAPIDJSON_WRITER_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_RAPIDJSON_WRITER_HPP

#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/json_style.hpp"

#include <string>

namespace cath::common {

	/// \brief Wrap rapidjson's Writer/PrettyWriter simple writing functionality
	///
	/// \todo The interface could be made a bit cleaner by many of the write_... methods
	///       being renamed to the overloads of the same name.
	template <json_style Style = json_style::PRETTY,
	          typename   OStrm = rapidjson::StringBuffer>
	class rapidjson_writer final {
	private:
		// /// \brief Type alias for the type of buffer to which the JSON should be written
		// using outstream_t = rapidjson::StringBuffer;

		/// \brief Type alias for the type of rapidjson writer
		using writer_t    = std::conditional_t< Style == json_style::PRETTY,
		                                        rapidjson::PrettyWriter< OStrm >,
		                                        rapidjson::Writer      < OStrm > >;

		/// \brief The buffer to which the JSON should be written
		OStrm outstream;

		/// \brief The rapidjson writer
		writer_t writer{ outstream };

	public:
		rapidjson_writer() = default;

		/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
		rapidjson_writer            (const rapidjson_writer &) = delete;
		/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
		rapidjson_writer            (rapidjson_writer &&     ) = delete;
		/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
		rapidjson_writer & operator=(const rapidjson_writer &) = delete;
		/// \brief Make this explicitly non-copyable/non-movable since the member rapidjson writer/outstream types are
		rapidjson_writer & operator=(rapidjson_writer &&     ) = delete;

		/// \brief Explicit ctor for when 1 or more arguments are specified, perfect-forwards all parameters to the outstream ctor
		template <typename T, typename... Ts >
		explicit rapidjson_writer(T  &&   prm_value, ///< The first parameter
		                          Ts &&...prm_values ///< Any further parameters
		                          ) : outstream{
		                              	std::forward< T  >( prm_value  ),
		                              	std::forward< Ts >( prm_values )...
		                              }  {
		}

		/// \brief Write a raw string in the JSON
		template <typename... Ts>
		rapidjson_writer & write_raw_string(const std::string &prm_string
		                                    ) {
			writer.RawValue( prm_string.c_str(), prm_string.length(), rapidjson::kNullType );
			return *this;
		}

		/// \brief Write a key to the JSON
		rapidjson_writer & write_key(const std::string &prm_key ///< The key string
		                             ) {
			writer.Key( prm_key.c_str(), debug_numeric_cast<unsigned int>( prm_key.length() ) );
			return *this;
		}

		/// \brief Write a key to the JSON
		rapidjson_writer & write_key(const char * const prm_key ///< The key string (can be std::string or char *)
		                             ) {
			writer.Key( prm_key );
			return *this;
		}

		/// \brief Write a string value to the JSON
		rapidjson_writer & write_value(const std::string &prm_value ///< The string value (can be std::string or char *)
		                               ) {
			writer.String( prm_value );
			return *this;
		}

		/// \brief Write a string value to the JSON
		rapidjson_writer & write_value(const char * const prm_value ///< The string value (can be std::string or char *)
		                               ) {
			writer.String( prm_value );
			return *this;
		}

		/// \brief Write a bool value to the JSON
		rapidjson_writer & write_value(const bool &prm_value ///< The bool value to write to the JSON
		                               ) {
			writer.Bool( prm_value );
			return *this;
		}

		/// \brief Write a double value to the JSON
		rapidjson_writer & write_value(const double &prm_value ///< The double value to write to the JSON
		                               ) {
			writer.Double( prm_value );
			return *this;
		}

		/// \brief Write a int value to the JSON
		rapidjson_writer & write_value(const int &prm_value ///< The int value to write to the JSON
		                               ) {
			writer.Int( prm_value );
			return *this;
		}

		/// \brief Write a int64_t value to the JSON
		rapidjson_writer & write_value(const int64_t &prm_value ///< The int64_t value to write to the JSON
		                               ) {
			writer.Int64( prm_value );
			return *this;
		}

		/// \brief Write a uint value to the JSON
		rapidjson_writer & write_value(const uint32_t &prm_value ///< The uint value to write to the JSON
		                               ) {
			writer.Uint( prm_value );
			return *this;
		}

		/// \brief Write a uint64_t value to the JSON
		rapidjson_writer & write_value(const uint64_t &prm_value ///< The uint64_t value to write to the JSON
		                               ) {
			writer.Uint64( prm_value );
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
		[[nodiscard]] const char *get_c_string() const {
			return outstream.GetString();
		}

		/// \brief Get a std::string of the JSON
		[[nodiscard]] std::string get_cpp_string() const {
			return { get_c_string() };
		}

		/// \brief Get whether the JSON is complete (ie finished the initial array/object/value)
		///        after which no more data can be written
		[[nodiscard]] bool is_complete() const {
			return writer.IsComplete();
		}

		/// \brief Write a key and a value
		template <typename   T>
		rapidjson_writer & write_key_value(const std::string              &prm_key,    ///< The key to write
		                                   const T                        &prm_value   ///< The value to write
		                                   ) {
			return write_key( prm_key ).write_value( prm_value );
		}
	};


} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_RAPIDJSON_WRITER_HPP
