/// \file
/// \brief The sub_string_parser class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_STRING_SUB_STRING_PARSER_H
#define _CATH_TOOLS_SOURCE_COMMON_STRING_SUB_STRING_PARSER_H

#include <string>

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		class sub_string_parser final {
		private:
			std::string the_substring;

		public:
			void reserve(const size_t &);
			double substr_as_double(const std::string &,
			                        const size_t &,
			                        const size_t &arg_count = std::string::npos);
			size_t substr_as_size_t(const std::string &,
			                        const size_t &,
			                        const size_t &arg_count = std::string::npos);
			int substr_as_int(const std::string &,
			                  const size_t &,
			                  const size_t &arg_count = std::string::npos);
			const std::string & substr_as_str_ref(const std::string &,
			                                      const size_t &,
			                                      const size_t &arg_count = std::string::npos);
		};

		/// \brief TODOCUMENT
		inline void sub_string_parser::reserve(const size_t &arg_size ///< TODOCUMENT
		                                       ) {
			the_substring.reserve( arg_size );
		}

		/// \brief TODOCUMENT
		inline double sub_string_parser::substr_as_double(const std::string &arg_source_string, ///< TODOCUMENT
		                                                  const size_t      &arg_from_pos,      ///< TODOCUMENT
		                                                  const size_t      &arg_count          ///< TODOCUMENT
		                                                  ) {
			the_substring.assign( arg_source_string, arg_from_pos, arg_count );
			return std::stod( the_substring );
		}

		/// \brief TODOCUMENT
		inline size_t sub_string_parser::substr_as_size_t(const std::string &arg_source_string, ///< TODOCUMENT
		                                                  const size_t      &arg_from_pos,      ///< TODOCUMENT
		                                                  const size_t      &arg_count          ///< TODOCUMENT
		                                                  ) {
			the_substring.assign( arg_source_string, arg_from_pos, arg_count );
			return std::stoul( the_substring );
		}

		/// \brief TODOCUMENT
		inline int sub_string_parser::substr_as_int(const std::string &arg_source_string, ///< TODOCUMENT
		                                            const size_t      &arg_from_pos,      ///< TODOCUMENT
		                                            const size_t      &arg_count          ///< TODOCUMENT
		                                            ) {
			the_substring.assign( arg_source_string, arg_from_pos, arg_count );
			return std::stoi( the_substring );
		}

		/// \brief TODOCUMENT
		inline const std::string & sub_string_parser::substr_as_str_ref(const std::string &arg_source_string, ///< TODOCUMENT
		                                                                const size_t      &arg_from_pos,      ///< TODOCUMENT
		                                                                const size_t      &arg_count          ///< TODOCUMENT
		                                                                ) {
			the_substring.assign( arg_source_string, arg_from_pos, arg_count );
			return the_substring;
		}

	} // namespace common
} // namespace cath

#endif
