/// \file
/// \brief The string parse tools header

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

#ifndef STRING_PARSE_TOOLS_H_INCLUDED
#define STRING_PARSE_TOOLS_H_INCLUDED

#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/c++14/cbegin_cend.h"
#include "common/debug_numeric_cast.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"

#include <string>

using namespace std::literals::string_literals;

namespace cath {
	/// \brief Type alias for boost::string_ref's const_iterator
	using str_ref_citr = boost::string_ref::const_iterator;

	namespace common {
		namespace detail {

			/// \brief Perform the actual spirit parse, throw if there's a problem and return the result
			///
			/// Note: please benchmark any changes to these functions to ensure they stay fast
			template <typename T, typename QiParse>
			inline T do_spirit_parse(str_citr         arg_begin_itr, ///< The iterator to the start of the strecth of string to parse (passed-by-value to allow efficient modification)
			                         const str_citr  &arg_end_itr,   ///< The iterator to the end of the strecth of string to parse
			                         QiParse        &&arg_qi_parse   ///< The boost::spirit parser
			                         ) {
				T value;
				const bool ok = boost::spirit::qi::parse(
					arg_begin_itr,
					arg_end_itr,
					std::forward<QiParse>( arg_qi_parse ),
					value
				);

				if ( ! ok || arg_begin_itr != arg_end_itr ) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception(
						"Failed to parse a number (of type "s
						+ typeid( T ).name()
						+ ") from "
						+ std::string{ arg_begin_itr, arg_end_itr }
					));
				}
				return value;
			}

		}

		/// \brief Find the iterator that points to (just before) the first non-whitespace character
		///        in the specified string_ref
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		///
		/// \todo Come C++17, switch from boost::string_ref to std::string_view
		inline str_ref_citr find_itr_before_first_non_space(const boost::string_ref &arg_substring ///< The string_ref in which to search 
		                                                    ) {
			const auto itr = boost::range::find_if(
				arg_substring,
				[] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); }
			);
			if ( itr == common::cend( arg_substring ) ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find any non-space chars in string"));
			}
			return itr;
		}

		/// \brief Find the iterator that points to (just after) the last non-whitespace character
		///        in the specified string_ref
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		///
		/// \todo Come C++17, switch from boost::string_ref to std::string_view
		inline str_ref_citr find_itr_after_last_non_space(const boost::string_ref &arg_substring ///< The string_ref in which to search 
		                                                  ) {
			const auto ritr = boost::range::find_if(
				arg_substring | boost::adaptors::reversed,
				[] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); }
			);
			if ( ritr.base() == common::cbegin( arg_substring ) ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find any non-space chars in string"));
			}
			return ritr.base();
		}

		/// \brief Return a string_ref to the section of the specified string_ref after trimming
		///
		/// This can sensibly be called with a string argument.
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		inline boost::string_ref dumb_trim_string_ref(const boost::string_ref &arg_substring ///< The string_ref to trim
		                                              ) {
			const auto before_first_non_space_itr = find_itr_before_first_non_space( arg_substring );
			const auto after_last_non_space_itr   = find_itr_after_last_non_space  ( arg_substring );
			return {
				&*before_first_non_space_itr,
				static_cast<size_t>( std::distance( before_first_non_space_itr, after_last_non_space_itr ) )
			};
		}

		/// \brief Parse a (possibly space-padded) double from the specified region of string
		///
		/// Note: please benchmark any changes to these functions to ensure they stay fast
		inline double parse_double_from_substring(const std::string &arg_string, ///< The string containing the region to parse
		                                          const size_t      &arg_start,  ///< The index of the start of the region to parse
		                                          const size_t      &arg_length  ///< The length of the region to parse
		                                          ) {
			const auto begin_itr = std::next( arg_string.begin(), static_cast<ptrdiff_t>( arg_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( arg_length ) );
			return detail::do_spirit_parse<double>(
				begin_itr,
				end_itr,
				   boost::spirit::omit[ *boost::spirit::qi::space ]
				>> boost::spirit::double_
				>> boost::spirit::omit[ *boost::spirit::qi::space ]
			);
		}

		/// \brief Parse a (possibly space-padded) double from the specified region of string
		///
		/// Note: please benchmark any changes to these functions to ensure they stay fast
		inline int parse_int_from_substring(const std::string &arg_string, ///< The string containing the region to parse
		                                    const size_t      &arg_start,  ///< The index of the start of the region to parse
		                                    const size_t      &arg_length  ///< The length of the region to parse
		                                    ) {
			const auto begin_itr = std::next( arg_string.begin(), static_cast<ptrdiff_t>( arg_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( arg_length ) );
			return detail::do_spirit_parse<int>(
				begin_itr,
				end_itr,
				   boost::spirit::omit[ *boost::spirit::qi::space ]
				>> boost::spirit::int_
				>> boost::spirit::omit[ *boost::spirit::qi::space ]
			);
		}

		/// \brief Parse a (possibly space-padded) double from the specified region of string
		///
		/// Note: please benchmark any changes to these functions to ensure they stay fast
		inline unsigned long int parse_ulong_from_substring(const std::string &arg_string, ///< The string containing the region to parse
		                                                    const size_t      &arg_start,  ///< The index of the start of the region to parse
		                                                    const size_t      &arg_length  ///< The length of the region to parse
		                                                    ) {
			const auto begin_itr = std::next( arg_string.begin(), static_cast<ptrdiff_t>( arg_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( arg_length ) );
			return detail::do_spirit_parse<unsigned long int>(
				begin_itr,
				end_itr,
				   boost::spirit::omit[ *boost::spirit::qi::space ]
				>> boost::spirit::ulong_
				>> boost::spirit::omit[ *boost::spirit::qi::space ]
			);
		}

	}
}

#endif
