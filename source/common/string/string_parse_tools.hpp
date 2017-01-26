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

#ifndef _CATH_TOOLS_SOURCE_COMMON_STRING_STRING_PARSE_TOOLS_H
#define _CATH_TOOLS_SOURCE_COMMON_STRING_STRING_PARSE_TOOLS_H

#include <boost/core/ignore_unused.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/boost_addenda/make_string_ref.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/runtime_error_exception.hpp"

#include <string>

using namespace std::literals::string_literals;

namespace cath {

	/// \brief Type alias for boost::string_ref's const_iterator
	using str_ref_citr = boost::string_ref::const_iterator;

	/// \brief A type alias for a pair of str_ctirs
	///
	/// \todo Should this be moved into common/type_aliases.hpp?
	using str_citr_str_citr_pair = std::pair<str_citr, str_citr>;

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

		} // namespace detail

		/// \brief Return an iterator pointing to the first point before a non-space character in the region between the specified
		///        iterators (or arg_end if none is found)
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		inline str_citr find_itr_before_first_non_space(const str_citr &arg_begin, ///< A  begin              iterator of the region of string to search
		                                                const str_citr &arg_end    ///< An end (one-past-end) iterator of the region of string to search
		                                                ) {
			return std::find_if(
				arg_begin,
				arg_end,
				[] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); }
			);
		}

		/// \brief Return an iterator pointing to the first point before a space character in the region between the specified
		///        iterators (or arg_end if none is found)
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		inline str_citr find_itr_before_first_space(const str_citr &arg_begin, ///< A  begin              iterator of the region of string to search
		                                            const str_citr &arg_end    ///< An end (one-past-end) iterator of the region of string to search
		                                            ) {
			return std::find_if(
				arg_begin,
				arg_end,
				[] (const auto &x) { return ( ( x == ' ' ) || ( x == '\t' ) ); }
			);
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
		///
		/// \todo Fix: this given a string of all whitespace chars, this would currently cross over the begin/end pointers
		///       so should maybe calculated end pointer first and then restrict calculation for start pointer within
		///       range up to that end.
		inline boost::string_ref dumb_trim_string_ref(const boost::string_ref &arg_substring ///< The string_ref to trim
		                                              ) {
			return make_string_ref(
				find_itr_before_first_non_space( arg_substring ),
				find_itr_after_last_non_space  ( arg_substring )
			);
		}

		/// \brief Return a pair of offsets to the section of the specified string_ref after trimming
		///
		/// This can sensibly be called with a string argument.
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		///
		/// \todo Fix: this given a string of all whitespace chars, this would currently cross over the begin/end pointers
		///       so should maybe calculated end pointer first and then restrict calculation for start pointer within
		///       range up to that end.
		template <typename T>
		std::pair<T, T> dumb_trim_string_ref_to_offsets(const boost::string_ref &arg_substring ///< The string_ref to trim
		                                                ) {
			return {
				static_cast<T>( std::distance( common::cbegin( arg_substring ), find_itr_before_first_non_space( arg_substring ) ) ),
				static_cast<T>( std::distance( common::cbegin( arg_substring ), find_itr_after_last_non_space  ( arg_substring ) ) )
			};
		}

		/// \brief Parse a double from the field between the two specified string iterators
		inline double parse_double_from_field(const str_citr &arg_begin_itr, ///< A const_iterator pointing to the begin              of the field to be parsed
		                                      const str_citr &arg_end_itr    ///< A const_iterator pointing to the end (one-past-end) of the field to be parsed
		                                      ) {
			return detail::do_spirit_parse<double>(
				arg_begin_itr,
				arg_end_itr,
				boost::spirit::double_
			);
		}

		/// \brief Parse a float from the field between the two specified string iterators
		inline float parse_float_from_field(const str_citr &arg_begin_itr, ///< A const_iterator pointing to the begin              of the field to be parsed
		                                    const str_citr &arg_end_itr    ///< A const_iterator pointing to the end (one-past-end) of the field to be parsed
		                                    ) {
			return detail::do_spirit_parse<float>(
				arg_begin_itr,
				arg_end_itr,
				boost::spirit::float_
			);
		}

		/// \brief Parse an unsigned int from the field between the two specified string iterators
		inline unsigned int parse_uint_from_field(const str_citr &arg_begin_itr, ///< A const_iterator pointing to the begin              of the field to be parsed
		                                          const str_citr &arg_end_itr    ///< A const_iterator pointing to the end (one-past-end) of the field to be parsed
		                                          ) {
			return detail::do_spirit_parse<unsigned int>(
				arg_begin_itr,
				arg_end_itr,
				boost::spirit::uint_
			);
		}

		/// \brief Parse a (possibly space-padded) float from the specified region of string
		///
		/// Note: please benchmark any changes to these functions to ensure they stay fast
		inline float parse_float_from_substring(const std::string &arg_string, ///< The string containing the region to parse
		                                        const size_t      &arg_start,  ///< The index of the start of the region to parse
		                                        const size_t      &arg_length  ///< The length of the region to parse
		                                        ) {
			const auto begin_itr = std::next( arg_string.begin(), static_cast<ptrdiff_t>( arg_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( arg_length ) );
			return detail::do_spirit_parse<float>(
				begin_itr,
				end_itr,
				   boost::spirit::omit[ *boost::spirit::qi::space ]
				>> boost::spirit::float_
				>> boost::spirit::omit[ *boost::spirit::qi::space ]
			);
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

		/// \brief Parse a (possibly space-padded) int from the specified region of string
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

		/// \brief Parse a (possibly space-padded) unsigned int from the specified region of string
		///
		/// Note: please benchmark any changes to these functions to ensure they stay fast
		inline unsigned int parse_uint_from_substring(const std::string &arg_string, ///< The string containing the region to parse
		                                              const size_t      &arg_start,  ///< The index of the start of the region to parse
		                                              const size_t      &arg_length  ///< The length of the region to parse
		                                              ) {
			const auto begin_itr = std::next( arg_string.begin(), static_cast<ptrdiff_t>( arg_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( arg_length ) );
			return detail::do_spirit_parse<unsigned int>(
				begin_itr,
				end_itr,
				   boost::spirit::omit[ *boost::spirit::qi::space ]
				>> boost::spirit::uint_
				>> boost::spirit::omit[ *boost::spirit::qi::space ]
			);
		}

		/// \brief Parse a (possibly space-padded) unsigned long from the specified region of string
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

		/// \brief Find the iterators wrapping the specified field in the specified in the specified string
		///        starting from the specified initial iterator at the specified index
		inline str_citr_str_citr_pair find_field_itrs(const std::string &arg_string,      ///< The string to search
		                                              const size_t      &arg_field_index, ///< The index of the field to find
		                                              const size_t      &arg_init_index,  ///< The index of the field from which the search should start
		                                              const str_citr    &arg_init_itr     ///< The iterator from which the search should start
		                                              ) {
			const auto end_itr = common::cend( arg_string );
			auto field_itr = find_itr_before_first_non_space( arg_init_itr, end_itr );
			for (const size_t field_ctr : boost::irange( arg_init_index, arg_field_index ) ) {
				boost::ignore_unused( field_ctr );
				field_itr = find_itr_before_first_space    ( field_itr, end_itr );
				field_itr = find_itr_before_first_non_space( field_itr, end_itr );
				if ( field_itr == end_itr ) {
					BOOST_THROW_EXCEPTION(runtime_error_exception(
						"Unable to find field "
						+ std::to_string( arg_field_index )
						+ " in line \""
						+ ( arg_string.size() > 103 ? ( arg_string.substr( 0, 100 ) + "[...]" ) : arg_string )
						+ "\""
					));
				}
			}
			return {
				field_itr,
				find_itr_before_first_space( field_itr, end_itr )
			};
		}

		/// \brief Find the iterators wrapping the specified field in the specified string
		inline str_citr_str_citr_pair find_field_itrs(const std::string &arg_string,      ///< The string to search
		                                              const size_t      &arg_field_index  ///< The index of the field to find
		                                              ) {
			return find_field_itrs(
				arg_string,
				arg_field_index,
				0,
				common::cbegin( arg_string )
			);
		}

		/// \brief Return an array<char, N> populated with N of the chars of the specified range of chars,
		///        (filling with 0s if the string isn't long enough)
		template <size_t N, typename Itr>
		inline std::array<char, N> get_char_arr_of_char_range(const Itr &arg_begin_itr, ///< The start of the range of characters to copy
		                                                      const Itr &arg_end_itr    ///< The end of the range of characters to copy
		                                                      ) {
			std::array<char, N> result;
			const auto end_at_n_itr = std::next( arg_begin_itr, N );
			const auto result_copied_itr = std::copy(
				arg_begin_itr,
				std::min( arg_end_itr, end_at_n_itr ),
				std::begin( result )
			);
			std::fill( result_copied_itr, std::end( result ), 0 );
			return result;
		}

		/// \brief Return an array<char, N> populated with N of the chars of the specified string starting from the
		///        specified index, (filling with 0s if the string isn't long enough)
		template <size_t N>
		inline std::array<char, N> get_char_arr_of_substring(const std::string &arg_string,     ///< The string from which to copy the characters
		                                                     const size_t      &arg_begin_index ///< The index at which to start copying characters from the string
		                                                     ) {
			const auto start_itr = std::next( common::cbegin( arg_string ), static_cast<ptrdiff_t>( arg_begin_index ) );
			const auto end_itr   = common::cend( arg_string );
			return get_char_arr_of_char_range<N>(
				std::min( start_itr, end_itr ),
				end_itr
			);
		}

	} // namespace common
} // namespace cath

#endif
