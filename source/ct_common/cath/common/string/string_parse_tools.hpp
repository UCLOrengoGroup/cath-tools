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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_STRING_PARSE_TOOLS_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_STRING_PARSE_TOOLS_HPP

#include <string>
#include <string_view>

#include <boost/core/demangle.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/utility/string_ref.hpp>

#include "cath/common/boost_addenda/make_string_ref.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::std::literals::string_literals;

namespace cath {

	/// \brief Type alias for boost::string_ref's const_iterator
	using str_ref_citr = boost::string_ref::const_iterator;

	/// \brief Type alias for std::string_view's const_iterator
	using str_view_citr = ::std::string_view::const_iterator;

	namespace common {
		namespace detail {

			/// \brief Perform the actual spirit parse, throw if there's a problem and return the result
			///
			/// Note: please benchmark any changes to these functions to ensure they stay fast
			template <typename T, typename QiParse>
			inline T do_spirit_parse(str_citr         prm_begin_itr, ///< The iterator to the start of the stretch of string to parse (passed-by-value to allow efficient modification)
			                         const str_citr  &prm_end_itr,   ///< The iterator to the end of the stretch of string to parse
			                         QiParse        &&prm_qi_parse   ///< The boost::spirit parser
			                         ) {
				T value;
				const bool ok = boost::spirit::qi::parse(
					prm_begin_itr,
					prm_end_itr,
					std::forward<QiParse>( prm_qi_parse ),
					value
				);

				if ( ! ok || prm_begin_itr != prm_end_itr ) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception(
						"Failed to parse a number (of type "s
						+ ::boost::core::demangle( typeid( T ).name() )
						+ R"() from ")"
						+ std::string{ prm_begin_itr, prm_end_itr }
						+ R"(")"
					));
				}
				return value;
			}

		} // namespace detail

		/// \brief Return an iterator pointing to the first point before a non-space character in the region between the specified
		///        iterators (or prm_end if none is found)
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		inline str_citr find_itr_before_first_non_space(const str_citr &prm_begin, ///< A  begin              iterator of the region of string to search
		                                                const str_citr &prm_end    ///< An end (one-past-end) iterator of the region of string to search
		                                                ) {
			return std::find_if(
				prm_begin,
				prm_end,
				[] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); }
			);
		}

		/// \brief Return an iterator pointing to the first point before a space character in the region between the specified
		///        iterators (or prm_end if none is found)
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		inline str_citr find_itr_before_first_space(const str_citr &prm_begin, ///< A  begin              iterator of the region of string to search
		                                            const str_citr &prm_end    ///< An end (one-past-end) iterator of the region of string to search
		                                            ) {
			return std::find_if(
				prm_begin,
				prm_end,
				[] (const auto &x) { return ( ( x == ' ' ) || ( x == '\t' ) ); }
			);
		}


		/// \brief Find the iterator that points to (just before) the first non-whitespace character
		///        in the specified string_ref
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		///
		/// \todo Come C++17, switch everything to use the ::std::string_view overload and then drop this one
		inline str_ref_citr find_itr_before_first_non_space(const boost::string_ref &prm_substring ///< The string_ref in which to search
		                                                    ) {
			const auto itr = boost::range::find_if(
				prm_substring,
				[] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); }
			);
			if ( itr == common::cend( prm_substring ) ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find any non-space chars in string"));
			}
			return itr;
		}

		/// \brief Find the iterator that points to (just before) the first non-whitespace character
		///        in the specified string_view
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		constexpr str_view_citr find_itr_before_first_non_space(const ::std::string_view &prm_substring ///< The string_view in which to search
		                                                        ) {
			for (auto itr = ::std::cbegin( prm_substring ); itr != ::std::cend( prm_substring ); ++itr) {
				if ( ( *itr != ' ' ) && ( *itr != '\t' ) ) {
					return itr;
				}
			}
			/// TODO: Come GCC >= 10, remove this silly dance to appeas it about the unconditionally non-constexpr throwing after here
			if ( !prm_substring.empty() && prm_substring.front() != ' ' && prm_substring.front() != '\t' ) {
				return ::std::cbegin( prm_substring );
			} else {
				BOOST_THROW_EXCEPTION( invalid_argument_exception( "Unable to find any non-space chars in string" ) );
			}
		}

		/// \brief Find the iterator that points to (just after) the last non-whitespace character
		///        in the specified string_ref
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		///
		/// \todo Come C++17, switch everything to use the ::std::string_view overload and then drop this one
		inline str_ref_citr find_itr_after_last_non_space(const boost::string_ref &prm_substring ///< The string_ref in which to search
		                                                  ) {
			const auto ritr = boost::range::find_if(
				prm_substring | boost::adaptors::reversed,
				[] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); }
			);
			if ( ritr.base() == common::cbegin( prm_substring ) ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find any non-space chars in string"));
			}
			return ritr.base();
		}

		/// \brief Find the iterator that points to (just after) the last non-whitespace character
		///        in the specified string_view
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		constexpr str_view_citr find_itr_after_last_non_space(const ::std::string_view &prm_substring ///< The string_view in which to search
		                                                      ) {
			for (auto itr = ::std::crbegin( prm_substring ); itr != ::std::crend( prm_substring ); ++itr) {
				if ( ( *itr != ' ' ) && ( *itr != '\t' ) ) {
					return itr.base();
				}
			}
			/// TODO: Come GCC >= 10, remove this silly dance to appeas it about the unconditionally non-constexpr throwing after here
			if ( !prm_substring.empty() && prm_substring.back() != ' ' && prm_substring.back() != '\t' ) {
				return ::std::crbegin( prm_substring ).base();
			} else {
				BOOST_THROW_EXCEPTION( invalid_argument_exception( "Unable to find any non-space chars in string" ) );
			}
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
		inline boost::string_ref dumb_trim_string_ref(const boost::string_ref &prm_substring ///< The string_ref to trim
		                                              ) {
			return make_string_ref(
				find_itr_before_first_non_space( prm_substring ),
				find_itr_after_last_non_space  ( prm_substring )
			);
		}

		/// \brief Return a pair of offsets to the section of the specified string_view after trimming
		///
		/// This can sensibly be called with a string argument.
		///
		/// This is dumb about whitespace (explicitly compares to ' ' and '\t'; ignores locale) for the sake of speed
		///
		/// \todo Fix: this given a string of all whitespace chars, this would currently cross over the begin/end pointers
		///       so should maybe calculated end pointer first and then restrict calculation for start pointer within
		///       range up to that end.
		template <typename T>
		constexpr ::std::pair<T, T> dumb_trim_string_view_to_offsets(const ::std::string_view &prm_substring ///< The string_view to trim
		                                                             ) {
			return {
				static_cast<T>( ::std::distance( ::std::cbegin( prm_substring ), find_itr_before_first_non_space( prm_substring ) ) ),
				static_cast<T>( ::std::distance( ::std::cbegin( prm_substring ), find_itr_after_last_non_space  ( prm_substring ) ) )
			};
		}

		/// \brief Parse a double from the field between the two specified string iterators
		inline double parse_double_from_field(const str_citr &prm_begin_itr, ///< A const_iterator pointing to the begin              of the field to be parsed
		                                      const str_citr &prm_end_itr    ///< A const_iterator pointing to the end (one-past-end) of the field to be parsed
		                                      ) {
			return detail::do_spirit_parse<double>(
				prm_begin_itr,
				prm_end_itr,
				boost::spirit::double_
			);
		}

		/// \brief Parse a float from the field between the two specified string iterators
		inline float parse_float_from_field(const str_citr &prm_begin_itr, ///< A const_iterator pointing to the begin              of the field to be parsed
		                                    const str_citr &prm_end_itr    ///< A const_iterator pointing to the end (one-past-end) of the field to be parsed
		                                    ) {
			return detail::do_spirit_parse<float>(
				prm_begin_itr,
				prm_end_itr,
				boost::spirit::float_
			);
		}

		/// \brief Parse an unsigned int from the field between the two specified string iterators
		inline unsigned int parse_uint_from_field(const str_citr &prm_begin_itr, ///< A const_iterator pointing to the begin              of the field to be parsed
		                                          const str_citr &prm_end_itr    ///< A const_iterator pointing to the end (one-past-end) of the field to be parsed
		                                          ) {
			return detail::do_spirit_parse<unsigned int>(
				prm_begin_itr,
				prm_end_itr,
				boost::spirit::uint_
			);
		}

		/// \brief Parse a (possibly space-padded) float from the specified region of string
		///
		/// Note: please benchmark any changes to these functions to ensure they stay fast
		inline float parse_float_from_substring(const std::string &prm_string, ///< The string containing the region to parse
		                                        const size_t      &prm_start,  ///< The index of the start of the region to parse
		                                        const size_t      &prm_length  ///< The length of the region to parse
		                                        ) {
			const auto begin_itr = std::next( prm_string.begin(), static_cast<ptrdiff_t>( prm_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( prm_length ) );
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
		inline double parse_double_from_substring(const std::string &prm_string, ///< The string containing the region to parse
		                                          const size_t      &prm_start,  ///< The index of the start of the region to parse
		                                          const size_t      &prm_length  ///< The length of the region to parse
		                                          ) {
			const auto begin_itr = std::next( prm_string.begin(), static_cast<ptrdiff_t>( prm_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( prm_length ) );
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
		inline int parse_int_from_substring(const std::string &prm_string, ///< The string containing the region to parse
		                                    const size_t      &prm_start,  ///< The index of the start of the region to parse
		                                    const size_t      &prm_length  ///< The length of the region to parse
		                                    ) {
			const auto begin_itr = std::next( prm_string.begin(), static_cast<ptrdiff_t>( prm_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( prm_length ) );
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
		inline unsigned int parse_uint_from_substring(const std::string &prm_string, ///< The string containing the region to parse
		                                              const size_t      &prm_start,  ///< The index of the start of the region to parse
		                                              const size_t      &prm_length  ///< The length of the region to parse
		                                              ) {
			const auto begin_itr = std::next( prm_string.begin(), static_cast<ptrdiff_t>( prm_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( prm_length ) );
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
		inline unsigned long int parse_ulong_from_substring(const std::string &prm_string, ///< The string containing the region to parse
		                                                    const size_t      &prm_start,  ///< The index of the start of the region to parse
		                                                    const size_t      &prm_length  ///< The length of the region to parse
		                                                    ) {
			const auto begin_itr = std::next( prm_string.begin(), static_cast<ptrdiff_t>( prm_start  ) );
			const auto end_itr   = std::next( begin_itr,          static_cast<ptrdiff_t>( prm_length ) );
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
		inline str_citr_str_citr_pair find_field_itrs(const std::string &prm_string,      ///< The string to search
		                                              const size_t      &prm_field_index, ///< The index of the field to find
		                                              const size_t      &prm_init_index,  ///< The index of the field from which the search should start
		                                              const str_citr    &prm_init_itr     ///< The iterator from which the search should start
		                                              ) {
			const auto end_itr = common::cend( prm_string );
			auto field_itr = find_itr_before_first_non_space( prm_init_itr, end_itr );
			for (const size_t field_ctr : boost::irange( prm_init_index, prm_field_index ) ) {
				boost::ignore_unused( field_ctr );
				field_itr = find_itr_before_first_space    ( field_itr, end_itr );
				field_itr = find_itr_before_first_non_space( field_itr, end_itr );
				if ( field_itr == end_itr ) {
					BOOST_THROW_EXCEPTION(runtime_error_exception(
						"Unable to find field "
						+ std::to_string( prm_field_index )
						+ " in line \""
						+ ( prm_string.size() > 103 ? ( prm_string.substr( 0, 100 ) + "[...]" ) : prm_string )
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
		inline str_citr_str_citr_pair find_field_itrs(const std::string &prm_string,      ///< The string to search
		                                              const size_t      &prm_field_index  ///< The index of the field to find
		                                              ) {
			return find_field_itrs(
				prm_string,
				prm_field_index,
				0,
				common::cbegin( prm_string )
			);
		}

		/// \brief Return an array<char, N> populated with N of the chars of the specified range of chars,
		///        (filling with 0s if the string isn't long enough)
		template <size_t N, typename Itr>
		inline std::array<char, N> get_char_arr_of_char_range(const Itr &prm_begin_itr, ///< The start of the range of characters to copy
		                                                      const Itr &prm_end_itr    ///< The end of the range of characters to copy
		                                                      ) {
			std::array<char, N> result;
			const auto end_at_n_itr = std::next( prm_begin_itr, N );
			const auto result_copied_itr = std::copy(
				prm_begin_itr,
				std::min( prm_end_itr, end_at_n_itr ),
				std::begin( result )
			);
			std::fill( result_copied_itr, std::end( result ), 0 );
			return result;
		}

		/// \brief Return an array<char, N> populated with N of the chars of the specified string starting from the
		///        specified index, (filling with 0s if the string isn't long enough)
		template <size_t N>
		inline std::array<char, N> get_char_arr_of_substring(const std::string &prm_string,     ///< The string from which to copy the characters
		                                                     const size_t      &prm_begin_index ///< The index at which to start copying characters from the string
		                                                     ) {
			const auto start_itr = std::next( common::cbegin( prm_string ), static_cast<ptrdiff_t>( prm_begin_index ) );
			const auto end_itr   = common::cend( prm_string );
			return get_char_arr_of_char_range<N>(
				std::min( start_itr, end_itr ),
				end_itr
			);
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_STRING_PARSE_TOOLS_HPP
