/// \file
/// \brief The string_of_rapidjson_write header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_STRING_OF_RAPIDJSON_WRITE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_STRING_OF_RAPIDJSON_WRITE_HPP

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/find.hpp>

#include "cath/common/algorithm/for_n.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/cpp17/invoke.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/rapidjson_addenda/rapidjson_writer.hpp"

namespace cath {
	namespace common {

		/// \brief Create a string from performing the specified write on a rapidjson_writer of the specified style
		///        and nest within the specified number of levels of depth
		///
		/// The extra levels of depth affect the indentation between lines of the resulting string
		template <json_style Style, typename Fn>
		std::string string_of_rapidjson_write(Fn           &&prm_fn,             ///< The write to perform on the writer (must be invokable with a single `rapidjson_writer<Style>&` argument)
		                                      const size_t  &prm_extra_depth = 0 ///< The number of levels of depth
		                                      ) {
			// Construct a writer
			rapidjson_writer<Style> the_writer;

			// Perform the specified function on the writer, within prm_extra_depth levels of array
			for_n( prm_extra_depth, [&] { the_writer.start_array(); } );
			invoke( std::forward<Fn>( prm_fn ), the_writer );
			for_n( prm_extra_depth, [&] { the_writer.end_array();   } );

			// Grab the C-style string
			const char * const c_string   = the_writer.get_c_string();

			// Find the end of the C-style string and make an iterator_range over the C-style string
			const auto c_string_end = std::next( c_string, static_cast<ptrdiff_t>( strlen( c_string ) ) );
			const auto string_range = boost::make_iterator_range( c_string, c_string_end );

			// Find the end of the any preceding prm_extra_depth array-open chars
			// and the start of any following prm_extra_depth array-close chars
			const auto end_of_preceding   = ( prm_extra_depth > 0 )
			                                ? common::cend  ( boost::algorithm::find_nth( string_range, "[",     debug_numeric_cast<int>( prm_extra_depth - 1 ) ) )
			                                : c_string;
			const auto begin_of_following = ( prm_extra_depth > 0 )
			                                ? common::cbegin( boost::algorithm::find_nth( string_range, "]", 0 - debug_numeric_cast<int>( prm_extra_depth     ) ) )
			                                : c_string_end;

			// Get an iterator_range of the trimmed section between those two points
			const auto trimmed_range = boost::algorithm::trim_copy( boost::make_iterator_range( end_of_preceding, begin_of_following ) );

			// Return a std::string of that range
			return {
				common::cbegin( trimmed_range ),
				common::cend  ( trimmed_range )
			};
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_STRING_OF_RAPIDJSON_WRITE_HPP
