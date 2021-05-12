/// \file
/// \brief The temp_file class definitions

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

#include "temp_file.hpp"

#include <filesystem>
#include <random>
#include <string>
#include <string_view>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/common/boost_addenda/range/to_vector.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::cath;
using namespace ::cath::common;

using ::boost::adaptors::transformed;
using ::boost::none;
using ::std::filesystem::path;
using ::std::filesystem::temp_directory_path;
using ::std::mt19937;
using ::std::random_device;
using ::std::string;
using ::std::string_view;
using ::std::uniform_int_distribution;

/// \brief Constructor for temp_file.
temp_file::temp_file(const string &prm_filename_pattern ///< A pattern for the filename to create with % symbols for characters to be altered (eg ".%%%%-%%%%-%%%%-%%%%.pml")
                     ) : filename( ! prm_filename_pattern.empty() ? path_opt( temp_filename_of_basename_pattern( prm_filename_pattern ) ) : none ) {
}

/// \brief A function to construct a temporary filename from a pattern.
///
/// \param prm_filename_pattern A pattern for the filename to create with % symbols for characters to be altered (eg ".%%%%-%%%%-%%%%-%%%%.pml")
path temp_file::temp_filename_of_basename_pattern( const string &prm_filename_pattern ) {
	if ( path( prm_filename_pattern ).has_parent_path() ) {
		BOOST_THROW_EXCEPTION( invalid_argument_exception( "temp_file pattern \"" + prm_filename_pattern
		                                                   + "\" should just be a basename, not a full path" ) );
	}

	constexpr string_view            PALETTE = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
	mt19937                          random_gen( random_device{}() );
	uniform_int_distribution<size_t> dist( 0, PALETTE.length() - 1 );

	// clang-format off
	const string populated_pattern = prm_filename_pattern
		| transformed( [&] ( const char &x ) {
			return ( x == '%' ) ? PALETTE[ dist( random_gen ) ] : x;
		} )
		| to_string;
	// clang-format on

	return temp_directory_path() / populated_pattern;
}

/// \brief Destructor for temp_file
temp_file::~temp_file() noexcept {
	try {
		if ( filename ) {
			if ( exists( *filename ) ) {
				remove( *filename );
			}
		}
	}
	catch (...) {
	}
}
/// \brief TODOCUMENT
const path_opt & temp_file::get_opt_filename() const {
	return filename;
}

/// \brief TODOCUMENT
bool cath::common::has_filename(const temp_file &prm_temp_file ///< TODOCUMENT
                                ) {
	return static_cast<bool>( prm_temp_file.get_opt_filename() );
}
/// \brief TODOCUMENT
path cath::common::get_filename(const temp_file &prm_temp_file ///< TODOCUMENT
                                ) {
	return *prm_temp_file.get_opt_filename();
}
