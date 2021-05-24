/// \file
/// \brief The data_file class header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DATA_FILE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DATA_FILE_HPP

#include <array>
#include <string>

#include "cath/common/algorithm/constexpr_is_uniq.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/lookup_arr_of_range_and_transform.hpp"

namespace cath::file {

	/// \brief Represent the different file types that a data_dirs_options_block manages
	enum class data_file : unsigned int {
		PDB,  ///< Enum value to denote PDB  files
		DSSP, ///< Enum value to denote DSSP files
		WOLF, ///< Enum value to denote wolf files
		SEC   ///< Enum value to denote sec  files
	};

	namespace detail {

		/// \brief All the data_file values
		static constexpr auto all_data_file_types =
		  common::make_array( data_file::PDB, data_file::DSSP, data_file::WOLF, data_file::SEC );

		static_assert( common::constexpr_is_uniq( all_data_file_types ),
		               "all_data_file_types shouldn't contain repeated values" );

	} // namespace detail

	/// \brief The names for each of the data file types
	///
	/// These names are used for:
	///  - referring to the file types in the options descriptions
	///  - constructing the options names (after being lower-cased and having spaces and underscores replaced with hyphens)
	constexpr ::std::string_view name_of_data_file( const data_file &prm_data_file ) {
		// clang-format off
		switch ( prm_data_file ) {
			case ( data_file::PDB  ) : { return ::std::string_view{ "PDB"  }; }
			case ( data_file::DSSP ) : { return ::std::string_view{ "DSSP" }; }
			case ( data_file::WOLF ) : { return ::std::string_view{ "wolf" }; }
			case ( data_file::SEC  ) : { return ::std::string_view{ "sec"  }; }
		}
		// clang-format on
		BOOST_THROW_EXCEPTION( common::out_of_range_exception( "Unhandled data_option" ) );
	}

	/// \brief The default sub-directory name to append to a cath-root-dir to get the file type's directory
	constexpr ::std::string_view default_subdir_name_of_data_file( const data_file &prm_data_file ) {
		// clang-format off
		switch ( prm_data_file ) {
			case ( data_file::PDB  ) : { return ::std::string_view{ "pdb"  }; }
			case ( data_file::DSSP ) : { return ::std::string_view{ "dssp" }; }
			case ( data_file::WOLF ) : { return ::std::string_view{ "wolf" }; }
			case ( data_file::SEC  ) : { return ::std::string_view{ "sec"  }; }
		}
		// clang-format on
		BOOST_THROW_EXCEPTION( common::out_of_range_exception( "Unhandled data_option" ) );
	}

	inline constexpr auto DATA_FILE_NAMES =
	  common::lookup_arr_of_range_and_transform( detail::all_data_file_types, name_of_data_file );

	::std::string to_string( const data_file & );

	::std::ostream &operator<<( ::std::ostream &, const data_file & );

	size_t str_length_of_data_file( const data_file & );

	size_t max_data_file_str_length();

} // namespace cath::file

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DATA_FILE_HPP
