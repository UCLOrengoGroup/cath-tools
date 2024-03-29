/// \file
/// \brief The data_dirs_spec class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_DIRS_SPEC_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_DIRS_SPEC_HPP

#include <array>
#include <filesystem>
#include <map>
#include <string_view>
#include <tuple>

#include "cath/common/path_type_aliases.hpp"
#include "cath/file/data_file.hpp"
#include "cath/file/options/data_option.hpp"

// clang-format off
namespace cath::file { class name_set; }
// clang-format on

namespace cath::opts {

	using data_file_str_map    = std::map<file::data_file, std::string>;
	// using data_file_str_pair   = data_file_str_map::value_type;
	using data_option_str_map  = std::map<detail::data_option, std::string>;
	// using data_option_str_pair = data_option_str_map::value_type;
	using file_option_str_map_map = std::map<file::data_file, data_option_str_map>;

	/// \brief Store the program options required to find various different types of file from a name
	///        (eg find a DSSP file from the name 1c0pA01)
	///
	/// Note that there are two different uses of the word 'path' around this code:
	///  * the Boost filesystem type that is used to store the location of a file or directory
	///  * a prioritised list of directories to search for some sort of file
	///
	/// For each file type, options are provided to specify
	///  - the prefix, which is prepended to the name when constructing the filename,
	///  - the suffix, which is appended to the name when constructing the filename and
	///  - the path, which is a colon-separated list of directories through which to search for the file (in descending order of preference).
	///
	/// \todo Add a (non-member, non-friend) function to do the reverse operation: ie use a data_dirs_spec and data_file to deduce a
	///       sensible name from a file (by stripping of the path and then any suitable leading prefix or trailing suffix that
	///       don't leave an empty name). For example, this would generate the name 1c0pA01 from a dssp filename "/cath/data/current/dssp/1c0pA01.dssp").
	class data_dirs_spec final {
	private:

		/// \brief Store the actual values by data_file and then data_option
		file_option_str_map_map values;

		/// \brief A root directory from which the other directories can be constructed
		::std::filesystem::path cath_root_dir;

	public:
		data_dirs_spec();

		[[nodiscard]] std::string get_value_of_option_and_data_file( const detail::data_option &,
		                                                             const file::data_file & ) const;

		static std::string_view get_name_of_data_file(const file::data_file &);

		[[nodiscard]] const ::std::filesystem::path &get_cath_root_dir() const;

		data_dirs_spec & set_value_of_option_and_data_file(const detail::data_option &,
		                                                   const file::data_file &,
		                                                   const std::string &);

		data_dirs_spec & set_path_of_data_file(const file::data_file &,
		                                       const std::string &);

		data_dirs_spec & set_cath_root_dir(const ::std::filesystem::path &);

		/// \brief Default values of each of the options (path, prefix, suffix) for each of the file types
		static constexpr ::std::array DATA_FILE_TYPE_OPTION_DEFAULTS = {
			// clang-format off
			::std::tuple{ file::data_file::PDB,  detail::data_option::PATH,   ::std::string_view{ "."     } },
			::std::tuple{ file::data_file::PDB,  detail::data_option::PREFIX, ::std::string_view{ ""      } },
			::std::tuple{ file::data_file::PDB,  detail::data_option::SUFFIX, ::std::string_view{ ""      } },
			::std::tuple{ file::data_file::DSSP, detail::data_option::PATH,   ::std::string_view{ "."     } },
			::std::tuple{ file::data_file::DSSP, detail::data_option::PREFIX, ::std::string_view{ ""      } },
			::std::tuple{ file::data_file::DSSP, detail::data_option::SUFFIX, ::std::string_view{ ".dssp" } },
			::std::tuple{ file::data_file::WOLF, detail::data_option::PATH,   ::std::string_view{ "."     } },
			::std::tuple{ file::data_file::WOLF, detail::data_option::PREFIX, ::std::string_view{ ""      } },
			::std::tuple{ file::data_file::WOLF, detail::data_option::SUFFIX, ::std::string_view{ ".wolf" } },
			::std::tuple{ file::data_file::SEC,  detail::data_option::PATH,   ::std::string_view{ "."     } },
			::std::tuple{ file::data_file::SEC,  detail::data_option::PREFIX, ::std::string_view{ ""      } },
			::std::tuple{ file::data_file::SEC,  detail::data_option::SUFFIX, ::std::string_view{ ".sec"  } },
			// clang-format on
		};
	};

	file_option_str_map_map data_file_type_option_defaults_map();

	std::string get_path_of_data_file(const data_dirs_spec &,
	                                  const file::data_file &);
	std::string get_prefix_of_data_file(const data_dirs_spec &,
	                                    const file::data_file &);
	std::string get_suffix_of_data_file(const data_dirs_spec &,
	                                    const file::data_file &);

	data_dirs_spec build_data_dirs_spec_of_path(const path_vec &);

	data_dirs_spec build_data_dirs_spec_of_dir(const ::std::filesystem::path &);

	path_vec get_paths_of_data_file(const data_dirs_spec &,
	                                const file::data_file &);

	::std::filesystem::path find_file(const data_dirs_spec &,
	                                  const file::data_file &,
	                                  const std::string &);

	::std::filesystem::path find_file(const data_dirs_spec &,
	                                  const file::data_file &,
	                                  const file::name_set &);

	path_vec split_path_into_directories(const std::string &);

	std::string join_directories_into_path(const path_vec &);

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_DIRS_SPEC_HPP
