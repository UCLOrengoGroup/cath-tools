/// \file
/// \brief The data_dirs_spec class definitions

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

#include "data_dirs_spec.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/algorithm/contains.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/file/find_file.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/runtime_error_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::file::detail;
using namespace cath::opts;
using namespace cath::opts::detail;
using namespace std;

using boost::adaptors::transformed;
using boost::algorithm::is_any_of;
using boost::algorithm::join;
using boost::algorithm::token_compress_on;
using boost::filesystem::path;
using boost::lexical_cast;

/// \brief Default values of each of the options (path, prefix, suffix) for each of the file types
const data_dirs_spec::file_option_str_map_map data_dirs_spec::DATA_FILE_TYPE_OPTION_DEFAULTS = {
	{ data_file::PDB,  data_option_str_map{ { data_option::PATH,   "."     },
	                                        { data_option::PREFIX, ""      },
	                                        { data_option::SUFFIX, ""      } } },
	{ data_file::DSSP, data_option_str_map{ { data_option::PATH,   "."     },
	                                        { data_option::PREFIX, ""      },
	                                        { data_option::SUFFIX, ".dssp" } } },
	{ data_file::WOLF, data_option_str_map{ { data_option::PATH,   "."     },
	                                        { data_option::PREFIX, ""      },
	                                        { data_option::SUFFIX, ".wolf" } } },
	{ data_file::SEC,  data_option_str_map{ { data_option::PATH,   "."     },
	                                        { data_option::PREFIX, ""      },
	                                        { data_option::SUFFIX, ".sec"  } } }
};

/// \brief The names for each of the data file types
///
/// These names are used for:
///  - referring to the file types in the options descriptions
///  - constructing the options names (after being lower-cased and having spaces and underscores replaced with hyphens)
const data_file_str_map data_dirs_spec::DATA_FILE_NAMES = {
	{ data_file::PDB,  "PDB"  },
	{ data_file::DSSP, "DSSP" },
	{ data_file::WOLF, "wolf" },
	{ data_file::SEC,  "sec"  }
};

/// \brief The default sub-directory name to append to a cath-root-dir to get the file type's directory
const data_file_str_map data_dirs_spec::DEFAULT_SUBDIR_NAME = {
	{ data_file::PDB,  "pdb"  },
	{ data_file::DSSP, "dssp" },
	{ data_file::WOLF, "wolf" },
	{ data_file::SEC,  "sec"  }
};

/// \brief Ctor for data_dirs_spec
data_dirs_spec::data_dirs_spec() : values( DATA_FILE_TYPE_OPTION_DEFAULTS ) {
}

/// \brief Get the value of an option type and data type
///
/// data_option is a private, implementation detail so this method provides a common interface for
/// the public getters to use to access values.
string data_dirs_spec::get_value_of_option_and_data_file(const data_option &arg_option_type, ///< The option type to query
                                                         const data_file   &arg_data_file    ///< The data type to query
                                                         ) const {
	if ( ! contains( values, arg_data_file ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to any value for data file type "
			+ lexical_cast<string>(arg_data_file)
			+ " in data_dirs_options_block"
		));
	}

	const data_option_str_map &values_of_arg_data_file = values.find( arg_data_file )->second;
	if ( ! contains( values_of_arg_data_file, arg_option_type ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to find value for option type "
			+ lexical_cast<string>(arg_option_type)
			+ " and data file type "
			+ lexical_cast<string>(arg_data_file)
			+ " in data_dirs_options_block"
		));
	}
	return values_of_arg_data_file.find( arg_option_type )->second;
}

/// \brief Getter for the name of a given data type
///
/// This is static because the name of each data type is static.
string data_dirs_spec::get_name_of_data_file(const data_file &arg_data_file ///< The data type to query
                                             ) {
	return DATA_FILE_NAMES.at( arg_data_file );
}

/// \brief Getter for cath_root_dir
const path & data_dirs_spec::get_cath_root_dir() const {
	return cath_root_dir;
}

/// \brief TODOCUMENT
data_dirs_spec & data_dirs_spec::set_value_of_option_and_data_file(const data_option &arg_option_type, ///< TODOCUMENT
                                                                   const data_file   &arg_data_file,   ///< TODOCUMENT
                                                                   const string      &arg_value        ///< TODOCUMENT
                                                                   ) {
	values[ arg_data_file ][ arg_option_type ] = arg_value;
	return *this;
}

/// \brief TODOCUMENT
data_dirs_spec & data_dirs_spec::set_path_of_data_file(const data_file &arg_data_file, ///< TODOCUMENT
                                                       const string    &arg_path       ///< TODOCUMENT
                                                       ) {
	set_value_of_option_and_data_file(
		data_option::PATH,
		arg_data_file,
		arg_path
	);
	return *this;
}

/// \brief Setter for cath_root_dir
data_dirs_spec & data_dirs_spec::set_cath_root_dir(const path &arg_cath_root_dir ///< The cath_root_dir value to set
                                                   ) {
	cath_root_dir = arg_cath_root_dir;
	return *this;
}

/// \brief Getter for the path value of a given data type
string cath::opts::get_path_of_data_file(const data_dirs_spec &arg_data_dirs_spec, ///< TODOCUMENT
                                         const data_file      &arg_data_file       ///< The data type to query
                                         ) {
	return arg_data_dirs_spec.get_value_of_option_and_data_file( data_option::PATH, arg_data_file );
}

/// \brief Getter for the prefix value of a given data type
string cath::opts::get_prefix_of_data_file(const data_dirs_spec &arg_data_dirs_spec, ///< TODOCUMENT
                                           const data_file      &arg_data_file ///< The data type to query
                                           ) {
	return arg_data_dirs_spec.get_value_of_option_and_data_file( data_option::PREFIX, arg_data_file );
}

/// \brief Getter for the suffix value of a given data type
string cath::opts::get_suffix_of_data_file(const data_dirs_spec &arg_data_dirs_spec, ///< TODOCUMENT
                                           const data_file      &arg_data_file ///< The data type to query
                                           ) {
	return arg_data_dirs_spec.get_value_of_option_and_data_file( data_option::SUFFIX, arg_data_file );
}

/// \brief TODOCUMENT
///
/// \relates data_dirs_spec
data_dirs_spec cath::opts::build_data_dirs_spec_of_path(const path_vec &arg_path ///< TODOCUMENT
                                                        ) {
	const string arg_path_string = join_directories_into_path(arg_path);
	data_dirs_spec new_data_dirs;
	for (const data_file &the_data_file : all_data_file_types) {
		new_data_dirs.set_path_of_data_file( the_data_file, arg_path_string );
	}
	return new_data_dirs;
}

/// \brief TODOCUMENT
///
/// \relates data_dirs_spec
data_dirs_spec cath::opts::build_data_dirs_spec_of_dir(const path &arg_path ///< TODOCUMENT
                                                       ) {
	return build_data_dirs_spec_of_path( { arg_path } );
}

/// \brief Convenience function to return a path_vec of directories for the path for a given data_file
///
/// \relates data_dirs_spec
path_vec cath::opts::get_paths_of_data_file(const data_dirs_spec &arg_data_dirs_options, ///< The data_dirs_spec from which to grab the path
                                            const data_file      &arg_data_file          ///< The data type to query
                                            ) {
	const path   cath_root_dir = arg_data_dirs_options.get_cath_root_dir();
	const string path_string   = get_path_of_data_file( arg_data_dirs_options, arg_data_file );
	path_vec     path_dirs     = split_path_into_directories( path_string );
	if ( ! cath_root_dir.empty() ) {
		path_dirs.push_back( cath_root_dir / data_dirs_spec::DEFAULT_SUBDIR_NAME.at( arg_data_file ) );
	}
	return path_dirs;
}

/// \brief Attempt to find a file of a given type for a given name using the preferences specified in a given data_dirs_spec
///
/// \returns The found file or an empty path object if one could not be found
///
/// \relates data_dirs_spec
path cath::opts::find_file(const data_dirs_spec &arg_data_dirs_options, ///< The data_dirs_spec specifying how to find the file
                           const data_file      &arg_data_file,         ///< The data type of file that's required
                           const string         &arg_name               ///< The name for which to find the file (eg 1c0pA01)
                           ) {
	const path_vec path_dirs       = get_paths_of_data_file( arg_data_dirs_options,  arg_data_file );

	const string   file_type_name  = arg_data_dirs_options.get_name_of_data_file  ( arg_data_file );
	const string   prefix          = get_prefix_of_data_file( arg_data_dirs_options, arg_data_file );
	const string   suffix          = get_suffix_of_data_file( arg_data_dirs_options, arg_data_file );

	const string   basename        = prefix + arg_name + suffix;

	const path     found_file      = common::find_file( path_dirs, basename );

	// If nothing was found, then throw an informative exception
	if ( found_file.empty() ) {
		ostringstream error_ss;
		error_ss << "Unable to find a "                 << file_type_name
		         << " file called "                     << basename
		         << " (built from the prefix \""        << prefix
		         << "\", the name \""                   << arg_name
		         << "\" and the suffix \""              << suffix
		         << "\") in any of directories in the " << file_type_name
		         << " path [";
		for (const path &path_dir : path_dirs) {
			error_ss << " " << path_dir.string();
		}
		error_ss << " ]";
		BOOST_THROW_EXCEPTION(runtime_error_exception(error_ss.str()));
	}

	// Otherwise, return the file that was found
	return found_file;
}

/// \brief Convenience function to split a string containing a path of colon-separated directories
///        into a vector of path objects
///
/// Empty paths are not added to the returned vector of paths
/// (eg input "/tmp::/var" would only return a vector with two entries)
///
/// \relates data_dirs_spec
path_vec cath::opts::split_path_into_directories(const string &arg_path_string ///< The path of colon-separated directory names
                                                 ) {
	str_vec dir_strings = split_build<str_vec>( arg_path_string, is_any_of( ":" ), token_compress_on );
	path_vec dirs;
	dirs.reserve(dir_strings.size());
	for (const string &dir_string : dir_strings) {
		if ( ! dir_string.empty() ) {
			dirs.push_back( path( dir_string ) );
		}
	}
	return dirs;
}

/// \brief Convenience function to join a vector of directories into a colon-separated string
///
/// \relates data_dirs_spec
string cath::opts::join_directories_into_path(const path_vec &arg_dirs ///< The vector of directories
                                              ) {
	// Return the result of joining the strings of the arg_dirs with colon characters
	return join(
		arg_dirs | transformed( [] (const path &x) { return x.string(); } ),
		":"
	);
}

