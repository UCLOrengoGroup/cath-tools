/// \file
/// \brief The data_dirs_options_block class definitions

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

#include "data_dirs_options_block.h"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include "common/algorithm/contains.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/c++14/make_unique.h"
#include "common/clone/make_uptr_clone.h"
#include "exception/invalid_argument_exception.h"
#include "exception/runtime_error_exception.h"
#include "file/data_file.h"

using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::file::detail;
using namespace cath::opts::detail;
using namespace cath::opts;
using namespace std;

using boost::adaptors::transformed;
using boost::algorithm::is_any_of;
using boost::algorithm::join;
using boost::algorithm::replace_all_copy;
using boost::algorithm::token_compress_on;
using boost::lexical_cast;
using boost::none;
using boost::range::max_element;

/// \brief Default values of each of the options (path, prefix, suffix) for each of the file types
const data_dirs_options_block::file_option_str_map_map data_dirs_options_block::DATA_FILE_TYPE_OPTION_DEFAULTS = {
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
const data_dirs_options_block::data_file_str_map data_dirs_options_block::DATA_FILE_NAMES = {
	{ data_file::PDB,  "PDB"  },
	{ data_file::DSSP, "DSSP" },
	{ data_file::WOLF, "wolf" },
	{ data_file::SEC,  "sec"  }
};

/// \brief For each of the options (path, prefix, suffix), the suffixes that are appended to the names of the file types to produce the option name
const data_dirs_options_block::data_option_str_map data_dirs_options_block::DATA_OPTION_SUFFIXES = {
	{ data_option::PATH,   "-path"   },
	{ data_option::PREFIX, "-prefix" },
	{ data_option::SUFFIX, "-suffix" }
};

/// \brief For each of the options (path, prefix, suffix), the start of the description up to the point of the data type name
const data_dirs_options_block::data_option_str_map data_dirs_options_block::DATA_OPTION_DESCRIPTION_START = {
	{ data_option::PATH,   "Search for "                                             },
	{ data_option::PREFIX, "Prepend the prefix arg to a protein's name to form its " },
	{ data_option::SUFFIX, "Append the suffix arg to a protein's name to form its "  }
};

/// \brief For each of the options (path, prefix, suffix), the end of the description from the point of the data type name
const data_dirs_options_block::data_option_str_map data_dirs_options_block::DATA_OPTION_DESCRIPTION_END = {
	{ data_option::PATH,   " files using the path arg" },
	{ data_option::PREFIX, " filename"                 },
	{ data_option::SUFFIX, " filename"                 }
};

/// \brief The default sub-directory name to append to a cath-root-dir to get the file type's directory
const data_dirs_options_block::data_file_str_map data_dirs_options_block::DEFAULT_SUBDIR_NAME = {
	{ data_file::PDB,  "pdb"  },
	{ data_file::DSSP, "dssp" },
	{ data_file::WOLF, "wolf" },
	{ data_file::SEC,  "sec"  }
};

/// \brief The option name for the CATH root directory
const string data_dirs_options_block::PO_CATH_ROOT_DIR{ "cath-root-dir" };

/// \brief A standard do_clone method
unique_ptr<options_block> data_dirs_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string data_dirs_options_block::do_get_block_name() const {
	return "Conversion between a protein's name and its data files";
}

/// \brief Add this block's options to the provided options_description
///
/// This adds one entry for each of the file types for each of the option types
void data_dirs_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                    ) {
	// Loop over each of the data option types (path, prefix, suffix) and grab the data_option and suffix string
	for (const data_option_str_pair &data_option_suffix_pair : DATA_OPTION_SUFFIXES) {
		const data_option &the_data_option    = data_option_suffix_pair.first;
		const string      &data_option_suffix = data_option_suffix_pair.second;

		// Loop over each of the data file types and grab the data_file and name string
		for (const data_file_str_pair &data_file_name_pair : DATA_FILE_NAMES) {
			const data_file &the_data_file  = data_file_name_pair.first;
			const string    &data_file_name = data_file_name_pair.second;

			// Tidy up the file type name and append the option suffix to get the full option name (eg "PDB" -> "pdb", "pdb" + "-path" -> "pdb-path")
			const string     tidy_file_name = replace_all_copy(replace_all_copy(to_lower_copy(data_file_name), "_", "-"), " ", "-");
			const string     option_name    = tidy_file_name + data_option_suffix;

			// Get a reference to the value for this file/option pair
			// (which will be created by the operator[]() call if it doesn't already exist)
			string          &str_value     = values[the_data_file][the_data_option];

			// Grab the default, description start and description end and join them together to get a description for the option
			const string     default_value = DATA_FILE_TYPE_OPTION_DEFAULTS.at( the_data_file   ).at( the_data_option );
			const string     desc_start    = DATA_OPTION_DESCRIPTION_START.at(  the_data_option );
			const string     desc_end      = DATA_OPTION_DESCRIPTION_END.at(    the_data_option );
			const string     description   = desc_start + data_file_name + desc_end;

			// Add a new option, using the accumulated data
			arg_desc.add_options()
				(option_name.c_str(), value<string>(&str_value)->default_value(default_value), description.c_str());
		}
	}
}

/// \brief Add any hidden options to the provided options_description
void data_dirs_options_block::do_add_hidden_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                   ) {
	arg_desc.add_options()
		( PO_CATH_ROOT_DIR.c_str(), value<path>( &cath_root_dir ), "A root directory from which CATH Tools programs can find sub-directories of standard names (\"pdb\", \"dssp\", \"wolf\", \"sec\")" );
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// At present, this always accepts
opt_str data_dirs_options_block::do_invalid_string() const {
	if (!cath_root_dir.empty() && !is_acceptable_input_dir(cath_root_dir)) {
		return "CATH root directory \"" + cath_root_dir.string() + "\" is not a valid input directory";
	}

	return none;
}

/// \brief TODOCUMENT
void data_dirs_options_block::set_value_of_option_and_data_file(const data_option &arg_option_type, ///< TODOCUMENT
                                                                const data_file   &arg_data_file,   ///< TODOCUMENT
                                                                const string      &arg_value        ///< TODOCUMENT
                                                                ) {
	values[arg_data_file][arg_option_type] = arg_value;
}

/// \brief Get the value of an option type and data type
///
/// data_option is a private, implementation detail so this method provides a common interface for
/// the public getters to use to access values.
string data_dirs_options_block::get_value_of_option_and_data_file(const data_option &arg_option_type, ///< The option type to query
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

/// \brief Ctor for data_dirs_options_block
data_dirs_options_block::data_dirs_options_block() {
	// Loop over each of the data option types (path, prefix, suffix) and the data file types
	for (const data_option_str_pair &data_option_suffix_pair : DATA_OPTION_SUFFIXES) {
		const data_option &the_data_option = data_option_suffix_pair.first;
		for (const data_file_str_pair &data_file_name_pair : DATA_FILE_NAMES) {
			const data_file &the_data_file     = data_file_name_pair.first;
			const string          default_value = DATA_FILE_TYPE_OPTION_DEFAULTS.find( the_data_file )->second.find(the_data_option)->second;
			values[the_data_file][the_data_option]      = default_value;
		}
	}
}

/// \brief TODOCUMENT
void data_dirs_options_block::set_path_of_data_file(const data_file &arg_data_file, ///< TODOCUMENT
                                                    const string    &arg_path       ///< TODOCUMENT
                                                    ) {
	set_value_of_option_and_data_file(
		data_option::PATH,
		arg_data_file,
		arg_path
	);
}

/// \brief Getter for the name of a given data type
///
/// This is static because the name of each data type is static.
string data_dirs_options_block::get_name_of_data_file(const data_file &arg_data_file ///< The data type to query
                                                      ) {
	return DATA_FILE_NAMES.at( arg_data_file );
}

/// \brief Getter for cath_root_dir
path data_dirs_options_block::get_cath_root_dir() const {
	return cath_root_dir;
}

/// \brief Getter for the path value of a given data type
string data_dirs_options_block::get_path_of_data_file(const data_file &arg_data_file ///< The data type to query
                                                      ) const {
	return get_value_of_option_and_data_file(data_option::PATH, arg_data_file);
}

/// \brief Getter for the prefix value of a given data type
string data_dirs_options_block::get_prefix_of_data_file(const data_file &arg_data_file ///< The data type to query
                                                        ) const {
	return get_value_of_option_and_data_file(data_option::PREFIX, arg_data_file);
}

/// \brief Getter for the suffix value of a given data type
string data_dirs_options_block::get_suffix_of_data_file(const data_file &arg_data_file ///< The data type to query
                                                        ) const {
	return get_value_of_option_and_data_file(data_option::SUFFIX, arg_data_file);
}

/// \brief TODOCUMENT
///
/// \relates data_dirs_options_block
data_dirs_options_block cath::opts::build_data_dirs_options_block_of_path(const path_vec &arg_path ///< TODOCUMENT
                                                                          ) {
	const string arg_path_string = join_directories_into_path(arg_path);
	data_dirs_options_block new_data_dirs;
	for (const data_file &the_data_file : all_data_file_types) {
		new_data_dirs.set_path_of_data_file( the_data_file, arg_path_string );
	}
	return new_data_dirs;
}

/// \brief TODOCUMENT
///
/// \relates data_dirs_options_block
data_dirs_options_block cath::opts::build_data_dirs_options_block_of_dir(const path &arg_path ///< TODOCUMENT
                                                                         ) {
	return build_data_dirs_options_block_of_path( { arg_path } );
}

/// \brief Convenience function to return a path_vec of directories for the path for a given data_file
///
/// \relates data_dirs_options_block
path_vec cath::opts::get_path_of_data_file(const data_dirs_options_block &arg_data_dirs_options, ///< The data_dirs_options_block from which to grab the path
                                           const data_file               &arg_data_file          ///< The data type to query
                                           ) {
	const path   cath_root_dir = arg_data_dirs_options.get_cath_root_dir();
	const string path_string   = arg_data_dirs_options.get_path_of_data_file( arg_data_file );
	path_vec     path_dirs     = split_path_into_directories( path_string );
	if ( ! cath_root_dir.empty() ) {
		path_dirs.push_back( cath_root_dir / data_dirs_options_block::DEFAULT_SUBDIR_NAME.at( arg_data_file ) );
	}
	return path_dirs;
}

/// \brief Attempt to find a file of a given type for a given name using the preferences specified in a given data_dirs_options_block
///
/// \returns The found file or an empty path object if one could not be found
///
/// \relates data_dirs_options_block
path cath::opts::find_file(const data_dirs_options_block &arg_data_dirs_options, ///< The data_dirs_options_block specifying how to find the file
                           const data_file               &arg_data_file,         ///< The data type of file that's required
                           const string                  &arg_name               ///< The name for which to find the file (eg 1c0pA01)
                           ) {
	const path_vec path_dirs       = get_path_of_data_file( arg_data_dirs_options,  arg_data_file );

	const string       file_type_name  = arg_data_dirs_options.get_name_of_data_file  ( arg_data_file );
	const string       prefix          = arg_data_dirs_options.get_prefix_of_data_file( arg_data_file );
	const string       suffix          = arg_data_dirs_options.get_suffix_of_data_file( arg_data_file );

	const string basename = prefix + arg_name + suffix;

	const path found_file = find_file(path_dirs, basename);

	// If nothing was found, then throw an informative exception
	if (found_file.empty()) {
		ostringstream error_ss;
		error_ss << "Unable to find a "         << file_type_name <<
			" file called "                     << basename       <<
			" (built from the prefix \""        << prefix         <<
			"\", the name \""                   << arg_name       <<
			"\" and the suffix \""              << suffix         <<
			"\") in any of directories in the " << file_type_name <<
			" path [";
		for (const path &path_dir : path_dirs) {
			error_ss << " " << path_dir.string();
		}
		error_ss << " ]";
		BOOST_THROW_EXCEPTION(runtime_error_exception(error_ss.str()));
	}

	// Otherwise, return the file that was found
	return found_file;
}

/// \brief Search for a particular file basename through a path of directories
///
/// \returns The found file or an empty path object if one could not be found
///
/// \relates data_dirs_options_block
path cath::opts::find_file(const path_vec &arg_path_dirs, ///< Directories through which to search for the file (in descending order of preference)
                           const string       &arg_basename   ///< The basename of the file to be located (eg 1c0pA01.dssp)
                           ) {
	for (const path &dir : arg_path_dirs) {
		const path potential_file = dir / arg_basename;
		if (exists(potential_file)) {
			return potential_file;
		}
	}
	return path();
}

/// \brief Convenience function to split a string containing a path of colon-separated directories
///        into a vector of path objects
///
/// Empty paths are not added to the returned vector of paths
/// (eg input "/tmp::/var" would only return a vector with two entries)
///
/// \relates data_dirs_options_block
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
/// \relates data_dirs_options_block
string cath::opts::join_directories_into_path(const path_vec &arg_dirs ///< The vector of directories
                                              ) {
	// Return the result of joining the strings of the arg_dirs with colon characters
	return join(
		arg_dirs | transformed( [] (const path &x) { return x.string(); } ),
		":"
	);
}

