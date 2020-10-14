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

#include "data_dirs_options_block.hpp"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/file/options/data_dirs_spec.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts::detail;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::algorithm::replace_all_copy;
using ::boost::algorithm::to_lower_copy;
using ::boost::filesystem::path;
using ::boost::none;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;

const string data_dirs_options_block::DATA_OPTION_PATH_VARNAME   = "<path>";
const string data_dirs_options_block::DATA_OPTION_PREFIX_VARNAME = "<pre>";
const string data_dirs_options_block::DATA_OPTION_SUFFIX_VARNAME = "<suf>";

/// \brief For each of the options (path, prefix, suffix), the suffixes that are appended to the names of the file types to produce the option name
const data_option_str_map data_dirs_options_block::DATA_OPTION_SUFFIXES = {
	{ data_option::PATH,   "-path"   },
	{ data_option::PREFIX, "-prefix" },
	{ data_option::SUFFIX, "-suffix" }
};

/// \brief For each of the options (path, prefix, suffix), the start of the description up to the point of the data type name
const data_option_str_map data_dirs_options_block::DATA_OPTION_VARNAME = {
	{ data_option::PATH,   DATA_OPTION_PATH_VARNAME   },
	{ data_option::PREFIX, DATA_OPTION_PREFIX_VARNAME },
	{ data_option::SUFFIX, DATA_OPTION_SUFFIX_VARNAME }
};

/// \brief For each of the options (path, prefix, suffix), the start of the description up to the point of the data type name
const data_option_str_map data_dirs_options_block::DATA_OPTION_DESCRIPTION_START = {
	{ data_option::PATH,   "Search for "                                                                            },
	{ data_option::PREFIX, "Prepend the prefix " + DATA_OPTION_PREFIX_VARNAME + " to a protein's name to form its " },
	{ data_option::SUFFIX, "Append the suffix "  + DATA_OPTION_SUFFIX_VARNAME + " to a protein's name to form its " }
};

/// \brief For each of the options (path, prefix, suffix), the end of the description from the point of the data type name
const data_option_str_map data_dirs_options_block::DATA_OPTION_DESCRIPTION_END = {
	{ data_option::PATH,   " files using the path " + DATA_OPTION_PATH_VARNAME },
	{ data_option::PREFIX, " filename"                                         },
	{ data_option::SUFFIX, " filename"                                         }
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
void data_dirs_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                    const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                    ) {
	// Loop over each of the data option types (path, prefix, suffix) and grab the data_option and suffix string
	for (const data_option_str_pair &data_option_suffix_pair : DATA_OPTION_SUFFIXES) {
		const data_option &the_data_option    = data_option_suffix_pair.first;
		const string      &data_option_suffix = data_option_suffix_pair.second;

		// Loop over each of the data file types and grab the data_file and name string
		for (const data_file_str_pair &data_file_name_pair : data_dirs_spec::DATA_FILE_NAMES) {
			const data_file &the_data_file  = data_file_name_pair.first;
			const string    &data_file_name = data_file_name_pair.second;

			// Tidy up the file type name and append the option suffix to get the full option name (eg "PDB" -> "pdb", "pdb" + "-path" -> "pdb-path")
			const string     tidy_file_name = replace_all_copy(replace_all_copy(to_lower_copy(data_file_name), "_", "-"), " ", "-");
			const string     option_name    = tidy_file_name + data_option_suffix;

			// Make a notifier function object that will take a string and set it in the correct part of the_data_dirs_spec
			const auto the_notifier = [&] (const string &x) { the_data_dirs_spec.set_value_of_option_and_data_file( the_data_option, the_data_file, x ); };

			// Grab the default, description start and description end and join them together to get a description for the option
			const string     default_value = data_dirs_spec::DATA_FILE_TYPE_OPTION_DEFAULTS.at( the_data_file ).at( the_data_option );
			const string    &desc_start    = DATA_OPTION_DESCRIPTION_START.at(  the_data_option );
			const string    &desc_end      = DATA_OPTION_DESCRIPTION_END.at(    the_data_option );
			const string     description   = desc_start + data_file_name + desc_end;
			const string    &varname       = DATA_OPTION_VARNAME.at( the_data_option );

			// Add a new option, using the accumulated data
			prm_desc.add_options()
				(
					option_name.c_str(),
					value<string>()
						->value_name   ( varname       )
						->notifier     ( the_notifier  )
						->default_value( default_value ),
					description.c_str()
				);
		}
	}
}

/// \brief Add any hidden options to the provided options_description
void data_dirs_options_block::do_add_hidden_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                   const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                   ) {
	const string rootdir_valname = "<dir>";
	const auto cath_root_dir_notifier = [&] (const path &x) { the_data_dirs_spec.set_cath_root_dir( x ); };
	prm_desc.add_options()
		(
			PO_CATH_ROOT_DIR.c_str(),
			value<path>()
				->notifier  ( cath_root_dir_notifier )
				->value_name( rootdir_valname        ),
			( "Find sub-directories of standard names (\"pdb\", \"dssp\", \"wolf\", \"sec\") in root directory " + rootdir_valname ).c_str()
		);
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// At present, this always accepts
str_opt data_dirs_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                   ) const {
	const auto &cath_root_dir = the_data_dirs_spec.get_cath_root_dir();
	if ( ! cath_root_dir.empty() &&  ! is_acceptable_input_dir( cath_root_dir ) ) {
		return "CATH root directory \"" + cath_root_dir.string() + "\" is not a valid input directory";
	}

	return none;
}

/// \brief Return all options names for this block
str_vec data_dirs_options_block::do_get_all_options_names() const {
	return {
		data_dirs_options_block::PO_CATH_ROOT_DIR,
	};
}

/// \brief Getter for the data_dirs_spec
const data_dirs_spec & data_dirs_options_block::get_data_dirs_spec() const {
	return the_data_dirs_spec;
}
