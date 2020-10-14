/// \file
/// \brief The options_block class definitions

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

#include "options_block.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/program_options/variables_map_contains.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::adaptors::filtered;
using ::boost::algorithm::any_of;
using ::boost::filesystem::is_empty;
using ::boost::filesystem::path;
using ::boost::numeric_cast;
using ::boost::program_options::options_description;
using ::boost::program_options::variables_map;

/// \brief A string to use to separate (valid values and their descriptions) from each other
const string options_block::SUB_DESC_SEPARATOR = "\n   ";

/// \brief The string to use to separate valid values' names from their descriptions
const string options_block::SUB_DESC_PAIR_SEPARATOR = " - ";

/// \brief TODOCUMENT
void options_block::do_add_hidden_options_to_description(options_description &/*prm_desc*/,       ///< TODOCUMENT
                                                         const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                         ) {
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<options_block> options_block::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief Add all this block's visible options to the specified options_description
void options_block::add_visible_options_to_description(options_description &prm_desc, ///< The options_description to which the block's hidden options should be added
                                                       const size_t &prm_line_length  ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                       ) {
	do_add_visible_options_to_description( prm_desc, prm_line_length );
}

/// \brief Add all this block's hidden options to the specified options_description
void options_block::add_hidden_options_to_description(options_description &prm_desc, ///< The options_description to which the block's hidden options should be added
                                                      const size_t &prm_line_length  ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                      ) {
	do_add_hidden_options_to_description( prm_desc, prm_line_length );
}

/// \brief A method that uses the concrete class's methods to construct an options description
///        (using the specified line_length) of the hidden options
options_description options_block::get_all_options_description(const size_t &prm_line_length ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                               ) {
	options_description desc{ do_get_block_name(), numeric_cast<unsigned int>( prm_line_length ) };
	do_add_visible_options_to_description( desc, prm_line_length );
	add_hidden_options_to_description ( desc, prm_line_length );
	return desc;
}

/// \brief A method that uses the concrete class's methods to construct an options description
///        (using the specified line_length) of the visible options
///
/// The Boost program_options documentation seems to be a bit out of sync with the actual code.
/// In particular, line_length isn't (well) documented.
options_description options_block::get_visible_options_description(const size_t &prm_line_length ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                   ) {
	const string block_name = do_get_block_name();
	options_description desc( block_name, numeric_cast<unsigned int>( prm_line_length ) );
	do_add_visible_options_to_description( desc, prm_line_length );
	return desc;
}

/// \brief A method that uses the concrete class's methods to construct an options description
///        (using the specified line_length) of the hidden options
///
/// The Boost program_options documentation seems to be a bit out of sync with the actual code.
/// In particular, line_length isn't (well) documented.
options_description options_block::get_hidden_options_description(const size_t &prm_line_length ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                  ) {
	const string block_name = do_get_block_name();
	options_description desc( block_name + " [hidden options]", numeric_cast<unsigned int>( prm_line_length )  );
	add_hidden_options_to_description( desc, prm_line_length );
	return desc;
}

/// \brief An NVI pass-through to a method that returns a string describing any problems or "" if none
str_opt options_block::invalid_string(const variables_map &prm_variables_map ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                      ) const {
	return do_invalid_string( prm_variables_map );
}

/// \brief An NVI pass-through to a method that returns the options_block's full list of options names
str_vec options_block::get_all_options_names() const {
	return do_get_all_options_names();
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_output_file(const path &prm_output_file ///< TODOCUMENT
                                              ) {
	if (exists(prm_output_file) && !is_regular_file(prm_output_file)) {
		return false;
	}
	if (prm_output_file.has_parent_path() && !is_directory(prm_output_file.parent_path())) {
		return false;
	}
	return true;
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_input_file(const path &prm_output_file, ///< TODOCUMENT
                                             const bool &prm_allow_empty  ///< TODOCUMENT
                                             ) {
	using ::boost::filesystem::is_empty;

	if (!exists(prm_output_file)) {
		return false;
	}
	if (!prm_allow_empty && is_empty(prm_output_file)) {
		return false;
	}
	return true;
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_output_dir(const path &prm_output_dir ///< TODOCUMENT
                                             ) {
	return ( exists( prm_output_dir ) && is_directory( prm_output_dir ) );
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_input_dir(const path &prm_input_dir ///< TODOCUMENT
                                            ) {
	return ( exists( prm_input_dir ) && is_directory( prm_input_dir ) );
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_executable(const path &prm_output_file ///< TODOCUMENT
                                             ) {
	if (prm_output_file.filename() == ".." || prm_output_file.filename() == ".") {
//		cerr << "Here0" << endl;
		return false;
	}
	if (prm_output_file.has_root_directory() || prm_output_file.has_root_name()) {
//		cerr << "Here1" << endl;
		if (!exists(prm_output_file) && !is_regular_file(prm_output_file)) {
//			cerr << "Here2" << endl;
			return false;
		}
	}
	return true;
}

/// \brief Return whether the specified variables_map has specified the option with the specified name
///
/// \relates options_block
bool cath::opts::specifies_option(const variables_map &prm_vm,      ///< The variables_map to examine
                                  const string        &prm_opt_name ///< The name of the option of interest
                                  ) {
	return (
		contains( prm_vm, prm_opt_name )
		&&
		! prm_vm[ prm_opt_name ].defaulted()
	);
}

/// \brief Return the names of any of the "specified" options that have been specified in the "specified" variables_map
///
/// \relates options_block
str_vec cath::opts::specified_options(const variables_map &prm_vm,        ///< The variables_map to examine,
                                      const str_vec       &prm_opts_names ///< The names of the options of interest
                                      ) {
	return copy_build<str_vec>(
		prm_opts_names
			| filtered( [&] (const string &x) { return specifies_option( prm_vm, x ); } )
	);
}

/// \brief Return the names of any of the options in the "specified" options_block that have been specified in the "specified" variables_map
///
/// \relates options_block
str_vec cath::opts::specified_options_from_block(const variables_map &prm_vm,           ///< The variables_map to examine,
                                                 const options_block &prm_options_block ///< The options block
                                                 ) {
	return specified_options( prm_vm, prm_options_block.get_all_options_names() );
}

/// \brief Return whether the specified variables_map has specified any options with any of the specified names
///
/// \relates options_block
bool cath::opts::specifies_any_of_options(const variables_map &prm_vm,        ///< The variables_map to examine
                                          const str_vec       &prm_opts_names ///< The names of the options of interest
                                          ) {
	return any_of(
		prm_opts_names,
		[&] (const string &x) { return specifies_option( prm_vm, x ); }
	);
}

/// \brief Return whether the specified variables_map has specified any of the options in this block
///
/// \relates options_block
bool cath::opts::specifies_options_from_block(const variables_map &prm_vm,           ///< The variables_map to examine
                                              const options_block &prm_options_block ///< The options block
                                              ) {
	return specifies_any_of_options( prm_vm, prm_options_block.get_all_options_names() );
}
