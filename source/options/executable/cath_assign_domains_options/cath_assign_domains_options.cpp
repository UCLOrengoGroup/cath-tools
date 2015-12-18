/// \file
/// \brief The cath_assign_domains_options class definitions

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

#include "cath_assign_domains_options.h"

#include <boost/program_options.hpp>

#include "alignment/alignment.h"
#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.h"
#include "common/logger.h"
#include "common/type_aliases.h"
#include "display/display_spec/display_spec.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "options/acquirer/alignment_acquirer/alignment_acquirer.h"
#include "options/acquirer/pdbs_acquirer/pdbs_acquirer.h"
#include "options/acquirer/selection_policy_acquirer/selection_policy_acquirer.h"
#include "options/acquirer/superposition_acquirer/align_based_superposition_acquirer.h"
#include "options/outputter/alignment_outputter/alignment_outputter.h"
#include "options/outputter/alignment_outputter/alignment_outputter_list.h"
#include "options/outputter/superposition_outputter/superposition_outputter.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"
#include "superposition/superposition_context.h"

#include <iostream>

namespace boost { namespace program_options { class variables_map; } }

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::lexical_cast;
using boost::ptr_vector;

/// The name of the program that uses this executable_options
const string cath_assign_domains_options::PROGRAM_NAME("cath-assign-domains");

/// Return the name of the program that uses this executable_options
string cath_assign_domains_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief Review all specified options and return a string containing any errors or a help string
///        (possibly using a description of all visible options)
///
/// This is a concrete definition of a virtual method that's pure in executable_options
///
/// This should only be called by executable_options, as the last step of the parse_options()
/// method, after all real parsing has completed.
///
/// \pre The options should have been parsed
///
/// \returns Any error/help string arising from the newly specified options
///          or an empty string if there aren't any
string cath_assign_domains_options::do_update_error_or_help_string(const options_description &arg_visible_program_options ///< The full options_description of visible options
                                                                   ) const {
	// If help was requested, then provide it
	if (the_misc_options_block.get_help()) {
		return the_misc_options_block.get_help_string(arg_visible_program_options, get_help_prefix_string(), get_help_suffix_string());
	}

	// If version information was requested, then provide it
	if (the_misc_options_block.get_version()) {
		return the_misc_options_block.get_version_string(get_program_name(), "This superposes protein structures.");
	}

	if ( get_rbf_svm_file().empty() ) {
		return "Please specify an SVM-light RBF model file\n" + get_standard_usage_error_string();
	}
	if ( get_data_data_file().empty() ) {
		return "Please specify a PRC/SSAP data files file\n" + get_standard_usage_error_string();
	}
	if ( get_sf_of_dom_file().empty() ) {
		return "Please specify a superfamily of domain file\n" + get_standard_usage_error_string();
	}

	return "";
}

/// \brief Get a string to prepend to the standard help
string cath_assign_domains_options::get_help_prefix_string() {
	ostringstream help_ss;
	help_ss << "Usage: " << PROGRAM_NAME << " [options]" << endl;
	help_ss << "Form plan of CATH assignments for query domains based on their PRC and SSAP results";
	return help_ss.str();
}

/// \brief Get a string to append to the standard help (just empty here)
string cath_assign_domains_options::get_help_suffix_string() {
	return "";
}

/// \brief Default-ctor that adds the required options_blocks to the parent executable_options class
cath_assign_domains_options::cath_assign_domains_options() {
	super::add_options_block( the_misc_options_block                );
	super::add_options_block( the_cath_assign_domains_options_block );
}

/// \brief Getter for the SVM-light RBF model file
const path & cath_assign_domains_options::get_rbf_svm_file() const {
	return the_cath_assign_domains_options_block.get_rbf_svm_file();
}

/// \brief Getter for the file containing the list of PRC/SSAP data files
const path & cath_assign_domains_options::get_data_data_file() const {
	return the_cath_assign_domains_options_block.get_data_data_file();
}

/// \brief Getter for the file containing superfamily of domain
const path & cath_assign_domains_options::get_sf_of_dom_file() const {
	return the_cath_assign_domains_options_block.get_sf_of_dom_file();
}

/// \brief Getter for the list of CATH nodes forbidden for assignment
const str_vec & cath_assign_domains_options::get_forbidden_nodes() const {
	return the_cath_assign_domains_options_block.get_forbidden_nodes();
}

