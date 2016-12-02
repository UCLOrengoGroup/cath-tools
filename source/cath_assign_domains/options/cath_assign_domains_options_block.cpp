/// \file
/// \brief The cath_assign_domains_options_block class definitions

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

#include "cath_assign_domains_options_block.hpp"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp"
#include "acquirer/pdbs_acquirer/istream_pdbs_acquirer.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "score/homcheck_tools/superfamily_of_domain.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::homcheck::detail;
using namespace cath::opts;
using namespace std;

using boost::algorithm::all_of;
using boost::algorithm::join;
using boost::assign::ptr_push_back;
using boost::filesystem::path;
using boost::none;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using boost::ptr_vector;

/// \brief The long option specifying the SVM-light RBF model file
const string cath_assign_domains_options_block::PO_SVMLIGHT_RBF_FILE( "svmlight-rbf-file" );

/// \brief The long option specifying the file containing the list of PRC/SSAP data files
const string cath_assign_domains_options_block::PO_FILELIST_FILE    ( "filelist-file"     );

/// \brief The long option specifying the superfamily of domain file
const string cath_assign_domains_options_block::PO_SF_OF_DOMAIN_FILE( "sf-of-domain-file" );

/// \brief The long option specifying the CATH nodes forbidden for assignment
const string cath_assign_domains_options_block::PO_FORBIDDEN_NODES  ( "forbidden-node"    );

/// \brief The default list of CATH nodes forbidden for assignment
const str_vec cath_assign_domains_options_block::DEFAULT_FORBIDDEN_NODES = {
	"2.105",
	"2.110",
	"2.115",
	"2.120",
	"2.130",
	"2.140"
};

/// \brief A standard do_clone method.
unique_ptr<options_block> cath_assign_domains_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string cath_assign_domains_options_block::do_get_block_name() const {
	return "CATH Assign Domains Specification";
}

/// \brief Add this block's options to the provided options_description
void cath_assign_domains_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                              ) {
	const string default_forbidden_nodes_str = join( DEFAULT_FORBIDDEN_NODES, ", " );
	arg_desc.add_options()
		( PO_SVMLIGHT_RBF_FILE.c_str(), value<path   >( &rbf_svm_file    ),                                                                        "File containing SVM-light RBF model for CATH assignment" )
		( PO_FILELIST_FILE.c_str(),     value<path   >( &data_data_file  ),                                                                        "File of data files (one line per query domain containing: ssap_results_file prc_results_file)" )
		( PO_SF_OF_DOMAIN_FILE.c_str(), value<path   >( &sf_of_dom_file  ),                                                                        "File containing up-to-date assignments (one line per domain containing: domain_id superfamily_id)" )
		( PO_FORBIDDEN_NODES.c_str(),   value<str_vec>( &forbidden_nodes )->default_value( DEFAULT_FORBIDDEN_NODES, default_forbidden_nodes_str ), "List of nodes to which automatic assignment is forbidden; specify option multiple times for multiple nodes\nRECOMMENDED: do not specify this option so that the default list of propeller architectures is used." );
}

str_opt cath_assign_domains_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                             ) const {
	if ( ! rbf_svm_file.empty()   && ! is_acceptable_input_file( rbf_svm_file   ) ) {
		return "SVM-light RBF model file " + rbf_svm_file.string() + " is not a valid input file";
	}
	if ( ! data_data_file.empty() && ! is_acceptable_input_file( data_data_file ) ) {
		return "List of PRC/SSAP data files file " + data_data_file.string() + " is not a valid input file";
	}
	if ( ! sf_of_dom_file.empty() && ! is_acceptable_input_file( sf_of_dom_file ) ) {
		return "Superfamily of domain file " + sf_of_dom_file.string() + " is not a valid input file";
	}
	if ( ! all_of( forbidden_nodes, is_valid_cath_node_id{} ) ) {
		return "Forbidden nodes list " + join( forbidden_nodes, ", ") + " is not a list of valid CATH nodes";
	}
	return none;
}

/// \brief Getter for the SVM-light RBF model file
const path & cath_assign_domains_options_block::get_rbf_svm_file() const {
	return rbf_svm_file;
}

/// \brief Getter for the file containing the list of PRC/SSAP data files
const path & cath_assign_domains_options_block::get_data_data_file() const {
	return data_data_file;
}

/// \brief Getter for the file containing superfamily of domain
const path & cath_assign_domains_options_block::get_sf_of_dom_file() const {
	return sf_of_dom_file;
}

/// \brief Getter for the list of CATH nodes forbidden for assignment
const str_vec & cath_assign_domains_options_block::get_forbidden_nodes() const {
	return forbidden_nodes;
}

