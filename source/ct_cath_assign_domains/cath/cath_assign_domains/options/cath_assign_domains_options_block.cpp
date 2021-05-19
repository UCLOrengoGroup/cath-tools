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

#include <array>
#include <filesystem>
#include <string_view>

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/string/join.hpp>

#include "cath/acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp"
#include "cath/acquirer/pdbs_acquirer/istream_pdbs_acquirer.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/score/homcheck_tools/superfamily_of_domain.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::homcheck::detail;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::algorithm::all_of;
using ::boost::algorithm::join;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::array;
using ::std::filesystem::path;
using ::std::literals::string_literals::operator""s;
using ::std::nullopt;
using ::std::string_view;

/// \brief The long option specifying the SVM-light RBF model file
constexpr string_view PO_SVMLIGHT_RBF_FILE( "svmlight-rbf-file" );

/// \brief The long option specifying the file containing the list of PRC/SSAP data files
constexpr string_view PO_FILELIST_FILE    ( "filelist-file"     );

/// \brief The long option specifying the superfamily of domain file
constexpr string_view PO_SF_OF_DOMAIN_FILE( "sf-of-domain-file" );

/// \brief The long option specifying the CATH nodes forbidden for assignment
constexpr string_view PO_FORBIDDEN_NODES  ( "forbidden-node"    );

/// \brief The default list of CATH nodes forbidden for assignment
str_vec DEFAULT_FORBIDDEN_NODES() {
	// clang-format off
	return {
		"2.105"s,
		"2.110"s,
		"2.115"s,
		"2.120"s,
		"2.130"s,
		"2.140"s
	};
	// clang-format on
}

/// \brief A standard do_clone method.
unique_ptr<options_block> cath_assign_domains_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string cath_assign_domains_options_block::do_get_block_name() const {
	return "CATH Assign Domains Specification";
}

/// \brief Add this block's options to the provided options_description
void cath_assign_domains_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                              const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                              ) {
	const str_vec default_forbidden_nodes = DEFAULT_FORBIDDEN_NODES();
	const string default_forbidden_nodes_str = join( default_forbidden_nodes, ", " );
	prm_desc.add_options()
		( string( PO_SVMLIGHT_RBF_FILE ).c_str(), value<path   >( &rbf_svm_file    ),                                                                        "File containing SVM-light RBF model for CATH assignment" )
		( string( PO_FILELIST_FILE ).c_str(),     value<path   >( &data_data_file  ),                                                                        "File of data files (one line per query domain containing: ssap_results_file prc_results_file)" )
		( string( PO_SF_OF_DOMAIN_FILE ).c_str(), value<path   >( &sf_of_dom_file  ),                                                                        "File containing up-to-date assignments (one line per domain containing: domain_id superfamily_id)" )
		( string( PO_FORBIDDEN_NODES ).c_str(),   value<str_vec>( &forbidden_nodes )->default_value( default_forbidden_nodes, default_forbidden_nodes_str ), "List of nodes to which automatic assignment is forbidden; specify option multiple times for multiple nodes\nRECOMMENDED: do not specify this option so that the default list of propeller architectures is used." );
}

/// \brief TODOCUMENT
str_opt cath_assign_domains_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
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
	return nullopt;
}

/// \brief Return all options names for this block
str_view_vec cath_assign_domains_options_block::do_get_all_options_names() const {
	return {
		PO_SVMLIGHT_RBF_FILE,
		PO_FILELIST_FILE,
		PO_SF_OF_DOMAIN_FILE,
		PO_FORBIDDEN_NODES,
	};
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

