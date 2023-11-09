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

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

#include <fmt/core.h>

#include "cath/cath_assign_domains/options/cath_assign_domains_options.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::std::filesystem::path;
using ::std::literals::string_literals::operator""s;
using ::std::nullopt;
using ::std::string;
using ::std::string_view;

/// \brief Get the name of the program that uses this executable_options
string_view cath_assign_domains_options::do_get_program_name() const {
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
str_opt cath_assign_domains_options::do_get_error_or_help_string() const {
	if ( get_rbf_svm_file().empty() ) {
		return "Please specify an SVM-light RBF model file"s;
	}
	if ( get_data_data_file().empty() ) {
		return "Please specify a PRC/SSAP data files file"s;
	}
	if ( get_sf_of_dom_file().empty() ) {
		return "Please specify a superfamily of domain file"s;
	}

	return nullopt;
}

/// \brief Get a string to prepend to the standard help
string cath_assign_domains_options::do_get_help_prefix_string() const {
	return ::fmt::format( "Usage: {} [options]\n\n{}", PROGRAM_NAME, get_overview_string() );
}

/// \brief Get a string to append to the standard help (just empty here)
string cath_assign_domains_options::do_get_help_suffix_string() const {
	return "";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string cath_assign_domains_options::do_get_overview_string() const {
	return R"(Use an SVM model on SSAP+PRC data to form a plan for assigning the domains to CATH superfamilies/folds)";
}

/// \brief Default-ctor that adds the required options_blocks to the parent executable_options class
cath_assign_domains_options::cath_assign_domains_options() {
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
