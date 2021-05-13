/// \file
/// \brief The clust_mapping_options_block class definitions

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

#include "clust_mapping_options_block.hpp"

#include <boost/lexical_cast.hpp>

#include "cath/common/clone/make_uptr_clone.hpp"

using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath;

using ::boost::lexical_cast;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::nullopt;
using ::std::string;
using ::std::unique_ptr;

/// \brief The option name for the fraction that the overlap over the longest of two domains must exceed for them to be considered equivalent
const string clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL   { "min_equiv_dom_ol"   };

/// \brief The option name for the fraction of the old cluster's entries that must map to a map-from cluster for them to be considered equivalent
const string clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL { "min_equiv_clust_ol" };

/// \brief A standard do_clone method
unique_ptr<options_block> clust_mapping_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string clust_mapping_options_block::do_get_block_name() const {
	return "Mapping";
}

/// \brief Add this block's options to the provided options_description
void clust_mapping_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                        const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                        ) {
	using ::std::to_string;
	const string percent_varname   { "<percent>" };

	const auto min_equiv_dom_ol_notifier   = [&] (const double &x) { the_spec.set_min_equiv_dom_ol  ( x / 100.0 ); };
	const auto min_equiv_clust_ol_notifier = [&] (const double &x) { the_spec.set_min_equiv_clust_ol( x / 100.0 ); };

	prm_desc.add_options()
		(
			PO_MIN_EQUIV_DOM_OL.c_str(),
			value<double>()
				->value_name   ( percent_varname                                        )
				->notifier     ( min_equiv_dom_ol_notifier                              )
				->default_value( 100.0 * clust_mapping_spec::DEFAULT_MIN_EQUIV_DOM_OL   ),
			( "Define domain equivalence as: sharing more than " + percent_varname + R"(% of residues (over the longest domain))"
				+ "\n"
				+ "(where "
				+ percent_varname
				+ R"( must be ≥ )"
				+ lexical_cast<string>( 100.0 * clust_mapping_spec::MIN_MIN_EQUIV_DOM_OL               )
				+ ")" ).c_str()
		)
		(
			PO_MIN_EQUIV_CLUST_OL.c_str(),
			value<double>()
				->value_name   ( percent_varname                                        )
				->notifier     ( min_equiv_clust_ol_notifier                            )
				->default_value( 100.0 * clust_mapping_spec::DEFAULT_MIN_EQUIV_CLUST_OL ),
			( "Define cluster equivalence as: more than " + percent_varname
				+ R"(% of the map-from cluster's members having equivalents in the working cluster)"
				+ "\n"
				R"([and them being equivalent to > )"
				+ lexical_cast<string>( 100.0 * clust_mapping_spec::MIN_EQUIV_FRAC_OF_NEW_CLUST        )
				+ R"(% of the working cluster's entries and > )"
				+ lexical_cast<string>( 100.0 * clust_mapping_spec::MIN_EQUIV_FRAC_OF_NEW_CLUST_EQUIVS )
				+ R"(% of those that have an equivalence])" "\n"
				"(where "
				+ percent_varname
				+ R"( must be ≥ )"
				+ lexical_cast<string>( 100.0 * clust_mapping_spec::MIN_MIN_EQUIV_CLUST_OL             )
				+ R"(%))" ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified clust_mapping_options_block invalid
///        or nullopt otherwise
str_opt clust_mapping_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                       ) const {
	return nullopt;
}

/// \brief Return all options names for this block
str_vec clust_mapping_options_block::do_get_all_options_names() const {
	return {
		clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL,
		clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL,
	};
}

/// \brief Getter for the clust_mapping_spec that the clust_mapping_options_block configures
const clust_mapping_spec & clust_mapping_options_block::get_clust_mapping_spec() const {
	return the_spec;
}

/// \brief Return the names of the threshold options
///
/// \relates clust_mapping_options_block
str_vec cath::clust::clust_thresh_option_names() {
	return {
		clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL,
		clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL,
	};
}

/// \brief Return whether the any of the threshold options have been specified in the specified variables_map
///
/// \relates clust_mapping_options_block
bool cath::clust::specified_clust_thresh_options(const variables_map &prm_vm ///< The variables_map to examine
                                                 ) {
	return specifies_any_of_options( prm_vm, clust_thresh_option_names() );
}
