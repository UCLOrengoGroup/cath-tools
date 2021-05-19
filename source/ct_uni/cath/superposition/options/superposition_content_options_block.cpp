/// \file
/// \brief The superposition_content_options_block class definitions

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

#include "superposition_content_options_block.hpp"

#include <boost/algorithm/string/join.hpp>

#include "cath/common/boost_addenda/program_options/layout_values_with_descs.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/superposition/supn_regions_context.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath::sup;

using ::boost::algorithm::join;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<options_block> superposition_content_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string superposition_content_options_block::do_get_block_name() const {
	return "Superposition content";
}

/// \brief Add this block's options to the provided options_description
void superposition_content_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                                const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                                ) {
	const auto &sep     = SUB_DESC_SEPARATOR;
	const auto &sub_sep = SUB_DESC_PAIR_SEPARATOR;

	const string context_varname { "<context>" };
	const string dist_varname    { "<dist>" };

	const auto regions_context_notifier                 = [&] (const supn_regions_context &x) { the_spec.set_regions_context                ( x ); };
	const auto include_dna_within_distance_notifier     = [&] (const double               &x) { the_spec.set_include_dna_within_distance    ( x ); };
	const auto include_organic_within_distance_notifier = [&] (const double               &x) { the_spec.set_include_organic_within_distance( x ); };

	const str_vec supn_region_context_descs = layout_values_with_descs(
		all_supn_regions_contexts,
		[] (const supn_regions_context &x) { return to_string( x ); },
		&description_of_supn_regions_context,
		sub_sep
	);

	prm_desc.add_options()
		(
			string( PO_REGIONS_CONTEXT ).c_str(),
			value<supn_regions_context>()
				->value_name   ( context_varname                                      )
				->notifier     ( regions_context_notifier                             )
				->default_value( superposition_content_spec::DEFAULT_REGIONS_CONTEXT  ),
			::fmt::format(
				"Show the alignment regions in the context {}, one of available options:{}{}",
				context_varname,
				sep,
				join( supn_region_context_descs, sep )
			).c_str()
		)
		(
			string( PO_INCLUDE_DNA_WITHIN_DISTANCE ).c_str(),
			value<double>()
				->value_name   ( dist_varname                                         )
				->notifier     ( include_dna_within_distance_notifier                 )
				->default_value( superposition_content_spec::DEFAULT_DNA_MAX_DIST     ),
			( "Show DNA within " + dist_varname + "Å of the alignment regions" ).c_str()
		)
		(
			string( PO_INCLUDE_ORGANIC_WITHIN_DISTANCE ).c_str(),
			value<double>()
				->value_name   ( dist_varname                                         )
				->notifier     ( include_organic_within_distance_notifier             )
				->default_value( superposition_content_spec::DEFAULT_ORGANIC_MAX_DIST ),
			( "Show organic molecules within " + dist_varname + "Å of the alignment regions" ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified superposition_content_options_block invalid
///        or nullopt otherwise
str_opt superposition_content_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                               ) const {
	return get_invalid_description( the_spec );
}

/// \brief Return all options names for this block
str_view_vec superposition_content_options_block::do_get_all_options_names() const {
	return {
		PO_REGIONS_CONTEXT,
		PO_INCLUDE_DNA_WITHIN_DISTANCE,
		PO_INCLUDE_ORGANIC_WITHIN_DISTANCE,
	};
}

/// \brief Getter for the superposition_content_spec that the superposition_content_options_block configures
const superposition_content_spec & superposition_content_options_block::get_superposition_content_spec() const {
	return the_spec;
}
