/// \file
/// \brief The supn_regions_context class definitions

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

#include "supn_regions_context.hpp"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/optional.hpp>

#include "chopping/region/region.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "common/program_options/validator.hpp"
#include "superposition/superposition_content_spec.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::sup;
using namespace cath::sup::detail;

using boost::algorithm::to_lower;
using boost::any;
using boost::none;
using std::istream;
using std::map;
using std::ostream;
using std::set;
using std::string;

/// \brief Getter for a map from name to supn_regions_context
///
/// \relates supn_regions_context
map<string, supn_regions_context> supn_regions_context_by_name::get() {
	map<string, supn_regions_context> result;
	for (const auto &x : all_supn_regions_contexts) {
		result.emplace( to_string( x ), x );
	}
	return result;
}


/// \brief Get a list of all the supn_regions_context names
str_vec all_supn_regions_context_names::get() {
	str_vec result;
	for (const auto &x : all_supn_regions_contexts) {
		result.emplace_back( to_string( x ) );
	}
	return result;
}

/// \brief Get a copy of the specified regions that have been expanded according to
///        the specified supn_regions_context
///
/// \relates supn_regions_context
region_vec_opt cath::sup::get_regions_expanded_for_context(const region_vec           &arg_regions,             ///< The regions to expand
                                                           const supn_regions_context &arg_supn_regions_context ///< The supn_regions_context describing how (if at all) the regions should be expanded
                                                           ) {
	switch ( arg_supn_regions_context ) {
		case ( supn_regions_context::ALONE    ) : {
			return arg_regions;
		}
		case ( supn_regions_context::IN_CHAIN ) : {
			// \todo Come a later GCC, just use this:
			//     return sort_uniq_build<region_vec>(
			//     	arg_regions
			//     		| transformed( [] (const region &x) {
			//     			return expand_to_chain( x );
			//     		} ),
			//     	[] (const region &x, const region &y) {
			//     		return get_chain_label( x ) < get_chain_label( y );
			//     	}
			//     );
			auto expanded_region_chain_labels = transform_build<set<chain_label>>(
				arg_regions,
				[] (const region &x) {
					return get_chain_label( x );
				}
			);
			return region_vec(
				common::cbegin( expanded_region_chain_labels ),
				common::cend  ( expanded_region_chain_labels )
			);
		}
		case ( supn_regions_context::IN_PDB   ) : {
			return none;
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of supn_regions_context not recognised whilst in get_regions_expanded_for_context()"));
}

/// \brief Get a copy of the specified regions that have been expanded according to
///        the supn_regions_context implied by the specified superposition_content_spec
///
/// \relates superposition_content_spec
region_vec_opt cath::sup::get_regions_expanded_for_context(const region_vec                 &arg_regions,     ///< The regions to expand
                                                           const superposition_content_spec &arg_content_spec ///< The superposition_content_spec to query for the supn_regions_context
                                                           ) {
	return get_regions_expanded_for_context( arg_regions, arg_content_spec.get_regions_context() );
}

/// \brief Generate a string describing the specified supn_regions_context
///
/// \relates supn_regions_context
string cath::sup::to_string(const supn_regions_context &arg_context ///< The supn_regions_context to describe
                            ) {
	switch ( arg_context ) {
		case ( supn_regions_context::ALONE    ) : { return "alone"    ; }
		case ( supn_regions_context::IN_CHAIN ) : { return "in_chain" ; }
		case ( supn_regions_context::IN_PDB   ) : { return "in_pdb"   ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("supn_regions_context value not recognised"));
}

/// \brief Insert a description of the specified supn_regions_context into the specified ostream
///
/// \relates supn_regions_context
ostream & cath::sup::operator<<(ostream                    &arg_os,     ///< The ostream into which the description should be inserted
                                const supn_regions_context &arg_context ///< The supn_regions_context to describe
                                ) {
	arg_os << to_string( arg_context );
	return arg_os;
}

/// \brief Extract into the specified supn_regions_context from the specified stream
///
/// \relates supn_regions_context
istream & cath::sup::operator>>(istream              &arg_is,        ///< The stream from which the supn_regions_context should be extracted
                                supn_regions_context &arg_format_tag ///< The supn_regions_context to populate from the specified stream
                                ) {
	string input_string;
	arg_is >> input_string;
	to_lower( input_string );

	const auto all_layout_tags_by_name = supn_regions_context_by_name::get();
	arg_format_tag = all_layout_tags_by_name.at( input_string );
	return arg_is;
}

/// \brief Generate a string containing a description of the specified supn_regions_context
///
/// \relates supn_regions_context
string cath::sup::description_of_supn_regions_context(const supn_regions_context &arg_format_tag ///< The supn_regions_context to describe
                                                      ) {
	switch ( arg_format_tag ) {
		case ( supn_regions_context::ALONE    ) : { return "alone"                                           ; }
		case ( supn_regions_context::IN_CHAIN ) : { return "within the chain(s) in which the regions appear" ; }
		case ( supn_regions_context::IN_PDB   ) : { return "within the PDB in which the regions appear"      ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("supn_regions_context value not recognised"));
}

/// \brief Provide Boost program_options validation for supn_regions_context
void cath::sup::validate(any           &arg_value,         ///< The value to populate
                         const str_vec &arg_value_strings, ///< The string values to validate
                         supn_regions_context *, int) {
	arg_value = lex_castable_validator<supn_regions_context>::perform_validate( arg_value, arg_value_strings );
}
