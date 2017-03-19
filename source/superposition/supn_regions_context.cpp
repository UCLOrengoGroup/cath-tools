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

#include <boost/optional.hpp>

#include "chopping/region/region.hpp"
#include "common/algorithm/transform_build.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "superposition/superposition_content_spec.hpp"

using namespace cath::chop;
using namespace cath::common;

using boost::none;
using std::set;

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
