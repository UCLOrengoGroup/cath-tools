/// \file
/// \brief The single_pair class definitions

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

#include "single_pair.hpp"

#include "common/clone/make_uptr_clone.hpp"
#include "scan/detail/scan_type_aliases.hpp"
#include "scan/res_pair_keyer/res_pair_keyer.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_from_phi_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_from_psi_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_index_dirn_keyer_part.hpp"
// #include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_orient_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_to_phi_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_to_psi_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_x_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_y_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_z_keyer_part.hpp"
#include "scan/scan_action/record_scores_scan_action.hpp"
#include "scan/scan_index.hpp"
#include "scan/scan_policy.hpp"
#include "scan/scan_query_set.hpp"
#include "scan/scan_stride.hpp"
#include "scan/scan_tools/scan_metrics.hpp"
#include "structure/geometry/angle.hpp"

#include <chrono>

using namespace cath::common;
using namespace cath::geom;
using namespace cath::scan;
using namespace cath::scan::detail;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<scan_type> single_pair::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
///
/// This can be used to scan all of one protein_list against all of another or all of a protein_list
/// against itself (by passing the same protein_list as arg_query_protein_list and arg_match_protein_list)

pair<record_scores_scan_action, scan_metrics> single_pair::do_perform_scan(const protein_list &arg_query_protein_list, ///< TODOCUMENT,
                                                                           const protein_list &arg_match_protein_list  ///< TODOCUMENT
                                                                           ) const {
	if ( arg_query_protein_list.size() != 1 || arg_match_protein_list.size() != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to perform single_pair scan because the two protein_lists don't both contain one entry"));
	}

	const auto angle_radius = make_angle_from_degrees<detail::angle_base_type>( 120 );
	const auto the_scan_policy = make_scan_policy(
		make_res_pair_keyer(
			res_pair_from_phi_keyer_part  { angle_radius },
			res_pair_from_psi_keyer_part  { angle_radius },
			res_pair_to_phi_keyer_part    { angle_radius },
			res_pair_to_psi_keyer_part    { angle_radius },
			res_pair_index_dirn_keyer_part{},
//			res_pair_orient_keyer_part    {},
			res_pair_view_x_keyer_part    { 12.65f },
			res_pair_view_y_keyer_part    { 12.65f },
			res_pair_view_z_keyer_part    { 12.65f }
		),
		make_default_quad_criteria(),
		scan_stride{ 4, 4, 2, 2 }
	);

	const auto the_query_set = make_scan_query_set( the_scan_policy, arg_query_protein_list );
	const auto the_index     = make_scan_index    ( the_scan_policy, arg_match_protein_list );

	record_scores_scan_action the_action(
		arg_query_protein_list.size(),
		arg_match_protein_list.size()
	);

	const auto scan_duration  = the_query_set.do_magic( the_index, the_action );
//	BOOST_LOG_TRIVIAL( warning ) << "Did magic - took " << do_magic_durn << "s (" << 1.0 / do_magic_durn << "/s)";

	const scan_metrics the_metrics{
		the_query_set.get_structures_build_durn_and_size(),
		the_query_set.get_index_build_durn_and_size(),
		the_index.get_structures_build_durn_and_size(),
		the_index.get_index_build_durn_and_size(),
		scan_duration
	};

	return make_pair( the_action, the_metrics );
}
