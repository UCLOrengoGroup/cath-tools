/// \file
/// \brief The vcie_match_criteria class definitions

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

#include "vcie_match_criteria.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "structure/geometry/quat_rot.hpp"

using namespace boost::algorithm;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::index::detail;
using namespace std;

using boost::algorithm::is_any_of;
using boost::numeric_cast;

/// \brief Ctor to fully populate a vcie_match_criteria
vcie_match_criteria::vcie_match_criteria(const bool            &arg_require_matching_directions,    ///< Whether the two vcies are required to have matching directions
                                         const index_type      &arg_minimum_index_distance,         ///< The minimum distance required between each vcie's from_index and its to_index
                                         const view_base_type  &arg_maximum_squared_distance,       ///< The maximum squared distance permissible between the two vcies' views
                                         const angle_type      &arg_maximum_frame_angle_difference, ///< The maximum angle permissible between the two vcies' frames
                                         angle_type             arg_maximum_phi_angle_difference,   ///< The maximum angle permissible between the two vcies' phi angles
                                         angle_type             arg_maximum_psi_angle_difference    ///< The maximum angle permissible between the two vcies' psi angles
                                         ) : require_matching_directions    ( arg_require_matching_directions               ),
                                             minimum_index_distance         ( arg_minimum_index_distance                    ),
                                             maximum_squared_distance       ( arg_maximum_squared_distance                  ),
                                             maximum_frame_angle_distance_1 (
                                             	distance_1_of_angle<frame_quat_rot_type>(
                                             		arg_maximum_frame_angle_difference
                                             	)
                                             ),
                                             maximum_frame_angle_difference ( arg_maximum_frame_angle_difference            ),
                                             maximum_phi_angle_difference   ( std::move( arg_maximum_phi_angle_difference ) ),
                                             maximum_psi_angle_difference   ( std::move( arg_maximum_psi_angle_difference ) ) {
}

/// \brief Factory function to construct the default vcie_match_criteria
vcie_match_criteria cath::index::detail::make_default_vcie_match_criteria() {
	return vcie_match_criteria(
		true,
		11,
		40.0,
		make_angle_from_degrees<angle_base_type>( 22.5 ),
		make_angle_from_degrees<angle_base_type>( 67.5 ),
		make_angle_from_degrees<angle_base_type>( 67.5 )
	);
}

/// \brief Factory function to construct a vcie_match_criteria from a string
///
/// \pre The string must be of format like:
///   `"dist_co=12,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"`
vcie_match_criteria cath::index::detail::parse_vcie_match_criteria(const string &arg_string ///< The string to specify the properties of the vcie_match_criteria to be built
                                                                   ) {
	// Build a map of field name to value (represented as double)
	const string  spaces_stripped_string = erase_all_copy( arg_string, " " );
	const str_vec parts = split_build<str_vec>( spaces_stripped_string, is_any_of( "," ) );
	str_doub_map values;
	for (const string &part : parts) {
		str_vec halves = split_build<str_vec>( part, is_any_of( "=" ) );
		if ( halves.size() != 2 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Part of vcie_match_criteria string didn't have two halves"));
		}
		const string name = ends_with( halves.front(), "_co" ) ? erase_tail_copy( halves.front(), 3 ) : halves.front();
		values[ name ] = stod( halves[ 1 ] );
	}

	// Extract the required values
	const double &dist       = values.at( "dist"       );
	const double &dirn       = values.at( "dirn"       );
	const double &index_dist = values.at( "index_dist" );
	const double &frame_ang  = values.at( "frame_ang"  );
	const double &phi_ang    = values.at( "phi_ang"    );
	const double &psi_ang    = values.at( "psi_ang"    );

	// Build a vcie_match_criteria from the specified values
	return {
		( dirn == 0.0 ),
		numeric_cast< index_type     >( -index_dist ),
		numeric_cast< view_base_type >( dist        ),
		make_angle_from_degrees<angle_base_type>( frame_ang ),
		make_angle_from_degrees<angle_base_type>( phi_ang   ),
		make_angle_from_degrees<angle_base_type>( psi_ang   )
	};
}

/// \brief Build a standard set of vcie_match_criteria objects
///
/// This is implemented by passing a bunch of strings through parse_vcie_match_criteria()
vcie_match_criteria_vec cath::index::detail::get_standard_vcie_match_criterias() {
	const str_vec vcie_match_criteria_strings = {
		"dist_co=12,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=12,dirn_co=0,index_dist_co=-16,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=12,dirn_co=0,index_dist_co=-21,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=12,dirn_co=0,index_dist_co=-31,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-16,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-21,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-26,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=16,dirn_co=0,index_dist_co=-31,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-16,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=20,dirn_co=0,index_dist_co=-16,frame_ang_co=45,   phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=24,dirn_co=0,index_dist_co=-16,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-16,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-16,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=24,dirn_co=0,index_dist_co=-21,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=28,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-16,frame_ang_co=45,   phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-16,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=28,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=28,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=45,  psi_ang_co=22.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=32,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=32,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=45,  psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=112.5" ,
		"dist_co=36,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=135"   ,
		"dist_co=36,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=36,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=36,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=22.5, phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=22.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=112.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=135"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-11,frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=112.5,phi_ang_co=180, psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=112.5,phi_ang_co=90,  psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=135,  phi_ang_co=180, psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=157.5,phi_ang_co=180, psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=180,  phi_ang_co=180, psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=112.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=135"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=157.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=67.5, phi_ang_co=67.5,psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=112.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=135"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=157.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=90,   phi_ang_co=67.5,psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-1, frame_ang_co=90,   phi_ang_co=90,  psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=67.5,psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=112.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=135"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=157.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=180"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=45"    ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=45,   phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=112.5" ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=135"   ,
		"dist_co=40,dirn_co=0,index_dist_co=-6, frame_ang_co=67.5, phi_ang_co=90,  psi_ang_co=67.5"  ,
		"dist_co=8,dirn_co=0,index_dist_co=-11, frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=8,dirn_co=0,index_dist_co=-16, frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
		"dist_co=8,dirn_co=0,index_dist_co=-31, frame_ang_co=22.5, phi_ang_co=22.5,psi_ang_co=22.5"  ,
	};

	// Return the result of performing parse_vcie_match_criteria() on each string
	return transform_build<vcie_match_criteria_vec>(
		vcie_match_criteria_strings,
		[] (const string &x) { return parse_vcie_match_criteria( x ); }
	);
}

/// \brief Simple insertion operator for vcie_match_criteria
///
/// \relates vcie_match_criteria
ostream & cath::index::detail::operator<<(ostream                   &arg_os,                 ///< The ostream to which the vcie_match_criteria should be output
                                          const vcie_match_criteria &arg_vcie_match_criteria ///< The vcie_match_criteria to output
                                          ) {
	arg_os << "vcie_match_criteria_vec[";
	arg_os << "require_matching_directions:"    << boolalpha << arg_vcie_match_criteria.get_require_matching_directions();
	arg_os << ",minimum_index_distance:"                     << arg_vcie_match_criteria.get_minimum_index_distance();
	arg_os << ",maximum_squared_distance:"                   << arg_vcie_match_criteria.get_maximum_squared_distance();
	arg_os << ",maximum_frame_angle_difference:"             << arg_vcie_match_criteria.get_maximum_frame_angle_difference();
	arg_os << ",maximum_phi_angle_difference:"               << arg_vcie_match_criteria.get_maximum_phi_angle_difference();
	arg_os << ",maximum_psi_angle_difference:"               << arg_vcie_match_criteria.get_maximum_psi_angle_difference();
	arg_os << "]";
	return arg_os;
}

