/// \file
/// \brief The dssp_accessibility class definitions

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

#include "dssp_accessibility.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/accumulate_proj.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/size_t_literal.hpp"
#include "file/pdb/pdb.hpp"
#include "structure/geometry/coord.hpp"

#include <cmath>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sec::detail;

using boost::irange;
using boost::math::constants::pi;
using std::pair;
using std::plus;
using std::sqrt;

constexpr size_t dssp_ball_constants::NUMBER;
constexpr double dssp_ball_constants::RADIUS_N;
constexpr double dssp_ball_constants::RADIUS_CA;
constexpr double dssp_ball_constants::RADIUS_C;
constexpr double dssp_ball_constants::RADIUS_O;
constexpr double dssp_ball_constants::RADIUS_SIDE_ATOM;
constexpr double dssp_ball_constants::RADIUS_WATER;

/// \brief Make a load of points on teh surface of a unit sphere for accessibility calculations
coord_vec cath::sec::make_dssp_ball_points(const size_t &arg_number ///< The input number (just copying DSSP code here; the actual number of points is 2 * this + 1)
                                           ) {
	const double golden_ratio = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
	const size_t num_points   = 2_z * arg_number + 1_z;
	const int    num_as_int   = debug_numeric_cast<int>( arg_number );

	return transform_build<coord_vec>(
		irange( -num_as_int, num_as_int + 1 ),
		[&] (const int &x) {
			const double i         = debug_numeric_cast<double>( x );
			const double latitude  = asin( ( 2.0 * i ) / debug_numeric_cast<double>( num_points ) );
			const double longitude = fmod( i, golden_ratio ) * 2 * pi<double>() / golden_ratio;
			return coord{
				sin( longitude ) * cos( latitude  ),
				cos( longitude ) * cos( latitude  ),
				sin( latitude  )
			};
		}
	);
}

/// \brief Get the radius of the DSSP accessibility sphere associated with the specified coarse_element_type
double cath::sec::get_dssp_access_radius_without_water(const coarse_element_type &arg_coarse_element_type ///< The coarse_element_type to query
                                                       ) {
	switch ( arg_coarse_element_type ) {
		case( coarse_element_type::NITROGEN     ) : { return dssp_ball_constants::RADIUS_N;         }
		case( coarse_element_type::CARBON_ALPHA ) : { return dssp_ball_constants::RADIUS_CA;        }
		case( coarse_element_type::CARBON       ) : { return dssp_ball_constants::RADIUS_C;         }
		case( coarse_element_type::OXYGEN       ) : { return dssp_ball_constants::RADIUS_O;         }
		case( coarse_element_type::CARBON_BETA  ) : { return dssp_ball_constants::RADIUS_SIDE_ATOM; }
		case( coarse_element_type::NON_CORE     ) : { return dssp_ball_constants::RADIUS_SIDE_ATOM; }
		default : {

			return dssp_ball_constants::RADIUS_SIDE_ATOM;
		}
	}
}

/// \brief Get the radius of the DSSP accessibility sphere associated with the coarse_element_type of the specified pdb_atom
double cath::sec::get_dssp_access_radius_without_water(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                                       ) {
	return get_dssp_access_radius_without_water( get_coarse_element_type( arg_pdb_atom ) );
}

/// \brief Get the radius of the DSSP accessibility sphere plus water sphere associated with the specified coarse_element_type
///
/// This is the distance from the centre of a water sphere to the centre of a just-touching atom sphere of the relevant type
double cath::sec::get_dssp_access_radius_with_water(const coarse_element_type &arg_coarse_element_type ///< The coarse_element_type to query
                                                    ) {
	return get_dssp_access_radius_without_water( arg_coarse_element_type ) + dssp_ball_constants::RADIUS_WATER;
}

/// \brief Get the radius of the DSSP accessibility sphere plus water sphere associated with the coarse_element_type of the specified pdb_atom
///
/// This is the distance from the centre of a water sphere to the centre of a just-touching atom sphere of the relevant type
double cath::sec::get_dssp_access_radius_with_water(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                                    ) {
	return get_dssp_access_radius_without_water( arg_pdb_atom ) + dssp_ball_constants::RADIUS_WATER;
}

/// \brief Return whether the two specified atoms overlap at the specified ball point of the first
bool cath::sec::access_overlap(const pdb_atom &arg_atom_lhs,   ///< The first atom to compare
                               const coord    &arg_ball_point, ///< The ball point to use for the first atom
                               const pdb_atom &arg_atom_rhs    ///< The second atom to compare
                               ) {
	const double radius_rhs = get_dssp_access_radius_with_water( arg_atom_rhs );
	return (
		squared_distance_between_points(
			arg_atom_lhs.get_coord() + arg_ball_point,
			arg_atom_rhs.get_coord()
		)
		<=
		( radius_rhs * radius_rhs )
	);
}

/// \brief Get the accessibility count for the specified atom in the context of the specified PDB
size_t cath::sec::get_accessibility_count(const pdb_atom &arg_pdb_atom, ///< The PDB atom for which the accessibility should be calculated
                                          const pdb      &arg_pdb,      ///< The PDB in which the accesibility should be calculated 
                                          const size_t   &arg_number    ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                          ) {
	const coord_vec ball_points = make_dssp_ball_points( arg_number );
	const double    radius      = get_dssp_access_radius_with_water( arg_pdb_atom );

	size_t count = 0;
	for (const coord &orig_ball_point : ball_points) {
		const coord ball_point = radius * orig_ball_point;
		bool found_overlap = false;
		for (const pdb_residue &the_res : arg_pdb) {
			for (const pdb_atom &the_atom : the_res) {
				if ( arg_pdb_atom .get_coord() != the_atom.get_coord() ) {
					found_overlap = found_overlap || access_overlap( arg_pdb_atom, ball_point, the_atom );
				}
			}
		}
		if ( ! found_overlap ) {
			++count;
		}
	}
	return count;
}

/// \brief Calculate the accessibility fraction for the specified atom in the specified PDB
///
/// The fraction is the fraction of the points on the sphere that are found to be accessible
double cath::sec::get_accessibility_fraction(const pdb_atom &arg_pdb_atom, ///< The PDB atom for which the accessibility should be calculated
                                             const pdb      &arg_pdb,      ///< The PDB in which the accesibility should be calculated
                                             const size_t   &arg_number    ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                             ) {
	const coord_vec ball_points = make_dssp_ball_points( arg_number );

	return debug_numeric_cast<double>( get_accessibility_count( arg_pdb_atom, arg_pdb, arg_number ) )
	     / debug_numeric_cast<double>( make_dssp_ball_points( arg_number ).size() );
}

/// \brief Calculate the accessibility surface-area for the specified atom in the specified PDB
double cath::sec::get_accessibility_surface_area(const pdb_atom &arg_pdb_atom, ///< The PDB atom for which the accessibility should be calculated
                                                 const pdb      &arg_pdb,      ///< The PDB in which the accesibility should be calculated
                                                 const size_t   &arg_number    ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                                 ) {
	const double radius              = get_dssp_access_radius_with_water( arg_pdb_atom );
	const double sphere_surface_area = 4.0 * pi<double>() * radius * radius;

	return get_accessibility_fraction( arg_pdb_atom, arg_pdb, arg_number ) * sphere_surface_area;
}

/// \brief Calculate the accessibility surface-area for the specified residue in the specified PDB
double cath::sec::get_accessibility_surface_area(const pdb_residue &arg_pdb_residue, ///< The PDB residue for which the accessibility should be calculated
                                                 const pdb         &arg_pdb,         ///< The PDB in which the accesibility should be calculated
                                                 const size_t      &arg_number       ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                                 ) {
	return accumulate_proj(
		arg_pdb_residue,
		0.0,
		plus<>{},
		[&] (const pdb_atom &x) { return get_accessibility_surface_area( x, arg_pdb, arg_number ); }
	);
}

/// \brief Calculate the per-residue surface-area accessibilities for the specified PDBs
doub_vec cath::sec::calc_accessibilities(const pdb    &arg_pdb,   ///< The PDB in which the accesibility should be calculated
                                         const size_t &arg_number ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                         ) {
	return transform_build<doub_vec>(
		arg_pdb,
		[&] (const pdb_residue &x) { return get_accessibility_surface_area( x, arg_pdb, arg_number ); }
	);
}
