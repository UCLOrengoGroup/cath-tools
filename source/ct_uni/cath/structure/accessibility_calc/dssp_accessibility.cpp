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
#include <boost/range/algorithm/count.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/accumulate_proj.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_store_helper.hpp"
#include "cath/scan/spatial_index/spatial_index.hpp"
#include "cath/structure/geometry/coord.hpp"

#include <cmath>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::geom;
using namespace ::cath::scan;
using namespace ::cath::scan::detail;
using namespace ::cath::sec::detail;

using ::boost::irange;
using ::boost::make_optional;
using ::boost::math::constants::pi;
using ::boost::none;
using ::boost::numeric_cast;
using ::boost::range::count;
using ::std::plus;
using ::std::sqrt;
using ::std::vector;

constexpr size_t dssp_ball_constants::NUMBER;
constexpr double dssp_ball_constants::RADIUS_N;
constexpr double dssp_ball_constants::RADIUS_CA;
constexpr double dssp_ball_constants::RADIUS_C;
constexpr double dssp_ball_constants::RADIUS_O;
constexpr double dssp_ball_constants::RADIUS_SIDE_ATOM;
constexpr double dssp_ball_constants::RADIUS_WATER;
constexpr double dssp_ball_constants::MAX_ATOM_DIST;

/// \brief Make a load of points on teh surface of a unit sphere for accessibility calculations
coord_vec cath::sec::make_dssp_ball_points(const size_t &prm_number ///< The input number (just copying DSSP code here; the actual number of points is 2 * this + 1)
                                           ) {
	const double golden_ratio = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
	const size_t num_points   = 2_z * prm_number + 1_z;
	const int    num_as_int   = debug_numeric_cast<int>( prm_number );

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
double cath::sec::get_dssp_access_radius_without_water(const coarse_element_type &prm_coarse_element_type ///< The coarse_element_type to query
                                                       ) {
	switch ( prm_coarse_element_type ) {
		case( coarse_element_type::NITROGEN     ) : { return dssp_ball_constants::RADIUS_N;         }
		case( coarse_element_type::CARBON_ALPHA ) : { return dssp_ball_constants::RADIUS_CA;        }
		case( coarse_element_type::CARBON       ) : { return dssp_ball_constants::RADIUS_C;         }
		case( coarse_element_type::OXYGEN       ) : { return dssp_ball_constants::RADIUS_O;         }
		case( coarse_element_type::CARBON_BETA  ) : { return dssp_ball_constants::RADIUS_SIDE_ATOM; }
		case( coarse_element_type::NON_CORE     ) : { return dssp_ball_constants::RADIUS_SIDE_ATOM; }
	}
	return dssp_ball_constants::RADIUS_SIDE_ATOM;
}

/// \brief Get the radius of the DSSP accessibility sphere associated with the coarse_element_type of the specified pdb_atom
double cath::sec::get_dssp_access_radius_without_water(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
                                                       ) {
	return get_dssp_access_radius_without_water( get_coarse_element_type( prm_pdb_atom ) );
}

/// \brief Get the radius of the DSSP accessibility sphere plus water sphere associated with the specified coarse_element_type
///
/// This is the distance from the centre of a water sphere to the centre of a just-touching atom sphere of the relevant type
double cath::sec::get_dssp_access_radius_with_water(const coarse_element_type &prm_coarse_element_type ///< The coarse_element_type to query
                                                    ) {
	return get_dssp_access_radius_without_water( prm_coarse_element_type ) + dssp_ball_constants::RADIUS_WATER;
}

/// \brief Get the radius of the DSSP accessibility sphere plus water sphere associated with the coarse_element_type of the specified pdb_atom
///
/// This is the distance from the centre of a water sphere to the centre of a just-touching atom sphere of the relevant type
double cath::sec::get_dssp_access_radius_with_water(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
                                                    ) {
	return get_dssp_access_radius_without_water( prm_pdb_atom ) + dssp_ball_constants::RADIUS_WATER;
}

/// \brief Return whether the two specified atoms overlap at the specified ball point of the first
bool cath::sec::access_overlap(const pdb_atom &prm_atom_lhs,   ///< The first atom to compare
                               const coord    &prm_ball_point, ///< The ball point to use for the first atom
                               const pdb_atom &prm_atom_rhs    ///< The second atom to compare
                               ) {
	const double radius_rhs = get_dssp_access_radius_with_water( prm_atom_rhs );
	return (
		squared_distance_between_points(
			prm_atom_lhs.get_coord() + prm_ball_point,
			prm_atom_rhs.get_coord()
		)
		<=
		( radius_rhs * radius_rhs )
	);
}

/// \brief Get the accessibility count for the specified atom in the context of the specified PDB
size_t cath::sec::get_accessibility_count(const pdb_atom &prm_pdb_atom, ///< The PDB atom for which the accessibility should be calculated
                                          const pdb      &prm_pdb,      ///< The PDB in which the accesibility should be calculated 
                                          const size_t   &prm_number    ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                          ) {
	const coord_vec ball_points = make_dssp_ball_points( prm_number );
	const double    radius      = get_dssp_access_radius_with_water( prm_pdb_atom );

	size_t count = 0;
	for (const coord &orig_ball_point : ball_points) {
		const coord ball_point = radius * orig_ball_point;
		bool found_overlap = false;
		for (const pdb_residue &the_res : prm_pdb) {
			for (const pdb_atom &the_atom : the_res) {
				if ( prm_pdb_atom .get_coord() != the_atom.get_coord() ) {
					found_overlap = found_overlap || access_overlap( prm_pdb_atom, ball_point, the_atom );
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
double cath::sec::get_accessibility_fraction(const pdb_atom &prm_pdb_atom, ///< The PDB atom for which the accessibility should be calculated
                                             const pdb      &prm_pdb,      ///< The PDB in which the accesibility should be calculated
                                             const size_t   &prm_number    ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                             ) {
	const coord_vec ball_points = make_dssp_ball_points( prm_number );

	return debug_numeric_cast<double>( get_accessibility_count( prm_pdb_atom, prm_pdb, prm_number ) )
	     / debug_numeric_cast<double>( make_dssp_ball_points( prm_number ).size() );
}

/// \brief Calculate the accessibility surface-area for the specified atom in the specified PDB
double cath::sec::get_accessibility_surface_area(const pdb_atom &prm_pdb_atom, ///< The PDB atom for which the accessibility should be calculated
                                                 const pdb      &prm_pdb,      ///< The PDB in which the accesibility should be calculated
                                                 const size_t   &prm_number    ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                                 ) {
	const double radius              = get_dssp_access_radius_with_water( prm_pdb_atom );
	const double sphere_surface_area = 4.0 * pi<double>() * radius * radius;

	return get_accessibility_fraction( prm_pdb_atom, prm_pdb, prm_number ) * sphere_surface_area;
}

/// \brief Calculate the accessibility surface-area for the specified residue in the specified PDB
double cath::sec::get_accessibility_surface_area(const pdb_residue &prm_pdb_residue, ///< The PDB residue for which the accessibility should be calculated
                                                 const pdb         &prm_pdb,         ///< The PDB in which the accesibility should be calculated
                                                 const size_t      &prm_number       ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                                 ) {
	return accumulate_proj(
		prm_pdb_residue,
		0.0,
		plus<>{},
		[&] (const pdb_atom &x) { return get_accessibility_surface_area( x, prm_pdb, prm_number ); }
	);
}

/// \brief Calculate the per-residue surface-area accessibilities for the specified PDBs
doub_vec cath::sec::calc_accessibilities(const pdb    &prm_pdb,   ///< The PDB in which the accesibility should be calculated
                                         const size_t &prm_number ///< The number to use to specify the sphere of points (tip: you should probably just use the default value)
                                         ) {
	return transform_build<doub_vec>(
		prm_pdb,
		[&] (const pdb_residue &x) { return get_accessibility_surface_area( x, prm_pdb, prm_number ); }
	);
}

// /// \brief TODOCUMENT
// enum class access_point_result : bool {
// 	VACANT, ///< TODOCUMENT
// 	OVERLAP ///< TODOCUMENT
// };

/// \brief Make a vector of all simple_locn_index values corresponding to the specified PDB
static vector<simple_locn_index> make_all_atom_entries(const pdb &prm_pdb ///< The PDB to query
                                                       ) {
	vector<simple_locn_index> results;
	for (const pdb_residue &the_residue : prm_pdb) {
		for (const pdb_atom &the_atom : the_residue) {
			results.push_back(
				make_simple_locn_index(
					the_atom.get_coord(),
					static_cast<unsigned int>( get_coarse_element_type( the_atom ) )
				)
			);
		}
	}
	return results;
}

/// \brief Build a lattice of the specified PDB's atoms using the specified cell size and maximum distance
template <sod Sod>
static auto make_access_atom_lattice(const pdb   &prm_pdb,       ///< The PDB for which the atom should be indexed
                                     const float &prm_cell_size, ///< The cell size of the lattice
                                     const float &prm_max_dist   ///< The maximum distance between points that lattic should be used to find
                                     ) {
	const auto keyer = make_res_pair_keyer(
		simple_locn_x_keyer_part{ prm_cell_size },
		simple_locn_y_keyer_part{ prm_cell_size },
		simple_locn_z_keyer_part{ prm_cell_size }
	);

	const auto all_atom_entries = make_all_atom_entries( prm_pdb );

	using cell_type  = vector< simple_locn_index >;
	using store_type = scan_index_lattice_store<decltype( keyer )::key_index_tuple_type, cell_type>;
	auto the_store = empty_store_maker<sod::SPARSE, store_type>{}( all_atom_entries, keyer, simple_locn_crit{ prm_max_dist * prm_max_dist } );
	for (const auto &data : all_atom_entries) {
		the_store.push_back_entry_to_cell(
			keyer.make_key( data ),
			data
		);
	}
	return the_store;
}

/// \brief Calculate the accessibilities using scanning
doub_vec cath::sec::calc_accessibilities_with_scanning(const pdb &prm_pdb ///< The PDB to query
                                                       ) {
	constexpr float MAX_DIST  = static_cast<float>( dssp_ball_constants::MAX_ATOM_DIST ); // 6.54
	constexpr float CELL_SIZE = 13.0;

	if ( prm_pdb.empty() ) {
		return {};
	}

	const auto the_store = make_access_atom_lattice<sod::SPARSE>( prm_pdb, CELL_SIZE, MAX_DIST );

	const auto keyer = make_res_pair_keyer(
		simple_locn_x_keyer_part{ CELL_SIZE },
		simple_locn_y_keyer_part{ CELL_SIZE },
		simple_locn_z_keyer_part{ CELL_SIZE }
	);

	const coord_vec dssp_ball_points = make_dssp_ball_points();
	using coord_opt = boost::optional<coord>;
	using coord_opt_vec = std::vector<coord_opt>;
	coord_opt_vec ball_points( dssp_ball_points.size() );

	return transform_build<doub_vec>(
		prm_pdb,
		[&] (const pdb_residue &this_residue) {
			return accumulate_proj(
				this_residue,
				0.0,
				plus<>{},
				[&] (const pdb_atom &this_atom) {
					const coord  &this_coord          = this_atom.get_coord();
					const double  this_radius         = get_dssp_access_radius_with_water( this_atom );
					const double  sphere_surface_area = 4.0 * pi<double>() * this_radius * this_radius;

					transform(
						dssp_ball_points,
						std::begin( ball_points ),
						[&] (const coord &x) { return make_optional( this_coord + ( this_radius * x ) ); }
					);

					const auto data = make_simple_locn_index(
						this_atom.get_coord(),
						static_cast<unsigned int>( get_coarse_element_type( this_atom ) )
					);

					for (const auto &key : common::cross( keyer.make_close_keys( data, simple_locn_crit{ MAX_DIST * MAX_DIST } ) ) ) {
						if ( the_store.has_matches( key ) ) {
							for (const simple_locn_index &eg : the_store.find_matches( key ) ) {
								if ( data == eg ) {
									continue;
								}
								const double  that_radius     = get_dssp_access_radius_with_water( static_cast<coarse_element_type>( eg.index ) );
								const double  that_radius_sq  = that_radius * that_radius;
								const double  total_radius    = this_radius + that_radius;
								const double  total_radius_sq = total_radius * total_radius;
								if ( are_within_distance_doub( data, eg, total_radius, total_radius_sq ) ) {
									const coord that_point = get_coord( eg );
									for (auto &ball_point : ball_points) {
										if ( ball_point ) {
											if ( squared_distance_between_points( *ball_point, that_point ) < that_radius_sq ) {
												ball_point = none;
											}
										}
									}
								}
							}
						}
					}

					return
						  sphere_surface_area
						* ( numeric_cast<double>( ball_points.size() ) - numeric_cast<double>( count( ball_points, none ) ) )
						/ numeric_cast<double>( ball_points.size() );
				}
			);
		}
	);
}