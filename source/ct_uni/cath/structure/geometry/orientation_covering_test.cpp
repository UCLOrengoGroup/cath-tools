/// \file
/// \brief The orientation_covering test suite

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

#include <boost/core/ignore_unused.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/boost_addenda/range/adaptor/limited.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/structure/geometry/orientation_covering.hpp"
#include "cath/structure/geometry/quat_rot.hpp"
#include "cath/test/global_test_constants.hpp"

#include <random>

using namespace ::cath::common;
using namespace ::cath::common::literals;
using namespace ::cath::geom;
using namespace ::std;

namespace cath {
	namespace test {

		/// \brief The orientation_covering_test_suite_fixture to assist in testing orientation_covering
		struct orientation_covering_test_suite_fixture : protected global_test_constants{
		protected:
			~orientation_covering_test_suite_fixture() noexcept = default;

			/// \brief TODOCUMENT
			template <typename T>
			void check_closest_neighbours_covers(const orientation_covering_impl<T> &prm_covering,      ///< TODOCUMENT
			                                     const quat_rot_impl<T>             &prm_orientation_a, ///< TODOCUMENT
			                                     const quat_rot_impl<T>             &prm_orientation_b  ///< TODOCUMENT
			                                     ) {
				const auto raw_angle   = angle_between_quat_rots( prm_orientation_a, prm_orientation_b );
				const auto roomy_angle = raw_angle * static_cast<T>( 1.001 );

				const auto general_neighbours   = calc_neighbours       ( prm_covering, roomy_angle );
				const auto closest_neighbour_a  = get_closest_neighbour ( prm_covering, prm_orientation_a );
				const auto closest_neighbour_b  = get_closest_neighbour ( prm_covering, prm_orientation_b );
				const auto candidates_for_a     = general_neighbours[ closest_neighbour_a ];
				const auto candidates_for_b     = general_neighbours[ closest_neighbour_b ];
				const auto closest_neighbours_a = get_closest_neighbours( prm_covering, prm_orientation_a, general_neighbours, roomy_angle );
				const auto closest_neighbours_b = get_closest_neighbours( prm_covering, prm_orientation_b, general_neighbours, roomy_angle );

//				if ( raw_angle > make_angle_from_degrees<T>( 21.5 ) && raw_angle < make_angle_from_degrees<T>( 23.5 ) ) {
//					const auto covering_size     = numeric_cast<double>( prm_covering.size() );
//					const auto neighbours_a_size = numeric_cast<double>( closest_neighbours_a.size() );
//					const auto neighbours_b_size = numeric_cast<double>( closest_neighbours_b.size() );
//					std::cerr << "***** The roomy angle is " << roomy_angle;
//					std::cerr << "\t closest neighbours a is " << closest_neighbours_a.size()
//					          << " (" << ( 100.0 * neighbours_a_size / covering_size ) << "%)";
//					std::cerr << "\tclosest neighbours b is "  << closest_neighbours_b.size()
//					          << " (" << ( 100.0 * neighbours_b_size / covering_size ) << "%)\n";
//				}

				const bool b_neighbours_contains_a = common::contains( closest_neighbours_b, closest_neighbour_a );
				if ( ! b_neighbours_contains_a ) {
					cerr << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n";
					std::cerr << "The roomy angle is " << roomy_angle << "\n";
					cerr << "prm_orientation_a    is    : " << prm_orientation_a   << "\n";
					cerr << "prm_orientation_b    is    : " << prm_orientation_b   << "\n";
					cerr << "closest_neighbour_a  is    : " << prm_covering[ closest_neighbour_a ] << "(" << angle_between_quat_rots( prm_orientation_a, prm_covering[ closest_neighbour_a ] ) << ")\n";
					cerr << "closest_neighbour_b  is    : " << prm_covering[ closest_neighbour_b ] << "(" << angle_between_quat_rots( prm_orientation_b, prm_covering[ closest_neighbour_b ] ) << ")\n";
					cerr << "angle twixt closest nghbrs : " << angle_between_quat_rots( prm_covering[ closest_neighbour_a ] , prm_covering[ closest_neighbour_b ] ) << "\n";
					for (const auto &x : closest_neighbours_b) {
						cerr << "closest_neighbours_b are.. : " << prm_covering[ x ] << "\n";
					}
					for (const auto &x : candidates_for_b) {
						cerr << "candidates for b are    .. : " << prm_covering[ x ] << "\n";
					}
					cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
					BOOST_CHECK( contains( closest_neighbours_b, closest_neighbour_a ) );
				}

				const bool a_neighbours_contains_b = common::contains( closest_neighbours_a, closest_neighbour_b );
				if ( ! a_neighbours_contains_b ) {
					cerr << "VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n";
					std::cerr << "The roomy angle is " << roomy_angle << "\n";
					cerr << "prm_orientation_a    is    : " << prm_orientation_a   << "\n";
					cerr << "prm_orientation_b    is    : " << prm_orientation_b   << "\n";
					cerr << "closest_neighbour_a  is    : " << prm_covering[ closest_neighbour_a ] << "(" << angle_between_quat_rots( prm_orientation_a, prm_covering[ closest_neighbour_a ] ) << ")\n";
					cerr << "closest_neighbour_b  is    : " << prm_covering[ closest_neighbour_b ] << "(" << angle_between_quat_rots( prm_orientation_b, prm_covering[ closest_neighbour_b ] ) << ")\n";
					cerr << "angle twixt closest nghbrs : " << angle_between_quat_rots( prm_covering[ closest_neighbour_a ] , prm_covering[ closest_neighbour_b ] ) << "\n";
					for (const auto &x : closest_neighbours_a) {
						cerr << "closest_neighbours_a are.. : " << prm_covering[ x ] << "\n";
					}
					for (const auto &x : candidates_for_a) {
						cerr << "candidates for a are    .. : " << prm_covering[ x ] << "\n";
					}
					cerr << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
					BOOST_CHECK( common::contains( closest_neighbours_a, closest_neighbour_b ) );
					// auto get_closest_neighbours(const orientation_covering_impl<T> &prm_orientations,   ///< TODOCUMENT
				}
			}
		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(orientation_covering_test_suite, cath::test::orientation_covering_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(mid_point_is_halfway, quat_rot_type, all_quat_rot_types) {
	const orientation_covering_impl<quat_rot_type> the_covering;
	for (const auto &x : cross( the_covering, the_covering ) | limited( 1000_z ) ) {
		const auto &orientation_a = get<0>( x );
		const auto &orientation_b = get<1>( x );
		const auto the_mid_point   = mid_point              ( orientation_a, orientation_b );
		const auto abs_of_midpoint = abs( the_mid_point );
		const auto angle_a_to_b    = angle_in_degrees( angle_between_quat_rots( orientation_a, orientation_b ) );
		const auto angle_mid_to_a  = angle_in_degrees( angle_between_quat_rots( orientation_a, the_mid_point ) );
		const auto angle_mid_to_b  = angle_in_degrees( angle_between_quat_rots( orientation_b, the_mid_point ) );
		BOOST_CHECK_CLOSE( abs_of_midpoint, 1.0,                             LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
		BOOST_CHECK_CLOSE( angle_mid_to_a,  angle_mid_to_b,                  LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
		BOOST_CHECK_CLOSE( angle_a_to_b,    angle_mid_to_a + angle_mid_to_b, LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
//		cerr << "a : "                  << orientation_a
//		     << "; mp : "               << the_mid_point
//		     << "; b : "                << orientation_b
//		     << "; mp : "               << the_mid_point
//		     << "; abs(mp) : "          << abs_of_midpoint
//		     << "; angle_a_to_b : "     << angle_a_to_b
//		     << "; angle_mid_to_a : "   << angle_mid_to_a
//		     << "; angle_mid_to_b : "   << angle_mid_to_b
//		     << "; norm_of_midpoint : " << abs( the_mid_point )
//		     << "\n";
	}
}

// Failure case:
// prm_orientation_a    is    : quat_rot[-0.201507,-0.609419,-0.744985,0.18166]
// prm_orientation_b    is    : quat_rot[0.0380258,-0.590544,-0.645591,0.482726]
// closest_neighbour_a  is    : quat_rot[      0,0.707107,0.707107,      0](33.445409°)
// closest_neighbour_b  is    : quat_rot[    0.5,   -0.5,   -0.5,    0.5](57.089721°)
// angle twixt closest nghbrs : 90.000000°
//
// - closest_neighbour_b is in candidates for a but isn't in closest_neighbours_a

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE_TEMPLATE(gubbins, quat_rot_type, all_quat_rot_types) {
//////	// ERROR SET 1 for c48u1.quat
////	const auto radius = make_angle_from_degrees<quat_rot_type>( 20.571917 );
////	const quat_rot_impl<quat_rot_type> A {  1.000000f,  0.000000f, 0.000000f,  0.000000f };
////	const quat_rot_impl<quat_rot_type> B {  0.707107f,  0.000000f, 0.707107f,  0.000000f };
////	const quat_rot_impl<quat_rot_type> a {  0.930342f, -0.198706f, 0.197123f, -0.236900f };
////	const quat_rot_impl<quat_rot_type> b {  0.875164f, -0.230578f, 0.363707f, -0.220543f };
//
////	// ERROR SET 2 for c48u1.quat
////	const auto radius = make_angle_from_degrees<quat_rot_type>( 45.937460 );
////	const quat_rot_impl<quat_rot_type> A {  0.0000000f,  0.707107f,  0.707107f,  0.000000f };
////	const quat_rot_impl<quat_rot_type> B {  0.5000000f, -0.500000f, -0.500000f,  0.500000f };
////	const quat_rot_impl<quat_rot_type> a { -0.2015070f, -0.609419f, -0.744985f,  0.181660f };
////	const quat_rot_impl<quat_rot_type> b {  0.0380258f, -0.590544f, -0.645591f,  0.482726f };
//
//	// ERROR SET 3 for c600v.quat
//	const auto radius = make_angle_from_degrees<quat_rot_type>( 21.268267 );
//	const quat_rot_impl<quat_rot_type> A {  0.500000f,  0.809017f, 0.309017f,  0.000000f };
//	const quat_rot_impl<quat_rot_type> B {  0.809017f,  0.309017f, 0.500000f,  0.000000f };
//	const quat_rot_impl<quat_rot_type> a {  0.478242f,  0.699823f, 0.510871f, -0.143330f };
//	const quat_rot_impl<quat_rot_type> b {  0.576349f,  0.565181f, 0.581116f, -0.103421f };
//
//	const quat_rot_impl<quat_rot_type> p1               = from_first_toward_second_at_angle(  a,  B, radius );
//	const quat_rot_impl<quat_rot_type> A_to_B_by_r      = from_first_toward_second_at_angle(  A,  B, radius );
//	const quat_rot_impl<quat_rot_type> dirn_A_to_b_by_r = rotation_between_rotations( A, from_first_toward_second_at_angle(  A,  B, radius ) );
//	const quat_rot_impl<quat_rot_type> apply_to_A       = dirn_A_to_b_by_r * A;
//	const quat_rot_impl<quat_rot_type> p2               = dirn_A_to_b_by_r * a;
//
//	cerr << "#################################" << "\n";
//	cerr << "A                   is " << A                                 << "\n";
//	cerr << "B                   is " << B                                 << "\n";
//	cerr << "a                   is " << a                                 << "\n";
//	cerr << "b                   is " << b                                 << "\n";
//	cerr << "A_to_B_by_r         is " << A_to_B_by_r                       << "\n";
//	cerr << "dirn_A_to_b_by_r    is " << dirn_A_to_b_by_r                  << "\n";
//	cerr << "apply_to_A          is " << apply_to_A                        << "\n";
//	cerr << "p2                  is " << p2                                << "\n";
//	cerr << "#################################" << "\n";
//	cerr << "angle from  A to  B is " << angle_between_quat_rots(  A,  B ) << "\n";
//	cerr << "angle from  a to  b is " << angle_between_quat_rots(  a,  b ) << "\n";
//	cerr << "angle from  a to  A is " << angle_between_quat_rots(  a,  A ) << "\n";
//	cerr << "angle from  a to  B is " << angle_between_quat_rots(  a,  B ) << "\n";
//	cerr << "angle from  b to  A is " << angle_between_quat_rots(  b,  A ) << "\n";
//	cerr << "angle from  b to  B is " << angle_between_quat_rots(  b,  B ) << "\n";
//	cerr << "#################################" << "\n";
//	cerr << "angle from  a to p1 is " << angle_between_quat_rots(  a, p1 ) << "\n";
//	cerr << "angle from p1 to  A is " << angle_between_quat_rots( p1,  A ) << "\n";
//	cerr << "angle from p1 to  B is " << angle_between_quat_rots( p1,  B ) << "\n";
//	cerr << "#################################" << "\n";
//	cerr << "angle from  a to p2 is " << angle_between_quat_rots(  a, p2 ) << "\n";
//	cerr << "angle from p2 to  A is " << angle_between_quat_rots( p2,  A ) << "\n";
//	cerr << "angle from p2 to  B is " << angle_between_quat_rots( p2,  B ) << "\n";
//	cerr << "#################################" << "\n";
//}

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE_TEMPLATE(neighbours, quat_rot_type, all_quat_rot_types) {
//	const orientation_covering_impl<quat_rot_type> the_covering;
////	constexpr size_t num_repeats = 400;
////	constexpr size_t num_repeats = 2000;
//	constexpr size_t num_repeats = 50000;
//
//	default_random_engine the_rng{ random_device{}() };
////	default_random_engine the_rng{};
//
//	size_t repeat_ctr = 0;
//	while ( repeat_ctr < num_repeats ) {
////	for (const auto &repeat_ctr : indices( num_repeats ) ) {
////		ignore_unused( repeat_ctr );
//		const auto quat_rot_1 = make_random_quat_rot<quat_rot_type>( the_rng );
//		const auto quat_rot_2 = make_random_quat_rot<quat_rot_type>( the_rng );
//		const auto the_angle  = angle_between_quat_rots( quat_rot_1, quat_rot_2 );
//		if ( the_angle > make_angle_from_degrees<quat_rot_type>( 21.5 ) && the_angle < make_angle_from_degrees<quat_rot_type>( 23.5 ) ) {
//			if ( repeat_ctr % 100 == 0 ) {
//				cerr << "repeat_ctr : " << repeat_ctr << "\n";
//			}
////			cerr << quat_rot_1 << ", " << quat_rot_2 << "\n";
//			check_closest_neighbours_covers( the_covering, quat_rot_1, quat_rot_2 );
//			++repeat_ctr;
//		}
//	}
//
////	const auto the_neighbours = calc_neighbours( the_covering, make_angle_from_degrees<quat_rot_type>( 22.5 ) );
//
////size_vec_vec calc_neighbours(const orientation_covering_impl<T> &prm_orientations, ///< TODOCUMENT
////		                             const geom::angle<A>               &prm_search_radius ///< TODOCUMENT
////		                             ) {
//}

BOOST_AUTO_TEST_SUITE_END()
