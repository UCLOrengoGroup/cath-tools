/// \file
/// \brief The quat_rot test suite

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

#include <boost/test/unit_test.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include <gnuplot-iostream.h>

#include "common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "common/boost_addenda/range/adaptor/limited.hpp"
#include "common/size_t_literal.hpp"
#include "scan/detail/scan_type_aliases.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/quat_rot.hpp"

#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::geometry::cs::cartesian;
using boost::geometry::model::point;
using boost::mpl::contains;
using boost::numeric_cast;
using boost::range::for_each;

/// \todo There are problems with quaternion tests when using double (or long double).
///       For now, the tests are only run on floats because that's what's used in the
///       scan code but if the scan code starts using other types, they must be tested
///       here and problems must be fixed.
static_assert(
	contains< all_quat_rot_types, scan::detail::frame_quat_rot_type>::value,
	"quaternion tests must be applied to the type used in scan code"
);

namespace cath {
	namespace test {

		/// \brief The quat_rot_test_suite_fixture to assist in testing quat_rot
		struct quat_rot_test_suite_fixture : protected global_test_constants {
		protected:
			~quat_rot_test_suite_fixture() noexcept = default;

			rotation_vec make_all_rotations_between_coords() const;

			template <typename T>
			vector<quat_rot_impl<T> > speed_test() const;

			template <typename T>
			vector<quat_rot_impl<T> > make_all_quat_rots_between_coords() const;
			doub_doub_pair_vec get_num_degrees_and_distance_1s() const;

			const coord_list example_coords{ coord_vec{
				{  1.0000000000000,  0.0000000000000,  0.0000000000000 },
				{  0.0000000000000,  1.0000000000000,  0.0000000000000 },
				{  0.0000000000000,  0.0000000000000,  1.0000000000000 },
				{ 65.9858249434368, 63.9621735856966, 45.1548984813602 },
		//		{ 66.4112292070971, 67.5423802415938, 31.4274822901101 },
				{ 39.6225754478113, 96.5897020491870, 32.9884091233001 },
		//		{ 43.8852218601788, 80.6365120744690, 59.7468007831253 },
				{ 94.6450202323984, 82.4461683283836, 25.0011455817333 },
		//		{ 58.3341608029187, 24.4223077180205, 41.2949279615525 },
				{  9.9823139478932, 65.4564285568000,  9.6871266526090 },
		//		{ 77.3463733867068, 74.2178851637863, 27.8122348907729 },
				{ 87.8739995225285, 63.0594667719400, 11.7214853422045 },
		//		{ 68.2096774515401, 50.1443867910421, 13.0103777419585 },
				{ 87.3584607495264, 36.3240031030912, 25.2551886993555 },
		//		{ 59.9520610852828, 91.5203650873085, 25.2456351180143 },
				{ 24.3720382058324, 21.5176590983852, 84.5428300026544 },
		//		{ 88.1958578037469, 98.8738986208173, 23.1059680458280 },
				{ 84.8254493891918, 64.8183702718370, 60.8994132103174 },
		//		{ 99.8331707398826, 35.5956607016186, 53.4450486217580 },
				{ 37.0202959607465, 29.1249939081460, 77.1123204960364 },
		//		{ 24.4687681399256, 15.5748340689978, 72.9988060733383 },
				{ 94.9420279756755, 91.7826904592900, 74.7498247973137 },
		//		{ 76.8498353657513, 84.6880307171180, 52.0741712352276 },
				{ 53.6690233607356, 60.1760632757340,  0.9460322257655 },
		//		{ 75.0650480194352, 35.0499657463669, 88.7813460811394 },
				{ 57.8541927710493, 17.6036819051291, 23.2608820482493 },
		//		{ 71.2773368028202, 42.7593543317659, 28.4354064240787 },
				{ 23.9898632893496, 22.9157333000522, 48.4138086899325 },
				{ 72.6671393488573,  4.6865190429724, 73.7396687671755 }
			} };

		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
rotation_vec cath::test::quat_rot_test_suite_fixture::make_all_rotations_between_coords() const {
	rotation_vec rotations;
	for (const coord &coord_1 : example_coords) {
		for (const coord &coord_2 : example_coords) {
			if ( coord_1 != coord_2 ) {
				rotations.push_back( rotation_to_x_axis_and_x_y_plane( coord_1, coord_2 ) );
			}
		}
	}
	return rotations;
}

/// \brief TODOCUMENT
template <typename T>
vector<quat_rot_impl<T> > cath::test::quat_rot_test_suite_fixture::make_all_quat_rots_between_coords() const {
	const rotation_vec rotations = make_all_rotations_between_coords();
	vector<quat_rot_impl<T> > quat_rots;
	for (const rotation &the_rotn : rotations) {
		quat_rots.push_back( make_quat_rot_from_rotation<T>( the_rotn ) );
	}
	return quat_rots;
}

/// \brief TODOCUMENT
doub_doub_pair_vec cath::test::quat_rot_test_suite_fixture::get_num_degrees_and_distance_1s() const {
	return { {   0.0, 0.00000000000000000000000000000000 },
	         {  10.0, 0.00380530190825447729848982070244 },
	         {  20.0, 0.01519224698779192943157503870030 },
	         {  30.0, 0.03407417371093169748989468170740 },
	         {  40.0, 0.06030737921409161828197043053730 },
	         {  50.0, 0.09369221296335002668008440362970 },
	         {  60.0, 0.13397459621556132117496315525600 },
	         {  70.0, 0.18084795571100817981055955407900 },
	         {  80.0, 0.23395555688102194149971882475300 },
	         {  90.0, 0.29289321881345245393406251377400 },
	         { 100.0, 0.35721239031346066359002069945400 },
	         { 110.0, 0.42642356364895389511226325707200 },
	         { 120.0, 0.49999999999999992367216705702000 },
	         { 130.0, 0.57738173825930058549304318971100 },
	         { 140.0, 0.65797985667433116988966207427900 },
	         { 150.0, 0.74118095489747929641068807660200 },
	         { 160.0, 0.82635182233306971912034644134300 },
	         { 170.0, 0.91284425725234203530386242753100 } };
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(quat_rot_test_suite, cath::test::quat_rot_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(simple, quat_rot_type, all_quat_rot_types) {
	const rotation                     simple_rot_matx = rotation_to_x_axis_and_x_y_plane( coord( 0, 0, 1), coord( 0, 1, 0) );
	const quat_rot_impl<quat_rot_type> simple_quat_rot = make_quat_rot_from_rotation<quat_rot_type>( simple_rot_matx );
	BOOST_CHECK_EQUAL( simple_rot_matx, make_rotation_from_quat_rot( simple_quat_rot ) );
	BOOST_CHECK_CLOSE( angle_in_degrees( angle_of_rotation( simple_rot_matx ) ), angle_in_degrees( angle_of_quat_rot( simple_quat_rot ) ), LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(distance, quat_rot_type, all_quat_rot_types) {
	const doub_doub_pair_vec num_degrees_and_distance_1s = get_num_degrees_and_distance_1s();
	for (const doub_doub_pair &num_degrees_and_distance_1 : num_degrees_and_distance_1s) {
		const double &num_degrees = num_degrees_and_distance_1.first;
		const double &distance_1  = num_degrees_and_distance_1.second;
		BOOST_CHECK_CLOSE(
			boost::numeric_cast<double>( distance_1_of_angle<quat_rot_type>( make_angle_from_degrees<double>( num_degrees ) ) ),
			distance_1,
			LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>()
		);
	}
}
/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(unchanged_by_conversion_to_rotation_and_back, quat_rot_type, all_quat_rot_types) {
	using quat_rot_t = quat_rot_impl<quat_rot_type>;
	const auto quat_rot_1 = quat_rot_t{  1.0,  0.0,  0.0,  0.0 };
	const auto quat_rot_2 = quat_rot_t{  0.0,  1.0,  0.0,  0.0 };
	const auto quat_rot_3 = quat_rot_t{  0.0,  0.0,  1.0,  0.0 };
	const auto quat_rot_4 = quat_rot_t{  0.0,  0.0,  0.0,  1.0 };
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot(  quat_rot_1 ) ), quat_rot_1 );
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot(  quat_rot_2 ) ), quat_rot_2 );
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot(  quat_rot_3 ) ), quat_rot_3 );
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot(  quat_rot_4 ) ), quat_rot_4 );
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot( -quat_rot_1 ) ), quat_rot_1 );
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot( -quat_rot_2 ) ), quat_rot_2 );
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot( -quat_rot_3 ) ), quat_rot_3 );
	BOOST_CHECK_EQUAL( make_quat_rot_from_rotation<quat_rot_type>( make_rotation_from_quat_rot( -quat_rot_4 ) ), quat_rot_4 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(from_first_toward_second_at_angle_works, quat_rot_type, all_quat_rot_types) {
	const auto quat_rots = make_all_quat_rots_between_coords<quat_rot_type>();
	for (const auto &x : cross( quat_rots, quat_rots ) | limited( 1000_z ) ) {
		const auto &quat_rot_a   = get<0>( x );
		const auto &quat_rot_b   = get<1>( x );
		const auto angle_between = angle_between_quat_rots( quat_rot_a, quat_rot_b );
		const bool small_angle   = ( angle_between <= make_angle_from_degrees<quat_rot_type>( 2.0 ) );
		const auto accuracy      = ( small_angle ? 100.0 : 1.0 ) * LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>();
		for (const double &frac : { 0.1, 0.25, 0.5, 0.75, 0.9 } ) {
			const auto expected_angle = static_cast<quat_rot_type>( frac ) * angle_between;
			const auto new_quat_rot   = from_first_toward_second_at_angle( quat_rot_a, quat_rot_b, expected_angle );

			const auto got_angle      = angle_between_quat_rots( quat_rot_a, new_quat_rot );
			if ( angle_in_degrees( got_angle ) == 0.0 ) {
				BOOST_CHECK( abs( angle_in_degrees( expected_angle ) - angle_in_degrees( got_angle ) ) < accuracy );
			}
			else {
				BOOST_CHECK_CLOSE( angle_in_degrees( expected_angle ), angle_in_degrees( got_angle ), accuracy );
			}

			const auto expected_rev_angle = static_cast<quat_rot_type>( 1.0 - frac ) * angle_between;
			const auto got_rev_angle      = angle_between_quat_rots( new_quat_rot, quat_rot_b );
			if ( angle_in_degrees( got_rev_angle ) == 0.0 ) {
				BOOST_CHECK( abs( angle_in_degrees( expected_rev_angle ) - angle_in_degrees( got_rev_angle ) ) < accuracy );
			}
			else {
				BOOST_CHECK_CLOSE( angle_in_degrees( expected_rev_angle ), angle_in_degrees( got_rev_angle ), accuracy );
			}
		}

	}
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(rotation_between_quat_rots, quat_rot_type, all_quat_rot_types) {
	using point_type = point<quat_rot_type, 3, cartesian>;
	const auto rotations = make_all_rotations_between_coords();
	const auto quat_rots = make_all_quat_rots_between_coords<quat_rot_type>();
	for_each(
		cross( rotations, rotations ),
		cross( quat_rots, quat_rots ),
		[&] (const tuple<rotation, rotation> &x, const tuple<quat_rot_impl<quat_rot_type>, quat_rot_impl<quat_rot_type>> &y) {
			const auto rotn_answer = rotation_between_rotations( get<0>( x ), get<1>( x ) );
			const auto quat_answer = rotation_between_rotations( get<0>( y ), get<1>( y ) );
			const auto rotq_answer = make_quat_rot_from_rotation<quat_rot_type>( rotn_answer );
			const auto distance    = distance_1_between_quat_rots( rotq_answer, quat_answer);
			if ( distance >= LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() ) {
				cerr << "Source.rotn.first  : " << get<0>( x ) << "\n";
				cerr << "Source.rotn.second : " << get<1>( x ) << "\n";
				cerr << "Source.quat.first  : " << get<0>( y ) << " - length : " << abs( get<0>( y ) ) << "\n";
				cerr << "Source.quat.second : " << get<1>( y ) << " - length : " << abs( get<1>( y ) ) << "\n";
				cerr << "Rotn direct answer : " << rotn_answer << "\n";
				cerr << "Rotn answer        : " << rotq_answer << " - length : " << abs( rotq_answer ) << "\n";
				cerr << "Quat answer        : " << quat_answer << " - length : " << abs( quat_answer ) << "\n";

				cerr << "\n";

				cerr << "Rot eg first       : " <<                                              rotate_copy( get<0>( x ), UNIT_X_COORD )       << "\t"
				                                <<                                              rotate_copy( get<0>( x ), UNIT_Y_COORD )       << "\t"
				                                <<                                              rotate_copy( get<0>( x ), UNIT_Z_COORD )       << "\n";

				cerr << "Rot eg second      : " <<                                              rotate_copy( get<1>( x ), UNIT_X_COORD )       << "\t"
				                                <<                                              rotate_copy( get<1>( x ), UNIT_Y_COORD )       << "\t"
				                                <<                                              rotate_copy( get<1>( x ), UNIT_Z_COORD )       << "\n";

				cerr << "Rot answer direct  : " <<        rotate_copy( rotn_answer,             rotate_copy( get<0>( x ), UNIT_X_COORD )   )   << "\t"
				                                <<        rotate_copy( rotn_answer,             rotate_copy( get<0>( x ), UNIT_Y_COORD )   )   << "\t"
				                                <<        rotate_copy( rotn_answer,             rotate_copy( get<0>( x ), UNIT_Z_COORD )   )   << "\n";

				cerr << "Rotn answer        : " << coord{ rotate_copy( rotq_answer, point_type( rotate_copy( get<0>( x ), UNIT_X_COORD ) ) ) } << "\t"
				                                << coord{ rotate_copy( rotq_answer, point_type( rotate_copy( get<0>( x ), UNIT_Y_COORD ) ) ) } << "\t"
				                                << coord{ rotate_copy( rotq_answer, point_type( rotate_copy( get<0>( x ), UNIT_Z_COORD ) ) ) } << "\n";

				cerr << "\n";

				cerr << "Rot eg first       : " <<                                              rotate_copy( get<0>( x ), UNIT_X_COORD )       << "\t"
				                                <<                                              rotate_copy( get<0>( x ), UNIT_Y_COORD )       << "\t"
				                                <<                                              rotate_copy( get<0>( x ), UNIT_Z_COORD )       << "\n";

				cerr << "Rot eg second      : " <<                                              rotate_copy( get<1>( x ), UNIT_X_COORD )       << "\t"
				                                <<                                              rotate_copy( get<1>( x ), UNIT_Y_COORD )       << "\t"
				                                <<                                              rotate_copy( get<1>( x ), UNIT_Z_COORD )       << "\n";

				cerr << "Quat answer        : " << coord{ rotate_copy( quat_answer, point_type( rotate_copy( get<0>( x ), UNIT_X_COORD ) ) ) } << "\t"
				                                << coord{ rotate_copy( quat_answer, point_type( rotate_copy( get<0>( x ), UNIT_Y_COORD ) ) ) } << "\t"
				                                << coord{ rotate_copy( quat_answer, point_type( rotate_copy( get<0>( x ), UNIT_Z_COORD ) ) ) } << "\n";

				cerr << "Distance           : " << distance << "\n";
				cerr << "\n";
			}
			BOOST_CHECK_LT( distance, LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
		}
	);
}

// example_coords.size() : 17
// rotn_a           : rotation[      1,       0,       0;      -0, 0.94633, 0.323202;       0, -0.323202, 0.94633]
// rotn_b           : rotation[0.361888, 0.882189, 0.301295; 0.932221, -0.342466, -0.116963; -7.45058e-08, 0.323202, -0.94633]
// quat_a           : quat_rot[0.986491,-0.163814,      0,     -0]
// quat_b           : quat_rot[0.135178,0.814046,0.55722,0.0925302]
// quat_quat_answer : quat_rot[-2.5332e-07,0.825193,0.56485,-2.23517e-07]
// rotn_quat_answer : quat_rot[0.000272957,-0.825193,-0.56485,-0.00012207]
// eg_coord         : coord[  87.3585,    36.324,   25.2552]
// target           : coord[  71.2679,   66.0438,  -12.1598]
// answer_quat      : coord[  71.2679,   66.0438,  -12.1598]
// answer_rotn      : coord[  71.2666,   66.0509,  -12.1285]
// wrong rotn answer         : rotation[0.361888, 0.932222, -5.04558e-07; 0.932221, -0.361888, 9.67772e-08; -7.45058e-08, -4.98712e-07,      -1] **** ROTN_SID
// should-a-been rotn answer : rotation[0.361888, 0.932221, -6.55066e-07; 0.932221, -0.361888, 1.65568e-07; -8.27146e-08, -6.70583e-07,      -1]

// Source.rotn.first  : rotation[      1,       0,       0;      -0, 0.94633, 0.323201;       0, -0.323201, 0.94633]
// Source.rotn.second : rotation[0.361888, 0.88219, 0.301295; 0.932222, -0.342465, -0.116963;       0, 0.323201, -0.94633]
// Source.quat.first  : quat_rot[0.986491,-0.163814,      0,     -0]
// Source.quat.second : quat_rot[0.135178,0.814046,0.55722,0.0925302]
// Rotn direct answer : rotation[0.361888, 0.932222, -5.55112e-17; 0.932222, -0.361888,       0;       0,       0,      -1] **** ROTN_BOB
// Rotn answer        : quat_rot[      0,0.825193,-0.564851,      0]
// Quat answer        : quat_rot[-2.98023e-08,0.825193,0.564851,-3.72529e-08]
// Distance           : 0.638112

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE_TEMPLATE(rotation_between_quat_rots_example_2, quat_rot_type, all_quat_rot_types) {
//
//	cerr << "example_coords.size() : " << example_coords.size() << endl;
//
//	const auto quat_a   = quat_rot_impl<quat_rot_type>(     0.986491f, -0.163814f,       0.0f,         -0.0f );
//	const auto quat_b   = quat_rot_impl<quat_rot_type>(     0.135178f,  0.814046f,   0.55722f,    0.0925302f );
//	const auto rotn_a   = make_rotation_from_quat_rot( quat_a );
//	const auto rotn_b   = make_rotation_from_quat_rot( quat_b );
//
//	const auto rotn_answer      = rotation_between_rotations( rotn_a, rotn_b );
//	const auto quat_quat_answer = rotation_between_rotations( quat_a, quat_b );
//	const auto rotn_quat_answer = make_quat_rot_from_rotation<quat_rot_type>( rotn_answer );
//
//	// const auto guess_a  = quat_rot_impl<quat_rot_type>(          0.0f,  0.825193f, -0.564851f,          0.0f );
//	// const auto guess_b  = quat_rot_impl<quat_rot_type>( -2.98023e-08f,  0.825193f,  0.564851f, -3.72529e-08f );
//	const auto eg_coord = point<quat_rot_type, 3, cartesian>{ 87.3584607495264f, 36.3240031030912f, 25.2551886993555f };
//	const auto target   = rotate_copy( quat_b, eg_coord );
//	const auto answer_quat = rotate_copy( quat_quat_answer, rotate_copy( quat_a, eg_coord ) );
//	const auto answer_rotn = rotate_copy( rotn_quat_answer, rotate_copy( quat_a, eg_coord ) );
//
//	// Rotn answer   : quat_rot[          0.0f, 0.825193f, -0.564851f,          0.0f ]
//	// Quat answer   : quat_rot[ -2.98023e-08f, 0.825193f,  0.564851f, -3.72529e-08f ]
//
//	cerr << "ROTN_BOB         : " << make_quat_rot_from_rotation<quat_rot_type>( rotation{ 0.361888, 0.932222, -5.55112e-17, 0.932222, -0.361888,           0,            0,            0,      -1 } ) << "\n";
//	cerr << "ROTN_SID         : " << make_quat_rot_from_rotation<quat_rot_type>( rotation{ 0.361888, 0.932222, -5.04558e-07, 0.932221, -0.361888, 9.67772e-08, -7.45058e-08, -4.98712e-07,      -1 } ) << "\n";
//
//	cerr << "rotn_a           : " << rotn_a           << "\n";
//	cerr << "rotn_b           : " << rotn_b           << "\n";
//	cerr << "quat_a           : " << quat_a           << "\n";
//	cerr << "quat_b           : " << quat_b           << "\n";
//	cerr << "quat_quat_answer : " << quat_quat_answer << "\n";
//	cerr << "rotn_quat_answer : " << rotn_quat_answer << "\n";
//	cerr << "eg_coord         : " << eg_coord         << "\n";
//	cerr << "target           : " << target           << "\n";
//	cerr << "answer_quat      : " << answer_quat      << "\n";
//	cerr << "answer_rotn      : " << answer_rotn      << "\n";
//	cerr << "\n";
//
//	cerr << "wrong rotn answer         : " << rotn_answer << "\n";
//	cerr << "should-a-been rotn answer : " << make_rotation_from_quat_rot( quat_quat_answer ) << "\n";
//}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(rotation_between_quat_rots_example_1, quat_rot_type, all_quat_rot_types) {
	const auto oort     = static_cast  <quat_rot_type>( 1.0 / sqrt( 2.0 ) );
	const auto quat_a   = quat_rot_impl<quat_rot_type>(  0.5,  0.5,  0.5,  0.5 );
	const auto quat_b   = quat_rot_impl<quat_rot_type>(  0.0, oort, oort,  0.0 );
	const auto expected = quat_rot_impl<quat_rot_type>( oort,  0.0, oort,  0.0 );
	const auto got      = rotation_between_rotations  ( quat_a, quat_b );
	const auto distance = distance_1_between_quat_rots( expected, got  );
	if ( distance >= ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() ) {
		cerr << "Quat_a   : " << quat_a   << "\n";
		cerr << "Quat_b   : " << quat_b   << "\n";
		cerr << "Expected : " << expected << "\n";
		cerr << "Got      : " << got      << "\n";
		cerr << "Distance : " << distance << "\n";
		cerr << "\n";
	}
	BOOST_CHECK_LT( distance, ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(rotate_point_by_quat_rots, quat_rot_type, all_quat_rot_types) {
	const auto rotations = make_all_rotations_between_coords();
	const auto quat_rots = make_all_quat_rots_between_coords<quat_rot_type>();
	for_each(
		rotations,
		quat_rots,
		[&] (const rotation &x, const quat_rot_impl<quat_rot_type> &y) {
			for (const coord &the_coord : example_coords) {
				const auto bg_coord = point<quat_rot_type, 3, cartesian>( the_coord );
				const auto rotn_answer = rotate_copy( x, the_coord );
				const auto quat_answer = rotate_copy( y, bg_coord  );
				const auto distance    = distance_between_points( rotn_answer, coord{ quat_answer } );
				if ( distance >= LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() ) {
					cerr << rotn_answer << "\n";
					cerr << coord{ quat_answer } << "\n";
					cerr << distance_between_points( rotn_answer, coord{ quat_answer } ) << "\n";
					cerr << "\n";
				}
				BOOST_CHECK_LT( distance, LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
			}
		}
	);
}

//interpolate

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(speed_test) {
//	BOOST_CHECK( true );
//	const auto        cutoff_angle = make_angle_from_degrees<double>( 67.5 );
//	const float       dist_1_float = distance_1_of_angle< float       >( cutoff_angle );
//	const double      dist_1_doubl = distance_1_of_angle< double      >( cutoff_angle );
//	const long double dist_1_longd = distance_1_of_angle< long double >( cutoff_angle );
//
//	const rotation_vec                          all_rots       = make_all_rotations_between_coords();
//	const vector<quat_rot_impl< float       > > all_qrs_float  = make_all_quat_rots_between_coords< float       >();
//	const vector<quat_rot_impl< double      > > all_qrs_doubl  = make_all_quat_rots_between_coords< double      >();
//	const vector<quat_rot_impl< long double > > all_qrs_longd  = make_all_quat_rots_between_coords< long double >();
////	const size_t                                num_rots       = all_rots.size();
////	const size_t                                num_pairs      = num_rots * num_rots;
////	const double                                num_pairs_doub = numeric_cast<double>( num_pairs );
//
////	cerr << endl;
////	cerr << num_rots << " * " << num_rots << " = " << num_pairs << endl;
//
////	const auto rotn_begin_time = high_resolution_clock::now();
//	size_t      rotn_count = 0;
//	for (const rotation &rot_a : all_rots) {
//		for (const rotation &rot_b : all_rots) {
//			if ( angle_between_rotations( rot_a, rot_b ) < cutoff_angle ) {
//				++rotn_count;
//			}
//		}
//	}
////	const duration rotn_duration = high_resolution_clock::now() - rotn_begin_time;
////	cerr << fixed;
////	cerr << right;
////	cerr << "rotation    : " << rotn_count << "\t" << num_pairs_doub * 1000000000.0 / numeric_cast<double>( rotn_duration.total_nanoseconds() ) << " / second" << endl;
//
//
////	const auto qr_float_angle_begin_time = high_resolution_clock::now();
//	size_t      qr_float_angle_count = 0;
//	for (const quat_rot_impl<float> &rot_a : all_qrs_float) {
//		for (const quat_rot_impl<float> &rot_b : all_qrs_float) {
//			if ( angle_between_quat_rots( rot_a, rot_b ) < cutoff_angle ) {
//				++qr_float_angle_count;
//			}
//		}
//	}
////	const duration qr_float_angle_duration = high_resolution_clock::now() - qr_float_angle_begin_time;
////	cerr << "float/angle : " << qr_float_angle_count << "\t" << num_pairs_doub * 1000000000.0 / numeric_cast<double>( qr_float_angle_duration.total_nanoseconds() ) << " / second" << endl;
//
//
////	const auto qr_float_inner_begin_time = high_resolution_clock::now();
//	size_t      qr_float_inner_count = 0;
//	for (const quat_rot_impl<float> &rot_a : all_qrs_float) {
//		for (const quat_rot_impl<float> &rot_b : all_qrs_float) {
//			if ( distance_1_between_quat_rots( rot_a, rot_b ) < dist_1_float ) {
//				++qr_float_inner_count;
//			}
//		}
//	}
////	const duration qr_float_inner_duration = high_resolution_clock::now() - qr_float_inner_begin_time;
////	cerr << "float/inner : " << qr_float_inner_count << "\t" << num_pairs_doub * 1000000000.0 / numeric_cast<double>( qr_float_inner_duration.total_nanoseconds() ) << " / second" << endl;
//
////	const auto qr_doubl_angle_begin_time = high_resolution_clock::now();
//	size_t      qr_doubl_angle_count = 0;
//	for (const quat_rot_impl<double> &rot_a : all_qrs_doubl) {
//		for (const quat_rot_impl<double> &rot_b : all_qrs_doubl) {
//			if ( angle_between_quat_rots( rot_a, rot_b ) < cutoff_angle ) {
//				++qr_doubl_angle_count;
//			}
//		}
//	}
////	const duration qr_doubl_angle_duration = high_resolution_clock::now() - qr_doubl_angle_begin_time;
////	cerr << "doubl/angle : " << qr_doubl_angle_count << "\t" << num_pairs_doub * 1000000000.0 / numeric_cast<double>( qr_doubl_angle_duration.total_nanoseconds() ) << " / second" << endl;
//
//
////	const auto qr_doubl_inner_begin_time = high_resolution_clock::now();
//	size_t      qr_doubl_inner_count = 0;
//	for (const quat_rot_impl<double> &rot_a : all_qrs_doubl) {
//		for (const quat_rot_impl<double> &rot_b : all_qrs_doubl) {
//			if ( distance_1_between_quat_rots( rot_a, rot_b ) < dist_1_doubl ) {
//				++qr_doubl_inner_count;
//			}
//		}
//	}
////	const duration qr_doubl_inner_duration = high_resolution_clock::now() - qr_doubl_inner_begin_time;
////	cerr << "doubl/inner : " << qr_doubl_inner_count << "\t" << num_pairs_doub * 1000000000.0 / numeric_cast<double>( qr_doubl_inner_duration.total_nanoseconds() ) << " / second" << endl;
//
////	const auto qr_lngd_angle_begin_time = high_resolution_clock::now();
//	size_t      qr_lngd_angle_count = 0;
//	for (const quat_rot_impl<long double> &rot_a : all_qrs_longd) {
//		for (const quat_rot_impl<long double> &rot_b : all_qrs_longd) {
//			if ( angle_between_quat_rots( rot_a, rot_b ) < cutoff_angle ) {
//				++qr_lngd_angle_count;
//			}
//		}
//	}
////	const duration qr_longd_angle_duration = high_resolution_clock::now() - qr_lngd_angle_begin_time;
////	cerr << "longd/angle : " << qr_lngd_angle_count << "\t" << num_pairs_doub * 1000000000.0 / numeric_cast<double>( qr_longd_angle_duration.total_nanoseconds() ) << " / second" << endl;
//
//
////	const auto qr_lngd_inner_begin_time = high_resolution_clock::now();
//	size_t      qr_lngd_inner_count = 0;
//	for (const quat_rot_impl<long double> &rot_a : all_qrs_longd) {
//		for (const quat_rot_impl<long double> &rot_b : all_qrs_longd) {
//			if ( distance_1_between_quat_rots( rot_a, rot_b ) < dist_1_longd ) {
//				++qr_lngd_inner_count;
//			}
//		}
//	}
////	const duration qr_longd_inner_duration = high_resolution_clock::now() - qr_lngd_inner_begin_time;
////	cerr << "longd/inner : " << qr_lngd_inner_count << "\t" << num_pairs_doub * 1000000000.0 / numeric_cast<double>( qr_longd_inner_duration.total_nanoseconds() ) << " / second" << endl;
//}

/// \todo Make these tests work
///
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE_TEMPLATE(conversions_work_and_angles_match, quat_rot_type, all_quat_rot_types) {
//	const rotation_vec                          rotns = make_all_rotations_between_coords();
//	for (const rotation &the_rot_matx : rotns) {
//		const quat_rot_impl<quat_rot_type> the_quat_rot = make_quat_rot_from_rotation<quat_rot_type>( the_rot_matx );
//		BOOST_CHECK_EQUAL( the_rot_matx, make_rotation_from_quat_rot( the_quat_rot ) );
//		BOOST_CHECK_CLOSE( angle_in_degrees( angle_of_rotation( the_rot_matx ) ), angle_in_degrees( angle_of_quat_rot( the_quat_rot ) ), ACCURACY_PERCENTAGE() );
//	}
//}

/// \todo Make these tests work
///
///// \brief TODOCUMENT
/////
///// \todo Should also test that both produce the same answer when applied as
/////       a rotation to each of the coords
//BOOST_AUTO_TEST_CASE_TEMPLATE(angle_between_rotations_works, quat_rot_type, all_quat_rot_types) {
//	const rotation_vec                          rotns = make_all_rotations_between_coords();
//    const vector<quat_rot_impl<quat_rot_type> > quats = make_all_quat_rots_between_coords<quat_rot_type>();
//
//    BOOST_REQUIRE_EQUAL( rotns.size(), quats.size() );
//
//    for (const size_t &ctr_a : indices( rotns.size() ) ) {
//    	for (const size_t &ctr_b : indices( rotns.size() ) ) {
//    		const angle rotn_angle_a_to_b = angle_between_rotations( rotns[ ctr_a ], rotns[ ctr_b ] );
//    		const angle rotn_angle_b_to_a = angle_between_rotations( rotns[ ctr_b ], rotns[ ctr_a ] );
//    		const angle quat_angle_a_to_b = angle_between_quat_rots( quats[ ctr_a ], quats[ ctr_b ] );
//    		const angle quat_angle_b_to_a = angle_between_quat_rots( quats[ ctr_b ], quats[ ctr_a ] );
//			BOOST_CHECK_CLOSE( angle_in_degrees( rotn_angle_a_to_b ), angle_in_degrees( quat_angle_a_to_b ), LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
//			BOOST_CHECK_CLOSE( angle_in_degrees( rotn_angle_b_to_a ), angle_in_degrees( quat_angle_b_to_a ), LOOSER_ACCURACY_PERCENTAGE_TMPL<quat_rot_type>() );
//    	}
//	}
//}

BOOST_AUTO_TEST_SUITE_END()
