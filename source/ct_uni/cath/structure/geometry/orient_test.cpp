/// \file
/// \brief The orient test suite

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
/// MERCHANTABILITY or ORIENTNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <boost/test/unit_test.hpp>

#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/orient.hpp"
#include "cath/structure/geometry/rotation.hpp"
#include "cath/structure/structure_type_aliases.hpp"

using namespace ::cath::geom;
using namespace ::cath::geom::detail;

BOOST_AUTO_TEST_SUITE(orient_test_suite)

BOOST_AUTO_TEST_SUITE(x_and_y_of_later_weighted_cog_works)

BOOST_AUTO_TEST_CASE(x_and_y_of_later_weighted_cog_works_if_all_same) {
	constexpr coord the_coord{ 2.0, 4.0, 6.0 };
	const auto got = x_and_y_of_later_weighted_cog( coord_list{ coord_vec{ {
		the_coord,
		the_coord,
		the_coord
	} } } );
	BOOST_TEST( got.first  == the_coord.get_x() );
	BOOST_TEST( got.second == the_coord.get_y() );
}

BOOST_AUTO_TEST_CASE(x_and_y_of_later_weighted_cog_works_if_some_different) {
	constexpr coord the_coord{ 2.0, 4.0, 6.0 };
	const auto got = x_and_y_of_later_weighted_cog( coord_list{ coord_vec{ {
		the_coord,
		the_coord,
		ORIGIN_COORD
	} } } );
	BOOST_TEST( got.first  == the_coord.get_x() / 2.0 );
	BOOST_TEST( got.second == the_coord.get_y() / 2.0 );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(neg_neg) {
	const auto x = get_orienting_transformation( coord_list{ coord_vec{ {
		coord{  -32.2901,  -18.9915,   3.3418 },
		coord{  -26.2388,  -18.5794,  -2.9048 },
		coord{   33.2988,   15.2591,  -7.5908 },
		coord{    7.2787,    3.1599,  -3.6327 },
		coord{  -44.2691,  -26.1551,   3.3464 },
		coord{  -14.8504,   -7.9242,   1.5012 },
		coord{   -7.3366,    0.0089,   6.0726 },
		coord{   71.8374,   37.5763, -10.1231 },
		coord{  -31.7500,  -17.4188,   4.4648 },
		coord{   23.3887,   11.8416,  -3.2778 },
	} } } );
	// Initial orientation's x_and_y_of_later_weighted_cog : -6.83283 -0.933984
	// Flipped orientation's x_and_y_of_later_weighted_cog :  6.83283  0.933984
	BOOST_TEST( x.first  == coord( 2.09314, 2.12232, 0.88024 ) );
	BOOST_TEST( x.second == rotation(
		 0.87361,   0.474223, -0.10917,
		-0.170357,  0.50818,   0.844234,
		 0.455834, -0.718934,  0.524738
	) );
}

BOOST_AUTO_TEST_CASE(neg_pos) {
	const auto x = get_orienting_transformation( coord_list{ coord_vec{ {
		coord{   33.1842,   19.8389,  -0.6638 },
		coord{    1.0632,    0.3512,  -0.2627 },
		coord{  -46.5568,  -24.9492,   6.3272 },
		coord{  -81.6825,  -40.9361,  15.0224 },
		coord{   46.1277,   28.1789,   0.3956 },
		coord{  -72.6674,  -37.9115,  13.4391 },
		coord{   53.4100,   31.6977,  -4.1862 },
		coord{  -49.9652,  -27.4963,   6.9813 },
		coord{   31.5065,   18.2541,  -1.1976 },
		coord{   29.8915,   13.7668,  -6.5958 },
	} } } );
	// Initial orientation's x_and_y_of_later_weighted_cog : -4.60324  0.492968
	// Flipped orientation's x_and_y_of_later_weighted_cog :  4.60324  0.492968
	BOOST_TEST( x.first  == coord( 5.5688800000000036050096242, 1.9205500000000004234834705, -2.9259499999999998287592007 ) );
	BOOST_TEST( x.second == rotation(
		 0.8733088305730067890664259,  0.4739233843227982911905372, -0.1128198219961049353354809,
		 0.1474278279587467321842809, -0.4778232235515197379172037, -0.8659965372784123038840676,
		-0.4643239407854969913458376,  0.7396496419433054025915908, -0.4871567357499705930301559
	) );
}

BOOST_AUTO_TEST_CASE(pos_neg) {
	const auto x = get_orienting_transformation( coord_list{ coord_vec{ {
		coord{    6.2975,    4.6808,   1.0715 },
		coord{  -55.8766,  -31.7599,   4.2710 },
		coord{    6.0873,    3.5345,  -0.1998 },
		coord{   11.5486,    9.1254,   2.3147 },
		coord{  -58.1234,  -36.1967,   0.5316 },
		coord{   78.5637,   43.6476, -10.3017 },
		coord{   68.5768,   36.9250, -11.0521 },
		coord{  -15.7869,   -4.5773,   7.9909 },
		coord{  -50.3832,  -26.5726,   9.0812 },
		coord{  -26.1118,  -16.6516,  -0.3351 },
	} } } );
	// Initial orientation's x_and_y_of_later_weighted_cog :  0.58850 -0.030858
	// Flipped orientation's x_and_y_of_later_weighted_cog :  0.58850  0.030858
	BOOST_TEST( x.first  == coord( 3.5208, 1.78448, -0.33722 ) );
	BOOST_TEST( x.second == rotation(
		-0.867108, -0.488554,  0.0971536,
		-0.149826,  0.441811,  0.884508,
		-0.475054,  0.752408, -0.456296
	) );
}

BOOST_AUTO_TEST_CASE(pos_pos) {
	const auto x = get_orienting_transformation( coord_list{ coord_vec{ {
		coord{   32.7212,   18.2694,  -3.8942 },
		coord{    8.8652,    7.0475,   2.5630 },
		coord{   68.2275,   36.7794,  -8.8260 },
		coord{  -59.5060,  -30.6612,  11.7655 },
		coord{   19.2661,   13.6014,   1.9597 },
		coord{  -19.5175,   -8.3081,   7.5434 },
		coord{  -87.6071,  -45.6423,  13.5100 },
		coord{  -57.8828,  -31.4848,   6.1705 },
		coord{   77.3277,   40.2418, -12.5506 },
		coord{   58.2665,   32.4735,  -7.2997 },
	} } } );
	// Initial orientation's x_and_y_of_later_weighted_cog :  0.49738  0.563029
	// Flipped orientation's x_and_y_of_later_weighted_cog :  0.49738  0.563029
	BOOST_TEST( x.first  == coord( -4.01608,  -3.23166,  -1.09416 ) );
	BOOST_TEST( x.second == rotation(
		-0.8756499359656251257888471, -0.4652124569033634360337714,  0.1296709666245027814390767,
		 0.0869147098253301580994545, -0.4159177532145612032898896, -0.9052393361851495123815425,
		 0.4750610727678866718193262, -0.7814024523271523303691311,  0.4046321596681818899554628
	) );
}

BOOST_AUTO_TEST_SUITE_END()
