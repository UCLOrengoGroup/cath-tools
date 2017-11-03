/// \file
/// \brief The proximity_calculator test suite

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

#include <boost/test/auto_unit_test.hpp>

#include "chopping/region/region.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/proximity_calculator.hpp"
#include "structure/geometry/coord.hpp"

using namespace cath::file;
using namespace cath::geom;

using std::string;
using std::istringstream;

BOOST_AUTO_TEST_SUITE(proximity_calculator_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	// Example is from chain B of 4c9a
	const string input_string = R"(ATOM   1861  N   THR B 138     -33.417  42.721 103.639  1.00142.96           N  
ATOM   1862  CA  THR B 138     -33.726  41.382 103.094  1.00124.05           C  
ATOM   1863  C   THR B 138     -34.247  41.593 101.682  1.00122.85           C  
ATOM   1864  O   THR B 138     -34.642  42.709 101.356  1.00126.34           O  
ATOM   1865  CB  THR B 138     -34.816  40.679 103.927  1.00108.88           C  
ATOM   1866  N   MSE B 139     -34.243  40.568 100.830  1.00114.60           N  
ATOM   1867  CA  MSE B 139     -34.924  40.699  99.523  1.00114.24           C  
ATOM   1868  C   MSE B 139     -36.170  39.811  99.512  1.00107.49           C  
ATOM   1869  O   MSE B 139     -36.514  39.225  98.475  1.00116.82           O  
ATOM   1870  CB  MSE B 139     -33.991  40.377  98.339  1.00105.83           C  
ATOM   1871  N   VAL B 140     -36.867  39.765 100.654  1.00 94.64           N  
ATOM   1872  CA  VAL B 140     -37.904  38.745 100.917  1.00101.80           C  
ATOM   1873  C   VAL B 140     -39.124  39.279 101.672  1.00100.30           C  
ATOM   1874  O   VAL B 140     -39.041  40.306 102.330  1.00 97.38           O  
ATOM   1875  CB  VAL B 140     -37.331  37.579 101.757  1.00 93.10           C)";
	istringstream input_ss{ input_string };
	const pdb the_pdb = read_pdb_file( input_ss );
	BOOST_REQUIRE_EQUAL( the_pdb.get_num_residues(), 3 );

	const proximity_calculator the_prox{ the_pdb, {} };

	BOOST_CHECK( ! the_prox.is_within_distance( coord{ -33.726, 41.382, 104.094 }, 0.9 ) );
	BOOST_CHECK(   the_prox.is_within_distance( coord{ -33.726, 41.382, 104.094 }, 1.1 ) );
	BOOST_CHECK( ! the_prox.is_within_distance( coord{ -33.991, 40.377,  97.339 }, 0.9 ) );
	BOOST_CHECK(   the_prox.is_within_distance( coord{ -33.991, 40.377,  97.339 }, 1.1 ) );
}
BOOST_AUTO_TEST_SUITE_END()
