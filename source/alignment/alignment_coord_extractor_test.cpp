/// \file


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

#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "test/global_test_constants.hpp"
#include "structure/bioplib_facade/bioplib_pdb.hpp"
#include "structure/geometry/coord.hpp"

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The alignment_coord_extractor_test_suite_fixture to assist in testing alignment_coord_extractor
		struct alignment_coord_extractor_test_suite_fixture : protected global_test_constants {
		protected:
			~alignment_coord_extractor_test_suite_fixture() noexcept = default;

			const size_t    NUM_RESIDUES_TO_CHECK = { 5 };
			const coord_vec EXPECTED_FIRST_CA_COORDS = {
				{ 0.374, 147.376, 55.166 },
				{ 1.193, 145.597, 51.902 },
				{ 1.645, 141.900, 51.218 },
				{ 5.337, 141.415, 50.367 },
				{ 7.108, 138.545, 48.773 },
			};
		};

	}  // namespace test
}  // namespace cath

/// \brief Check that alignment_coord_extractor does as expected
///
/// At present, this doesn't really test alignment_coord_extractor but just checks that manually
/// manually grabbing CA coordinates in the same way as these two subroutines does as expected.
///
/// \todo Add a test for cath::protein (which will probably require tidying up code to load data
///       into a cath::protein object.
///
/// \todo Abstract out CA-grabbing functionality, perhaps into a decorator, so there is just
///       one alignment_coord_extractor subroutine that is given a grabber as an argument.
///
/// ...then...
///
/// \todo Move the existing test(s) here to test that new functionality
///
/// ...then...
///
/// \todo Write a proper test that checks At present, this doesn't really test alignment_coord_extractor but that the checks
///       manually

BOOST_FIXTURE_TEST_SUITE(alignment_coord_extractor_test_suite, cath::test::alignment_coord_extractor_test_suite_fixture)

/// \brief Check that grabbing the first few carbon alpha coordinates from a test file produces the expected results
BOOST_AUTO_TEST_CASE(bioplib_pdb_ca_coord_querying) {
	bioplib_pdb my_pdb;
	my_pdb.read_file(EXAMPLE_A_PDB_FILENAME().string());
	coord_vec first_ca_coords;
	first_ca_coords.reserve(EXPECTED_FIRST_CA_COORDS.size());
	for (size_t residue_ctr = 0; residue_ctr < EXPECTED_FIRST_CA_COORDS.size(); ++residue_ctr) {
		first_ca_coords.push_back(my_pdb.get_residue_ca_coord_of_index__backbone_unchecked(residue_ctr));
	}
	BOOST_CHECK_EQUAL_RANGES( EXPECTED_FIRST_CA_COORDS, first_ca_coords );
}

BOOST_AUTO_TEST_SUITE_END()
