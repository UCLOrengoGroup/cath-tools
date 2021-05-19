/// \file
/// \brief The display_colour_spec test suite

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

#include "cath/common/size_t_literal.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"
#include "cath/test/global_test_constants.hpp"
#include "display_colour_spec.hpp"

using namespace ::cath;
using namespace ::cath::common;

namespace cath {
	namespace test {

		/// \brief The display_colour_spec_test_suite_fixture to assist in testing display_colour_spec
		struct display_colour_spec_test_suite_fixture : protected global_test_constants {
		protected:
			~display_colour_spec_test_suite_fixture() noexcept = default;

		  public:
			[[nodiscard]] display_colour_spec get_spec() const;

			const display_colour_spec the_spec{ get_spec() };
			const display_colour_vec  pdb_colours{ { BLUE, RED } };
			const display_colour_vec  residue_colours{ { YELLOW } };
			const display_colour_vec all_colours{ { BLUE, RED, YELLOW } };
			const size_vec           red_pdbs{ 1 };
			const size_vec           blue_pdbs{ 3 };
			const size_vec           pdb_two_residues{ { 4, 6 } };
		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
display_colour_spec cath::test::display_colour_spec_test_suite_fixture::get_spec() const {
	display_colour_spec temp_spec;
	temp_spec.colour_pdb        ( 1,    RED    );
	temp_spec.colour_pdb        ( 3,    BLUE   );
	temp_spec.colour_pdb_residue( 2, 4, YELLOW );
	temp_spec.colour_pdb_residue( 2, 6, YELLOW );
	return temp_spec;
}


/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(display_colour_spec_test_suite, cath::test::display_colour_spec_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(pdb_colours_works) {
	const display_colour_vec got_pdb_colours = get_pdb_colours( the_spec );
	BOOST_CHECK_EQUAL_RANGES( pdb_colours, got_pdb_colours );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(residue_colours_works) {
	const display_colour_vec got_residue_colours = get_residue_colours( the_spec );
	BOOST_CHECK_EQUAL_RANGES( residue_colours, got_residue_colours );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(all_colours_works) {
	const display_colour_vec got_all_colours = get_all_colours( the_spec );
	BOOST_CHECK_EQUAL_RANGES( all_colours, got_all_colours );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(pdbs_of_colour_works) {
	const size_vec got_red_pdbs = get_pdbs_of_colour( the_spec, RED );
	BOOST_CHECK_EQUAL_RANGES( red_pdbs, got_red_pdbs );
	const size_vec got_blue_pdbs = get_pdbs_of_colour( the_spec, BLUE );
	BOOST_CHECK_EQUAL_RANGES( blue_pdbs, got_blue_pdbs );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(residues_of_colour_works) {
	const size_size_vec_map got_residues = get_residues_of_colour( the_spec, YELLOW );
	BOOST_REQUIRE_EQUAL( 1_z, got_residues.size()         );
	BOOST_CHECK_EQUAL  ( 2_z, cbegin( got_residues )->first  );
	const size_vec &got_pdb_two_residues = cbegin( got_residues )->second;
	BOOST_CHECK_EQUAL_RANGES( pdb_two_residues, got_pdb_two_residues );
}

BOOST_AUTO_TEST_SUITE_END()
