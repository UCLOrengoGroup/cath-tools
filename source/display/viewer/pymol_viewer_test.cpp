/// \file
/// \brief The pymol_viewer_test test suite

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

#include "display/viewer/pymol_viewer.hpp"
#include "structure/residue_name.hpp"

using namespace cath;

BOOST_AUTO_TEST_SUITE(pymol_viewer_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL(
		pymol_viewer{}.get_colour_pdb_residues_str( "red", "1cukA", { residue_name{ 1 }, residue_name{ 2 } } ),
		"colour red, /1cukA///1+2/\n"
	);
}
BOOST_AUTO_TEST_SUITE_END()

