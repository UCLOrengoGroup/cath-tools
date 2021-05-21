/// \file
/// \brief The vcie_match_criteria test suite

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

#include "vcie_match_criteria.hpp"

#include "cath/common/size_t_literal.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath::common;
using namespace ::cath::index;
using namespace ::cath::index::detail;
using namespace ::cath::geom;

namespace cath {
	namespace test {

		/// \brief The vcie_match_criteria_test_suite_fixture to assist in testing vcie_match_criteria
		struct vcie_match_criteria_test_suite_fixture {
		protected:
			~vcie_match_criteria_test_suite_fixture() noexcept = default;

		    const angle_type the_zero_angle = ZERO_ANGLE<angle_base_type>;
		};

	}  // namespace test
}  // namespace cath

/// \brief Unit test vcie_match_criteria
BOOST_FIXTURE_TEST_SUITE(vcie_match_criteria_test_suite, cath::test::vcie_match_criteria_test_suite_fixture)

/// \brief Check that the default vcie_match_criteria built by make_default_vcie_match_criteria() returns the expected values from the getters
BOOST_AUTO_TEST_CASE(make_default_works) {
	const auto default_vcie_match_criteria = make_default_vcie_match_criteria();
	BOOST_CHECK_EQUAL( default_vcie_match_criteria.get_require_matching_directions(),    true                                             );
	BOOST_CHECK_EQUAL( default_vcie_match_criteria.get_minimum_index_distance(),         11_z                                             );
	BOOST_CHECK_EQUAL( default_vcie_match_criteria.get_maximum_squared_distance(),       40.0                                             );
	BOOST_CHECK_EQUAL( default_vcie_match_criteria.get_maximum_frame_angle_difference(), make_angle_from_degrees<angle_base_type>( 22.5 ) );
	BOOST_CHECK_EQUAL( default_vcie_match_criteria.get_maximum_phi_angle_difference(),   make_angle_from_degrees<angle_base_type>( 67.5 ) );
	BOOST_CHECK_EQUAL( default_vcie_match_criteria.get_maximum_psi_angle_difference(),   make_angle_from_degrees<angle_base_type>( 67.5 ) );
}

/// \brief Check that get_standard_vcie_match_criterias() successfully returns 141 vcie_match_criteria objects
BOOST_AUTO_TEST_CASE(get_standard_works) {
	const vcie_match_criteria_vec standard_criterias = get_standard_vcie_match_criterias();
	BOOST_CHECK_EQUAL( standard_criterias.size(), 141_z );
}

/// \brief Check that the vcie_match_criteria function operator for a single view_cache_index_entry works as expected
BOOST_AUTO_TEST_CASE(singe_view_cache_index_entry_check_works) {
	const auto default_vcie_match_criteria = make_default_vcie_match_criteria();
	BOOST_CHECK_EQUAL( default_vcie_match_criteria( view_cache_index_entry( 10, 20, view_type( ORIGIN_COORD ), IDENTITY_ROTATION, the_zero_angle, the_zero_angle, the_zero_angle, the_zero_angle ) ), false );
	BOOST_CHECK_EQUAL( default_vcie_match_criteria( view_cache_index_entry( 10, 21, view_type( ORIGIN_COORD ), IDENTITY_ROTATION, the_zero_angle, the_zero_angle, the_zero_angle, the_zero_angle ) ), true  );
}

/// \todo Test the vcie_match_criteria function operator that tests pairs of view_cache_index_entries

BOOST_AUTO_TEST_SUITE_END()
