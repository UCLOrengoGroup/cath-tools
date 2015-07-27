/// \file
/// \brief The quad_criteria test suite

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

#include "quad_criteria.h"

#include "common/size_t_literal.h"
#include "scan/detail/quad_criteria_are_met_by.h"
#include "scan/detail/res_pair/single_struc_res_pair.h"
#include "structure/geometry/angle.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

using namespace cath::common;
using namespace cath::scan;
using namespace cath::scan::detail;
using namespace cath::geom;

namespace cath {
	namespace test {

		/// \brief The quad_criteria_test_suite_fixture to assist in testing quad_criteria
		struct quad_criteria_test_suite_fixture {
		protected:
			~quad_criteria_test_suite_fixture() noexcept = default;

		    const angle_type the_zero_angle = zero_angle<angle_base_type>();
		};

	}
}

/// \brief Unit test quad_criteria
BOOST_FIXTURE_TEST_SUITE(quad_criteria_test_suite, cath::test::quad_criteria_test_suite_fixture)

/// \brief Check that the default quad_criteria built by make_default_quad_criteria() returns the expected values from the getters
BOOST_AUTO_TEST_CASE(make_default_works) {
	const auto default_quad_criteria = make_default_quad_criteria();
	BOOST_CHECK_EQUAL( default_quad_criteria.get_index_direction_criterion(),      res_pair_index_dirn_criterion::MUST_MATCH        );
	BOOST_CHECK_EQUAL( default_quad_criteria.get_minimum_index_distance(),         11_z                                             );
	BOOST_CHECK_EQUAL( default_quad_criteria.get_maximum_squared_distance(),       40.0                                             );
	BOOST_CHECK_EQUAL( default_quad_criteria.get_maximum_frame_angle_difference(), make_angle_from_degrees<angle_base_type>( 22.5 ) );
	BOOST_CHECK_EQUAL( default_quad_criteria.get_maximum_phi_angle_difference(),   make_angle_from_degrees<angle_base_type>( 67.5 ) );
	BOOST_CHECK_EQUAL( default_quad_criteria.get_maximum_psi_angle_difference(),   make_angle_from_degrees<angle_base_type>( 67.5 ) );
}

/// \brief Check that get_standard_quad_criterias() successfully returns 141 quad_criteria objects
BOOST_AUTO_TEST_CASE(get_standard_works) {
	const quad_criteria_vec standard_criterias = get_standard_quad_criterias();
	BOOST_CHECK_EQUAL( standard_criterias.size(), 141_z );
}

/// \brief Check that the quad_criteria function operator for a single res_pair works as expected
BOOST_AUTO_TEST_CASE(singe_res_pair_check_works) {
	const auto default_quad_criteria = make_default_quad_criteria();
	BOOST_CHECK_EQUAL( are_not_violated_by( default_quad_criteria, single_struc_res_pair( coord::ORIGIN_COORD, make_identity_quat_rot<frame_quat_rot_type>(), the_zero_angle, the_zero_angle, the_zero_angle, the_zero_angle, 10, 20 ) ), false );
	BOOST_CHECK_EQUAL( are_not_violated_by( default_quad_criteria, single_struc_res_pair( coord::ORIGIN_COORD, make_identity_quat_rot<frame_quat_rot_type>(), the_zero_angle, the_zero_angle, the_zero_angle, the_zero_angle, 10, 21 ) ), true  );
}

/// \todo Test the quad_criteria function operator that tests pairs of res_pairs

BOOST_AUTO_TEST_SUITE_END()
