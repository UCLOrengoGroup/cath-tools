/// \file
/// \brief The restrict_to_single_linkage_extension test suite

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

#include "common/algorithm/sort_copy.hpp"
#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/restrict_to_single_linkage_extension.hpp"

using namespace cath::common;
using namespace cath::geom;

using std::tie;

BOOST_AUTO_TEST_SUITE(restrict_to_single_linkage_extension_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	const auto coord_less_fn = [] (const coord &lhs, const coord &rhs) {
		return (
			tie( lhs.get_x(), lhs.get_y(), lhs.get_z() )
			<
			tie( rhs.get_x(), rhs.get_y(), rhs.get_z() )
		);
	};

	const auto got = sort_copy<coord_vec>(
		restrict_to_single_linkage_extension_copy(
			coord_vec{
				coord{ 58.308,  3.632,  -9.326 },
				coord{ 58.697,  5.944,  -9.915 },
				coord{ 47.776, -0.017, -14.037 },
				coord{ 58.266,  6.731, -10.952 },
				coord{ 26.642,  7.586,   2.822 },
				coord{ 60.053,  4.597,   0.550 },
				coord{ 56.863,  2.610, -10.925 },
				coord{ 77.850, 10.964, -22.808 },
				coord{ 54.392, -3.073,  -5.950 },
				coord{ 56.649,  3.707, -11.743 },
				coord{ 65.828, 22.140,  -7.057 },
				coord{ 64.871,  9.909,   1.336 },
				coord{ 47.988,  6.662,   5.622 },
				coord{ 55.888,  3.583, -12.717 },
				coord{ 58.133,  4.709, -10.123 },
				coord{ 57.382,  4.823, -11.271 },
				coord{ 57.495,  6.104, -11.787 },
				coord{ 68.876, 26.299, -17.868 },
				coord{ 57.638,  2.594,  -9.790 },
				coord{ 88.223, 19.643, -13.873 },
			},
			1,
			4.2
		),
		coord_less_fn
	);

	const coord_vec expected = {
		coord{ 55.888, 3.583, -12.717 },
		coord{ 56.649, 3.707, -11.743 },
		coord{ 56.863, 2.610, -10.925 },
		coord{ 57.382, 4.823, -11.271 },
		coord{ 57.495, 6.104, -11.787 },
		coord{ 57.638, 2.594,  -9.790 },
		coord{ 58.133, 4.709, -10.123 },
		coord{ 58.266, 6.731, -10.952 },
		coord{ 58.308, 3.632,  -9.326 },
		coord{ 58.697, 5.944,  -9.915 },
	};

	BOOST_CHECK_EQUAL_RANGES( got, expected );
}

BOOST_AUTO_TEST_SUITE_END()

