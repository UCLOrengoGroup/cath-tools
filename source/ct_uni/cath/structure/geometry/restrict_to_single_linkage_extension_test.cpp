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

#include <boost/range/adaptor/transformed.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/algorithm/sort_uniq_build.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/restrict_to_single_linkage_extension.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"

namespace cath { namespace test { } }

using namespace cath::common;
using namespace cath::geom;
using namespace cath::test;

using boost::adaptors::transformed;
using std::make_pair;
using std::tie;

namespace cath {
	namespace test {

		/// \brief The restrict_to_single_linkage_extension_test_suite_fixture to assist in testing restrict_to_single_linkage_extension
		struct restrict_to_single_linkage_extension_test_suite_fixture {
		protected:
			~restrict_to_single_linkage_extension_test_suite_fixture() noexcept = default;

			/// \brief Check that the members of the specified restricted coord_coord_linkage_pair_vec match the
			///        the expected list (though not necessarily in the same order; they both get sorted before comparison)
			void compare_restricted(const coord_coord_linkage_pair_vec &prm_got_unsorted,     ///< The got restricted coord_coord_linkage_pair_vec (doesn't need to be pre-sorted)
			                        coord_vec                           prm_expected_unsorted ///< The expected coord_vec (doesn't need to be pre-sorted)
			                        ) {
				// A function object to sort coordinates to allow for easy equality testing of coord_vecs
				const auto coord_less_fn = [] (const coord &lhs, const coord &rhs) {
					return (
						tie( lhs.get_x(), lhs.get_y(), lhs.get_z() )
						<
						tie( rhs.get_x(), rhs.get_y(), rhs.get_z() )
					);
				};

				// Get a sorted version of the two lists of coords
				const auto got_sorted = sort_build<coord_vec>(
					prm_got_unsorted
						| transformed( [] (const coord_coord_linkage_pair &x) { return x.first; } ),
					coord_less_fn
				);
				const auto expected_sorted = sort_copy(
					std::move( prm_expected_unsorted ),
					coord_less_fn
				);

				// Check got matches expected
				BOOST_CHECK_EQUAL_RANGES( got_sorted, expected_sorted );
			}

		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(restrict_to_single_linkage_extension_test_suite, restrict_to_single_linkage_extension_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	compare_restricted(
		restrict_to_single_linkage_extension_copy(
			coord_coord_linkage_pair_vec{ {
				make_pair( coord{ 58.308,  3.632,  -9.326 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 58.697,  5.944,  -9.915 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 47.776, -0.017, -14.037 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 58.266,  6.731, -10.952 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 26.642,  7.586,   2.822 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 60.053,  4.597,   0.550 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 56.863,  2.610, -10.925 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 77.850, 10.964, -22.808 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 54.392, -3.073,  -5.950 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 56.649,  3.707, -11.743 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 65.828, 22.140,  -7.057 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 64.871,  9.909,   1.336 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 47.988,  6.662,   5.622 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 55.888,  3.583, -12.717 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 58.133,  4.709, -10.123 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 57.382,  4.823, -11.271 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 57.495,  6.104, -11.787 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 68.876, 26.299, -17.868 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 57.638,  2.594,  -9.790 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 88.223, 19.643, -13.873 }, coord_linkage::ADD_AND_LINK ),
			} },
			1,
			4.2
		),
		coord_vec{
			coord{ 55.888,  3.583, -12.717 },
			coord{ 56.649,  3.707, -11.743 },
			coord{ 56.863,  2.610, -10.925 },
			coord{ 57.382,  4.823, -11.271 },
			coord{ 57.495,  6.104, -11.787 },
			coord{ 57.638,  2.594,  -9.790 },
			coord{ 58.133,  4.709, -10.123 },
			coord{ 58.266,  6.731, -10.952 },
			coord{ 58.308,  3.632,  -9.326 },
			coord{ 58.697,  5.944,  -9.915 },
		}
	);
}

BOOST_AUTO_TEST_CASE(add_only_does_not_daisy_chain) {
	compare_restricted(
		restrict_to_single_linkage_extension_copy(
			coord_coord_linkage_pair_vec{ {
				make_pair( coord{ 0.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 1.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 2.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 3.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 4.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 5.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 6.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 7.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 8.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
				make_pair( coord{ 9.0, 0.0, 0.0 }, coord_linkage::ADD_ONLY ),
			} },
			1,
			4.2
		),
		coord_vec{
			coord{ 0.0, 0.0, 0.0 },
			coord{ 1.0, 0.0, 0.0 },
			coord{ 2.0, 0.0, 0.0 },
			coord{ 3.0, 0.0, 0.0 },
			coord{ 4.0, 0.0, 0.0 },
		}
	);
}

BOOST_AUTO_TEST_CASE(add_and_link_does_daisy_chain) {
	compare_restricted(
		restrict_to_single_linkage_extension_copy(
			coord_coord_linkage_pair_vec{ {
				make_pair( coord{ 0.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 1.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 2.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 3.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 4.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 5.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 6.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 7.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 8.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
				make_pair( coord{ 9.0, 0.0, 0.0 }, coord_linkage::ADD_AND_LINK ),
			} },
			1,
			4.2
		),
		coord_vec{
			coord{ 0.0, 0.0, 0.0 },
			coord{ 1.0, 0.0, 0.0 },
			coord{ 2.0, 0.0, 0.0 },
			coord{ 3.0, 0.0, 0.0 },
			coord{ 4.0, 0.0, 0.0 },
			coord{ 5.0, 0.0, 0.0 },
			coord{ 6.0, 0.0, 0.0 },
			coord{ 7.0, 0.0, 0.0 },
			coord{ 8.0, 0.0, 0.0 },
			coord{ 9.0, 0.0, 0.0 },
		}
	);
}

BOOST_AUTO_TEST_SUITE_END()

