/// \file
/// \brief The superpose_fit test suite

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

namespace cath { namespace test { } }

#include "structure/bioplib_facade/bioplib_interface.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/geometry/rotation.hpp"
#include "superposition/superpose_fit.hpp"

using namespace cath;
using namespace cath::geom;
using namespace cath::sup;
using namespace cath::test;

namespace cath {
	namespace test {

		/// \brief The superpose_fit_test_suite_fixture to assist in testing superpose_fit
		struct superpose_fit_test_suite_fixture {
		protected:
			~superpose_fit_test_suite_fixture() noexcept = default;

		public:

		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(superpose_fit_test_suite, superpose_fit_test_suite_fixture)

BOOST_AUTO_TEST_CASE(issue_34_test) {
	// const coord_list a{ coord_vec{
	// 	coord{  -5.862210,  -3.998580, -11.220700    },
	// 	coord{  -6.901210,  -0.330583, -10.636700    },
	// 	coord{  -8.304210,  -1.434580,  -7.228670    },
	// 	coord{  -5.069210,  -3.148580,  -5.878670    },
	// 	coord{  -3.018210,  -0.227583,  -7.333670    },
	// 	coord{  -5.274210,   2.282420,  -5.409670    },
	// 	coord{  -4.637210,   0.101417,  -2.319670    },
	// 	coord{  -0.778208,   0.289417,  -2.449670    },
	// 	coord{  -0.849208,   3.974420,  -3.662670    },
	// 	coord{  -3.040210,   4.998420,  -0.662667    },
	// 	coord{  -0.589208,   3.006420,   1.567330    },
	// 	coord{   2.669790,   4.732420,   0.311333    },
	// 	coord{   0.963792,   8.236420,   0.423333    },
	// 	coord{  -0.196208,   7.353420,   3.913330    },
	// 	coord{   5.319790,   3.521420,   9.306330    },
	// 	coord{   7.179790,   0.206417,   9.983330    },
	// 	coord{   3.988790,  -1.972580,   9.930330    },
	// 	coord{   2.740790,  -0.206583,   6.797330    },
	// 	coord{   6.203790,  -1.177580,   5.304330    },
	// 	coord{   5.273790,  -4.796580,   6.322330    },
	// 	coord{   1.675790,  -4.548580,   4.877330    },
	// 	coord{   3.367790,  -3.157580,   1.752330    },
	// 	coord{   1.059790,  -6.866580,  -0.658667    },
	// 	coord{   4.075790,  -6.836580,  -3.027670    },
	// } };
	// const coord_list b{ coord_vec{
	// 	coord{   4.089940,  -6.717720,   1.227690    },
	// 	coord{   4.298780,  -5.580620,  -2.388410    },
	// 	coord{   1.480720,  -3.052560,  -1.916530    },
	// 	coord{   3.226990,  -1.661550,   1.177890    },
	// 	coord{   6.496840,  -1.252130,  -0.705125    },
	// 	coord{   4.744580,   0.561099,  -3.545520    },
	// 	coord{   2.831770,   2.687990,  -1.053820    },
	// 	coord{   6.053870,   3.870940,   0.543387    },
	// 	coord{   7.459100,   4.774870,  -2.873160    },
	// 	coord{   4.290960,   6.810250,  -3.353990    },
	// 	coord{   4.693920,   8.544250,   0.00280986  },
	// 	coord{   8.308070,   9.480560,  -0.766366    },
	// 	coord{   7.136960,  10.941100,  -4.067790    },
	// 	coord{   4.512280,  13.089600,  -2.382770    },
	// 	coord{  -4.282620,   2.562800,   0.210617    },
	// 	coord{  -7.637030,   1.398250,  -1.120400    },
	// 	coord{  -6.113330,  -1.881040,  -2.325910    },
	// 	coord{  -4.428670,  -2.337410,   1.054550    },
	// 	coord{  -7.774610,  -1.970400,   2.838860    },
	// 	coord{  -9.399960,  -4.480000,   0.516275    },
	// 	coord{  -6.416500,  -6.806050,   1.001890    },
	// 	coord{  -6.770300,  -6.451230,   4.775900    },
	// 	coord{  -7.174810, -11.553500,   5.146210    },
	// 	coord{  -9.626950, -10.977500,   8.003710    },
	// } };
	// const auto bioplib_rot = bioplib_fit( a, b );

	// const auto svd_rot     = superpose_fit( a, b );

	// std::cerr << "bioplib_rot : " << bioplib_rot << "\n";
	// // ignore_unused( svd_rot );
	// std::cerr << "svd_rot     : " << svd_rot     << "\n";

	// std::cerr << "a                         : " << a                         << "\n\n";
	// std::cerr << "rotate_copy( svd_rot, b ) : " << rotate_copy( svd_rot, b ) << "\n\n";
	// std::cerr << "RMSD ( bioplib_rot, b, a ) : " << calc_rmsd( rotate_copy( bioplib_rot, b ), a ) << "\n\n";
	// std::cerr << "RMSD ( svd_rot,     b, a ) : " << calc_rmsd( rotate_copy( svd_rot,     b ), a ) << "\n\n";

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()

