/// \file
/// \brief The pdbs_acquirer test suite

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

#include "acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp"
#include "chopping/domain/domain.hpp"
#include "file/name_set/name_set_list.hpp"
#include "file/strucs_context.hpp"
#include "test/global_test_constants.hpp"

namespace cath{ namespace test { } }

using namespace cath;
using namespace cath::chop;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::test;

using std::istringstream;

namespace cath {
	namespace test {

		/// \brief The pdbs_acquirer_test_suite_fixture to assist in testing pdbs_acquirer
		struct pdbs_acquirer_test_suite_fixture {
		protected:
			~pdbs_acquirer_test_suite_fixture() noexcept = default;

			istringstream input_ss;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(pdbs_acquirer_test_suite, pdbs_acquirer_test_suite_fixture)

BOOST_AUTO_TEST_CASE(regions_names_make_it_into_strucs_context) {
	constexpr bool remove_partial_residues = false;

	const name_set_list got_name_sets = get_strucs_context(
		file_list_pdbs_acquirer{ {
			global_test_constants::TEST_SOURCE_DATA_DIR() / "1c0pA01",
			global_test_constants::TEST_SOURCE_DATA_DIR() / "1hdoA00",
		} },
		input_ss,
		remove_partial_residues,
		str_vec{},
		domain_vec{ {
			domain{ {}, "domain_1c0pA01" },
			domain{ {}, "domain_1hdoA00" },
		} }
	).get_name_sets();

	const str_vec expected = { "domain_1c0pA01", "domain_1hdoA00" };

	BOOST_TEST( get_domain_or_specified_or_from_acq_names( got_name_sets ) == expected );
}

BOOST_AUTO_TEST_SUITE_END()
