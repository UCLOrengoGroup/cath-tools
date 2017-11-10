/// \file
/// \brief The rbf_model test suite

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

#include <boost/filesystem.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/algorithm/contains.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/file/open_fstream.hpp"
#include "file/prc_scores_file/prc_scores_entry.hpp"
#include "file/prc_scores_file/prc_scores_file.hpp"
#include "file/ssap_scores_file/ssap_scores_entry.hpp"
#include "file/ssap_scores_file/ssap_scores_file.hpp"
#include "score/score_classification/rbf_model.hpp"
#include "test/global_test_constants.hpp"

#include <fstream>
#include <utility>

using namespace cath::common;
using namespace cath::file;
using namespace cath::score;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The rbf_model_test_suite_fixture to assist in testing rbf_model
		struct rbf_model_test_suite_fixture : protected global_test_constants {
		protected:
			~rbf_model_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(rbf_model_test_suite, cath::test::rbf_model_test_suite_fixture)

BOOST_AUTO_TEST_CASE(parsed_model_calculates_correct_score) {
	const auto the_model = parse_rbf_model( TEST_SVM_DIR() / "cath_svm.rbf_gamma_1_c_5.model" );

	const auto the_prc_scores_entry  = prc_scores_entry_from_line ( "1by5A02     2       554     554     1       1xkhA02 7       543     543      223.8   215.5  3.1e-185" );
	const auto the_ssap_scores_entry = ssap_scores_entry_from_line( "1by5A02 1xkhA02 554 535 79.16 511 92 17 3.11" );
	BOOST_CHECK_CLOSE( get_score( the_model, the_prc_scores_entry, the_ssap_scores_entry ), 8.3668143327915736, ACCURACY_PERCENTAGE() );
}

BOOST_AUTO_TEST_SUITE_END()
