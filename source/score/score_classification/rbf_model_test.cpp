/// \file
/// \brief The rbf_model test suite

#include <boost/filesystem.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/irange.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/algorithm/contains.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/file/open_fstream.hpp"
#include "common/size_t_literal.hpp"
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

using boost::filesystem::directory_iterator;
using boost::filesystem::path;
using boost::irange;
using boost::range::sort;

namespace cath {
	namespace test {

		/// \brief The rbf_model_test_suite_fixture to assist in testing rbf_model
		struct rbf_model_test_suite_fixture : protected global_test_constants {
		protected:
			~rbf_model_test_suite_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(rbf_model_test_suite, cath::test::rbf_model_test_suite_fixture)

BOOST_AUTO_TEST_CASE(parsed_model_calculates_correct_score) {
	const auto the_model = parse_rbf_model( TEST_SVM_DIR() / "cath_svm.rbf_gamma_1_c_5.model" );

	const auto the_prc_scores_entry  = prc_scores_entry_from_line ( "1by5A02     2       554     554     1       1xkhA02 7       543     543      223.8   215.5  3.1e-185" );
	const auto the_ssap_scores_entry = ssap_scores_entry_from_line( "1by5A02 1xkhA02 554 535 79.16 511 92 17 3.11" );
	BOOST_CHECK_CLOSE( get_score( the_model, the_prc_scores_entry, the_ssap_scores_entry ), 8.3668143327915736, ACCURACY_PERCENTAGE() );
}

BOOST_AUTO_TEST_SUITE_END()
