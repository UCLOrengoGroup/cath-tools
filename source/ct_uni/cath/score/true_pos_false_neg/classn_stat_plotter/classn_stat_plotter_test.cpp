/// \file
/// \brief The classn_stat_plotter test suite

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

#include "cath/common/file/temp_file.hpp"
#include "cath/score/score_classification/score_classn_value_list.hpp"
#include "cath/score/true_pos_false_neg/classn_rate_stat.hpp"
#include "cath/score/true_pos_false_neg/classn_stat_pair_series.hpp"
#include "cath/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter.hpp"
#include "cath/score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_spec.hpp"
#include "cath/score/true_pos_false_neg/named_true_false_pos_neg_list.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::score;

namespace {

	/// \brief The classn_stat_plotter_test_suite_fixture to assist in testing classn_stat_plotter
	struct classn_stat_plotter_test_suite_fixture: protected global_test_constants {
	protected:
		~classn_stat_plotter_test_suite_fixture() noexcept = default;

		/// \brief TODOCUMENT
		temp_file temp_base_file{ ".classn_stat_plotter_test_file.%%%%-%%%%-%%%%-%%%%" };

		/// \brief TODOCUMENT
		const score_classn_value_list the_score_classn_values{ make_score_classn_value_list(
			{
				score_classn_value{ 3.0,  true,  "a" },
				score_classn_value{ 3.1, false,  "b" },
				score_classn_value{ 3.2,  true,  "c" },
				score_classn_value{ 3.3,  true,  "d" },
				score_classn_value{ 3.4,  true, "e1" },
				score_classn_value{ 3.4, false, "e2" },
				score_classn_value{ 3.5,  true,  "f" },
				score_classn_value{ 3.6, false,  "g" },
				score_classn_value{ 3.7, false,  "h" },
				score_classn_value{ 3.8,  true,  "i" },
				score_classn_value{ 3.9, false,  "j" }
			},
			false,
			"Test Series"
		) };

	};

} // namespace

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(classn_stat_plotter_test_suite, classn_stat_plotter_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(check_area_under_roc_curve) {
	BOOST_CHECK_EQUAL(
		area_under_roc_curve( the_score_classn_values ),
		0.6833333333333333333333333333333333333333
	);
}

// Commenting-out because this pulls in a bit of gnuplot-iostream.h that needs boost::filesystem
//
// BOOST_AUTO_TEST_CASE(plot_does_not_throw) {
// 	BOOST_CHECK_NO_THROW_DIAG(
// 		plot_roc(
// 			classn_stat_plotter( true ),
// 			*temp_base_file.get_opt_filename(),
// 			make_named_true_false_pos_neg_list( the_score_classn_values ),
// 			make_standard_score_roc_plotter_spec( { } )
// 		)
// 	);
// }

BOOST_AUTO_TEST_SUITE_END()
