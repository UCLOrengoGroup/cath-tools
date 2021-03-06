/// \file
/// \brief The classn_rate_stat test suite

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

#include "cath/common/algorithm/constexpr_for_n.hpp"
#include "cath/score/true_pos_false_neg/classn_rate_stat.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::cath::score::detail;
using namespace ::std;

namespace {

	/// \brief TODOCUMENT
	template <size_t I>
	class std_classn_rate_stat_tester final {
	public:
		/// \brief TODOCUMENT
		void operator()() {
			constexpr std_classn_rate_stat to_test = get<0>( get<I>( properties_of_classn_rate_stat::numerator_and_denominator_of_stat ) );
			BOOST_CHECK_NO_THROW_DIAG( classn_rate_stat<to_test>() );
			classn_rate_stat<to_test> stat;
			BOOST_CHECK_NO_THROW_DIAG( [[maybe_unused]] auto &&x = stat.get_name() );
			BOOST_CHECK_GT( stat.get_name().length(), 0 );
		}
	};

} // namespace

/// \brief TODOCUMENT
BOOST_AUTO_TEST_SUITE(classn_rate_stat_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ctor_does_not_throw_and_has_non_empty_name) {
	constexpr_for_n<std_classn_rate_stat_tester, properties_of_classn_rate_stat::num_std_classn_rate_stats>();
}

BOOST_AUTO_TEST_SUITE_END()
