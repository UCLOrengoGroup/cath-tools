/// \file
/// \brief The masked_bests_cache test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "resolve_hits/algo/discont_hits_index_by_start.hpp"
#include "resolve_hits/algo/masked_bests_cache.hpp"
#include "resolve_hits/algo/masked_bests_cacher.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;
using namespace cath::seq;

using std::string;

BOOST_AUTO_TEST_SUITE(masked_bests_cache_test_suite)

BOOST_AUTO_TEST_CASE(get_unmasked_regions_before_arrow_handles_simple_eg) {
	BOOST_CHECK_EQUAL_RANGES(
		get_unmasked_regions_before_arrow(
			calc_hit_vec{
				calc_hit{ arrow_before_res(  0 ), arrow_before_res( 20 ), 1.0, 0 },
				calc_hit{ arrow_before_res( 20 ), arrow_before_res( 30 ), 1.0, 1 },
				calc_hit{ arrow_before_res( 40 ), arrow_before_res( 50 ), 1.0, 2 },
				calc_hit{ arrow_before_res( 60 ), arrow_before_res( 80 ), 1.0, 3 },
			},
			arrow_after_res(  70 )
		),
		seq_seg_vec{
			{ arrow_before_res( 30 ), arrow_before_res( 40 ) },
			{ arrow_before_res( 50 ), arrow_before_res( 60 ) }
		}
	);
}


BOOST_AUTO_TEST_CASE(get_arrows_before_starts_of_doms_right_interspersed_with_all_of_handles_tricky_case) {
	const calc_hit_list discont_hits{
		full_hit_list{ full_hit_vec{
			full_hit( { seq_seg{ 10, 19 }, seq_seg{ 40, 49 }, }, "match_a", 1.0 ),
			full_hit( { seq_seg{ 30, 39 }, seq_seg{ 50, 59 }, }, "match_b", 1.0 ),
			full_hit( { seq_seg{  0,  9 }, seq_seg{ 60, 69 }, }, "match_c", 1.0 ),
		} },
		make_neutral_score_spec(),
		make_no_action_crh_segment_spec()
	};
	BOOST_CHECK_EQUAL_RANGES(
		get_arrows_before_starts_of_doms_right_interspersed_with_all_of(
			calc_hit_vec{
				calc_hit( { seq_seg{ 10, 19 }, seq_seg{ 40, 49 }, }, 1.0, 1 ),
				calc_hit( { seq_seg{  0,  9 }, seq_seg{ 60, 69 }, }, 1.0, 3 ),
			},
			discont_hits_index_by_start{ discont_hits },
			arrow_before_res( 20 )
		),
		res_arrow_vec{ arrow_before_res( 30 ) }
	);
}

BOOST_AUTO_TEST_SUITE_END()
