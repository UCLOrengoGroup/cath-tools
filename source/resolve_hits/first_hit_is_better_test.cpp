/// \file
/// \brief The first_hit_is_better test suite

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

#include "resolve_hits/calc_hit_list.hpp"
#include "resolve_hits/first_hit_is_better.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"
#include "seq/seq_type_aliases.hpp"

using namespace cath::rslv;
using namespace cath::seq;

using boost::logic::indeterminate;

BOOST_AUTO_TEST_SUITE(first_hit_is_better_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	const full_hit      a             { seq_seg_vec{ { 20, 79 } }, "betty",   1.0 };
	const full_hit      b             { seq_seg_vec{ { 10, 89 } }, "camilla", 2.0 };
	const full_hit_list the_full_list { { a, b } };
	const calc_hit_list the_calc_list { calc_hit_list( the_full_list, crh_score_spec{}, crh_segment_spec{} ) };
	BOOST_CHECK( indeterminate( first_hit_is_better( the_calc_list[ 0 ], the_calc_list[ 1 ], the_full_list ) ) );
	BOOST_CHECK( indeterminate( first_hit_is_better( the_calc_list[ 1 ], the_calc_list[ 0 ], the_full_list ) ) );
}


BOOST_AUTO_TEST_SUITE_END()
