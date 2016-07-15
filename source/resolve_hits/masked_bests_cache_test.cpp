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

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "resolve_hits/masked_bests_cache.h"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using std::string;

BOOST_AUTO_TEST_SUITE(masked_bests_cache_test_suite)

BOOST_AUTO_TEST_CASE(get_unmasked_regions_before_arrow_handles_simple_eg) {
	BOOST_CHECK_EQUAL_RANGES(
		get_unmasked_regions_before_arrow(
			hit_vec{
				hit{ arrow_before_res(  0 ), arrow_before_res( 20 ), 1.0, 0 },
				hit{ arrow_before_res( 20 ), arrow_before_res( 30 ), 1.0, 1 },
				hit{ arrow_before_res( 40 ), arrow_before_res( 50 ), 1.0, 2 },
				hit{ arrow_before_res( 60 ), arrow_before_res( 80 ), 1.0, 3 },
			},
			arrow_after_res(  70 )
		),
		hit_seg_vec{
			{ arrow_before_res( 30 ), arrow_before_res( 40 ) },
			{ arrow_before_res( 50 ), arrow_before_res( 60 ) }
		}
	);
}

BOOST_AUTO_TEST_SUITE_END()
