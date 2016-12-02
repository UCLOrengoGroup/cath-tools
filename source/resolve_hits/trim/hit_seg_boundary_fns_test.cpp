/// \file
/// \brief The hit_seg_boundary_fns test suite

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

#include <boost/optional/optional_io.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "resolve_hits/hit_seg.hpp"
#include "resolve_hits/trim/hit_seg_boundary_fns.hpp"
#include "resolve_hits/trim/trim_spec.hpp"

using namespace cath::rslv;

BOOST_AUTO_TEST_SUITE(hit_seg_boundary_fns_test_suite)

BOOST_AUTO_TEST_CASE(handles_example) {
	constexpr auto seg_a = hit_seg_of_res_idcs( 138, 168 );
	constexpr auto seg_b = hit_seg_of_res_idcs(  67, 295 );
	constexpr auto trim  = trim_spec{ 200, 190 };
	BOOST_CHECK_EQUAL( calc_resolved_boundary( seg_a, seg_b, trim ), make_optional( arrow_before_res( 156 ) ) );
}

BOOST_AUTO_TEST_SUITE_END()
