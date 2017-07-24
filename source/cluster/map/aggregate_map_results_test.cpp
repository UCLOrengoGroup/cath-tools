/// \map
/// \brief The aggregate_map_results test suite

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

#include "cluster/map/aggregate_map_results.hpp"
#include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"

using namespace cath::clust;

BOOST_AUTO_TEST_SUITE(aggregate_map_results_test_suite)

BOOST_AUTO_TEST_CASE(default_ctor_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( aggregate_map_results a{} );
}

BOOST_AUTO_TEST_SUITE_END()
