/// \file
/// \brief The domain_cluster_ids test suite

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
#include <boost/test/unit_test.hpp>

#include "cath/cluster/domain_cluster_ids.hpp"
#include "cath/common/boost_addenda/range/front.hpp"

using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::seq;

using ::boost::none;

BOOST_AUTO_TEST_SUITE(domain_cluster_ids_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	domain_cluster_ids the_ids;

	BOOST_CHECK      (        the_ids.empty()               );
	BOOST_CHECK_EQUAL(        the_ids.size(),          0    );

	the_ids.emplace_back(
		none,
		5ul
	);

	BOOST_CHECK      (      ! the_ids.empty()               );
	BOOST_CHECK_EQUAL(        the_ids.size(),          1    );

	BOOST_CHECK_EQUAL(        the_ids[ 0 ].segments,   none );
	BOOST_CHECK_EQUAL(        the_ids[ 0 ].cluster_id, 5    );

	BOOST_CHECK_EQUAL( front( the_ids ).segments,      none );
	BOOST_CHECK_EQUAL( front( the_ids ).cluster_id,    5    );
}

BOOST_AUTO_TEST_SUITE_END()
