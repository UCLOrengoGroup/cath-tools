/// \file
/// \brief The coarse_element_type test suite

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

#include <boost/lexical_cast.hpp>

#include "common/cpp14/make_unique.h"
#include "file/pdb/coarse_element_type.h"

#include <memory>

using namespace cath;
using namespace cath::file;

using boost::lexical_cast;
using std::string;

BOOST_AUTO_TEST_SUITE(coarse_element_type_test_suite)

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string           ( coarse_element_type::CARBON       ), "carbon"       );
	BOOST_CHECK_EQUAL( to_string           ( coarse_element_type::CARBON_ALPHA ), "carbon_alpha" );
	BOOST_CHECK_EQUAL( to_string           ( coarse_element_type::CARBON_BETA  ), "carbon_beta"  );
	BOOST_CHECK_EQUAL( to_string           ( coarse_element_type::NITROGEN     ), "nitrogen"     );
	BOOST_CHECK_EQUAL( to_string           ( coarse_element_type::OXYGEN       ), "oxygen"       );
	BOOST_CHECK_EQUAL( to_string           ( coarse_element_type::NON_CORE     ), "non_core"     );
}

BOOST_AUTO_TEST_CASE(insertion_operator_works) {
	BOOST_CHECK_EQUAL( lexical_cast<string>( coarse_element_type::CARBON       ), "carbon"       );
	BOOST_CHECK_EQUAL( lexical_cast<string>( coarse_element_type::CARBON_ALPHA ), "carbon_alpha" );
	BOOST_CHECK_EQUAL( lexical_cast<string>( coarse_element_type::CARBON_BETA  ), "carbon_beta"  );
	BOOST_CHECK_EQUAL( lexical_cast<string>( coarse_element_type::NITROGEN     ), "nitrogen"     );
	BOOST_CHECK_EQUAL( lexical_cast<string>( coarse_element_type::OXYGEN       ), "oxygen"       );
	BOOST_CHECK_EQUAL( lexical_cast<string>( coarse_element_type::NON_CORE     ), "non_core"     );
}

BOOST_AUTO_TEST_SUITE_END()