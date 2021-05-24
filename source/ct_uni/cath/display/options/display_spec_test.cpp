/// \file
/// \brief The display_spec test suite

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

#include "display_spec.hpp"

#include <boost/test/unit_test.hpp>

#include "cath/display/display_colourer/display_colourer.hpp"
#include "cath/display/options/display_spec.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace ::cath;

BOOST_AUTO_TEST_SUITE(display_spec_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	const display_spec a_display_spec{};
	BOOST_CHECK_NO_THROW_DIAG( get_display_colourer( a_display_spec ) );
}

BOOST_AUTO_TEST_SUITE_END()
