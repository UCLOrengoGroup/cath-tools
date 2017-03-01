/// \file
/// \brief The domain test suite

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

#include "chopping/chopping_type_aliases.hpp"
#include "chopping/domain/domain.hpp"
#include "chopping/region/region.hpp"
#include "common/test_tools.hpp"

using namespace cath::chop;
using namespace cath::common::test;

BOOST_AUTO_TEST_SUITE(domain_test_suite)

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string( domain{ {                                                                          }         } ), "domain[]"                                                                                                  );
	BOOST_CHECK_EQUAL( to_string( domain{ { make_simple_region(      121, 232 )                                      }, "man"  } ), "domain[name:man, region[ start_idx:121, stop_idx:232 ]]"                                                   );
	BOOST_CHECK_EQUAL( to_string( domain{ { make_simple_region(      121, 232 )                                      }, "road" } ), "domain[name:road, region[ start_idx:121, stop_idx:232 ]]"                                                  );
	BOOST_CHECK_EQUAL( to_string( domain{ { make_simple_region(      121, 232 )                                      }         } ), "domain[region[ start_idx:121, stop_idx:232 ]]"                                                             );
	BOOST_CHECK_EQUAL( to_string( domain{ { make_simple_region(      121, 232 ), make_simple_region(      234, 434 ) }         } ), "domain[region[ start_idx:121, stop_idx:232 ],region[ start_idx:234, stop_idx:434 ]]"                       );
	BOOST_CHECK_EQUAL( to_string( domain{ { make_simple_region( 'K', 121, 232 )                                      }         } ), "domain[region[ chain:K, start_name:121, stop_name:232 ]]"                                                  );
	BOOST_CHECK_EQUAL( to_string( domain{ { make_simple_region( 'K', 121, 232 ), make_simple_region( 'K', 234, 434 ) }         } ), "domain[region[ chain:K, start_name:121, stop_name:232 ],region[ chain:K, start_name:234, stop_name:434 ]]" );
}

BOOST_AUTO_TEST_CASE(equality_works) {
	check_equality_operators_on_diff_vals_range( domain_vec{
		domain{ {                                                                          }             },
		domain{ {                                                                          }, "archaea"  },
		domain{ {                                                                          }, "bacteria" },
		domain{ { make_simple_region(      121, 232 )                                      }, "archaea"  },
		domain{ { make_simple_region(      121, 232 )                                      }, "bacteria" },
		domain{ { make_simple_region(      121, 232 )                                      }             },
		domain{ { make_simple_region(      121, 232 ), make_simple_region(      121, 232 ) }             },
		domain{ { make_simple_region( 'K', 121, 232 )                                      }             },
		domain{ { make_simple_region( 'K', 121, 232 ), make_simple_region( 'K', 121, 232 ) }             }
	} );
}

BOOST_AUTO_TEST_SUITE_END()
