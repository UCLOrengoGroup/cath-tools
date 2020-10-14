/// \file
/// \brief The superfamily_of_domain test suite

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

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/score/homcheck_tools/superfamily_of_domain.hpp"

using namespace ::cath::common;
using namespace ::cath::homcheck;
using namespace ::cath::homcheck::detail;
using namespace ::std;

namespace cath {
	namespace test {

		/// \brief The superfamily_of_domain_test_suite_fixture to assist in testing superfamily_of_domain
		struct superfamily_of_domain_test_suite_fixture {
		protected:
			~superfamily_of_domain_test_suite_fixture() noexcept = default;

			/// \brief Some example superfamily_of_domain text for testing parsing
			const string superfamily_of_domain_text = R"(1qinA00 3.10.180.10
1vkeD00 1.20.1290.10
1yjwH00 3.90.1170.10
2nx5G00 2.60.40.10
2x2jC01 2.60.40.1760
3i3eB01 2.60.120.260
3s27A04 3.40.50.2000
3x1iA00 1.10.565.10
4n9f400 3.10.20.90
6cscA01 1.10.580.10)";

			/// \brief The superfamily of domain data implied by the above text
			const str_str_map superfamily_of_domain_map = {
				{ "1qinA00", "3.10.180.10"  },
				{ "1vkeD00", "1.20.1290.10" },
				{ "1yjwH00", "3.90.1170.10" },
				{ "2nx5G00", "2.60.40.10"   },
				{ "2x2jC01", "2.60.40.1760" },
				{ "3i3eB01", "2.60.120.260" },
				{ "3s27A04", "3.40.50.2000" },
				{ "3x1iA00", "1.10.565.10"  },
				{ "4n9f400", "3.10.20.90"   },
				{ "6cscA01", "1.10.580.10"  }
			};

		};
	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(superfamily_of_domain_test_suite, cath::test::superfamily_of_domain_test_suite_fixture)

BOOST_AUTO_TEST_CASE(is_valid_cath_node_id_works) {
	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "1"           ), true );
	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "1.10"        ), true );
	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "1.10.8"      ), true );
	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "1.10.8.10"   ), true );

	// Currently rejects S, O, L, I or D
	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "1.10.8.10.1" ), false );

	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "1.1a0.8"     ), false );

	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "1a"          ), false );
	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "a1"          ), false );
	BOOST_CHECK_EQUAL( is_valid_cath_node_id{}( "a"           ), false );
}

// this should probably be split up a bit
BOOST_AUTO_TEST_CASE(parses_correctly) {
	const auto sf_of_dom = parse_superfamily_of_domain( superfamily_of_domain_text );

	BOOST_CHECK_EQUAL(   sf_of_dom.size(), 10 );
	BOOST_CHECK      ( ! sf_of_dom.has_superfamily_of_domain( "1cukA01" ) );
	BOOST_CHECK_THROW(   sf_of_dom.get_superfamily_of_domain( "1cukA01" ), invalid_argument_exception );
	for (const auto &domain_superfamily_pair : superfamily_of_domain_map) {
		const string &domain_id      = domain_superfamily_pair.first;
		const string &superfamily_id = domain_superfamily_pair.second;
		BOOST_CHECK      ( sf_of_dom.has_superfamily_of_domain( domain_id ) );
		BOOST_CHECK_EQUAL( sf_of_dom.get_superfamily_of_domain( domain_id ), superfamily_id );
	}
}


BOOST_AUTO_TEST_CASE(parses_single_line) {
	BOOST_CHECK_EQUAL( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.10" } ).size(), 1 );
}


BOOST_AUTO_TEST_CASE(parses_two_lines) {
	BOOST_CHECK_EQUAL( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.10\n1vkeD00 1.20.1290.10" } ).size(), 2 );
}


BOOST_AUTO_TEST_CASE(parse_rejects_repeated_domain_ids) {
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.10\n1qinA00 1.20.1290.10" } ), runtime_error_exception );
}


BOOST_AUTO_TEST_CASE(parse_rejects_wrong_num_fields) {
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00"                         } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.10 3.10.180.10" } ), runtime_error_exception );
}


BOOST_AUTO_TEST_CASE(rejects_non_numeric_superfamily_field) {
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 a.10.180.10" } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.1b.180.10" } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.1c0.10" } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.d0" } ), runtime_error_exception );
}


BOOST_AUTO_TEST_CASE(rejects_wrong_num_superfamily_fields) {
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3"                 } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10"              } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.180"          } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.10.1"     } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.10.1.1"   } ), runtime_error_exception );
	BOOST_CHECK_THROW( parse_superfamily_of_domain( string{ "1qinA00 3.10.180.10.1.1.1" } ), runtime_error_exception );
}


BOOST_AUTO_TEST_CASE(correctly_adds_new_superfamily) {
	auto sf_of_dom = parse_superfamily_of_domain( superfamily_of_domain_text );

	sf_of_dom.add_domain_in_new_sf_in_fold_of_domain( "1cukA01", "2x2jC01" );

	BOOST_CHECK_EQUAL( sf_of_dom.size(),                                 11                                  );
	BOOST_CHECK      ( sf_of_dom.has_superfamily_of_domain( "1cukA01" )                                      );
	BOOST_CHECK_EQUAL( sf_of_dom.get_superfamily_of_domain( "1cukA01" ), "2.60.40.new_sf_in_fold_of_2x2jC01" );
}


BOOST_AUTO_TEST_CASE(correctly_adds_new_superfamily_of_new_superfamily) {
	auto sf_of_dom = parse_superfamily_of_domain( superfamily_of_domain_text );

	sf_of_dom.add_domain_in_new_sf_in_fold_of_domain( "1cukA01", "2x2jC01" );
	sf_of_dom.add_domain_in_new_sf_in_fold_of_domain( "1cg2A01", "1cukA01" );

	BOOST_CHECK_EQUAL( sf_of_dom.size(),                                 12                                  );
	BOOST_CHECK      ( sf_of_dom.has_superfamily_of_domain( "1cg2A01" )                                      );
	BOOST_CHECK_EQUAL( sf_of_dom.get_superfamily_of_domain( "1cg2A01" ), "2.60.40.new_sf_in_fold_of_1cukA01" );
}


BOOST_AUTO_TEST_CASE(throws_on_attempt_to_add_in_fold_of_absent_domain) {
	auto sf_of_dom = parse_superfamily_of_domain( superfamily_of_domain_text );

	BOOST_CHECK_THROW( sf_of_dom.add_domain_in_new_sf_in_fold_of_domain( "1cukA01", "1cg2A01" ), invalid_argument_exception );
}


BOOST_AUTO_TEST_CASE(throws_on_attempt_to_add_for_existing_domain) {
	auto sf_of_dom = parse_superfamily_of_domain( superfamily_of_domain_text );

	BOOST_CHECK_THROW( sf_of_dom.add_domain_in_new_sf_in_fold_of_domain( "4n9f400", "2x2jC01" ), invalid_argument_exception );
}


BOOST_AUTO_TEST_SUITE_END()
