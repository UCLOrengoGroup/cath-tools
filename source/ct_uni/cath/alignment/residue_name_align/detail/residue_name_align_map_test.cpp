/// \file


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

#include <boost/test/tools/output_test_stream.hpp>

#include "cath/alignment/residue_name_align/detail/residue_name_align_map.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

#include <vector>

using namespace ::cath;
using namespace ::cath::align::detail;
using namespace ::cath::common;
using namespace ::std;

using ::boost::test_tools::output_test_stream;

namespace {

	/// \brief The residue_name_align_map_test_suite_fixture to assist in testing residue_name_align_map
	struct residue_name_align_map_test_suite_fixture {
	protected:
		~residue_name_align_map_test_suite_fixture() noexcept = default;

		const str_vec valid_names = {
			"my_zero",
			"my_one",
			"my_two",
			"my_three",
			"my_four",
			"my_five",
		};
		const str_vec names_with_duplicate = {
			"my_zero",
			"my_one",
			"my_two",
			"my_three",
			"my_one",
			"my_four",
		};
	};

} // namespace

BOOST_FIXTURE_TEST_SUITE(residue_name_align_map_test_suite, residue_name_align_map_test_suite_fixture)

/// \brief Check that the residue_name_align_map does what you'd expect via its get_index_of_residue_name() and contains_residue_name() methods
BOOST_AUTO_TEST_CASE(indexing_works) {
	const residue_name_align_map my_map( valid_names );
	for (const size_t &name_ctr : indices( valid_names.size() ) ) {
		const string name = valid_names[name_ctr];
		BOOST_CHECK_EQUAL( name_ctr, my_map.get_index_of_residue_name_string( name ) );
		BOOST_CHECK_EQUAL( true,     my_map.contains_residue_name_string    ( name ) );
	}
}

/// \brief Check that the constructor throws if given a list containing duplicates
BOOST_AUTO_TEST_CASE(ctor_arg_with_dup_throws) {
	BOOST_CHECK_THROW( const residue_name_align_map my_map( names_with_duplicate ), invalid_argument_exception );
}

/// \brief Check that get_index_of_residue_name() throws and contains_residue_name() returns false if passed a previously unseen string
BOOST_AUTO_TEST_CASE( unseens_work ) {
	const residue_name_align_map my_map( valid_names );
	BOOST_CHECK_THROW( [[maybe_unused]] auto &&x =
	                     my_map.get_index_of_residue_name_string( "my_previously_unseen_string" ),
	                   invalid_argument_exception );
	BOOST_CHECK_EQUAL( false, my_map.contains_residue_name_string( "my_previously_unseen_string" ) );
}

/// \brief Check that the residue_name_align_map insertion operator correctly summarises the residue_name_align_map
BOOST_AUTO_TEST_CASE(insertion_operator) {
	const residue_name_align_map my_map( valid_names );
	output_test_stream my_output_test_stream;
	my_output_test_stream << my_map;
	BOOST_CHECK( my_output_test_stream.is_equal(
			"residue_name_align_map[ my_five -> 5 my_four -> 4 my_one -> 1 my_three -> 3 my_two -> 2 my_zero -> 0]"
	));
}

BOOST_AUTO_TEST_SUITE_END()
