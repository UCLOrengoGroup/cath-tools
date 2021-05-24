/// \file
/// \brief The options_block test suite

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

#include <boost/core/demangle.hpp>
#include <boost/mpl/list.hpp>

#include "cath/alignment/options_block/alignment_input_options_block.hpp"
#include "cath/options/options_block/detail_help_options_block.hpp"
#include "cath/options/options_block/options_block_tester.hpp"
#include "cath/options/options_block/pdb_input_options_block.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace ::cath;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::core::demangle;
using ::boost::program_options::options_description;
using ::boost::program_options::unknown_option;
using ::boost::program_options::variables_map;

BOOST_TEST_DONT_PRINT_LOG_VALUE(type_info)

/// \brief A type alias for a list of all the different types of options_blocks, so they can all be tested
using all_options_block_types = boost::mpl::list<alignment_input_options_block,
                                                 detail_help_options_block,
                                                 pdb_input_options_block>;

namespace {

	/// \brief The options_block_test_suite_fixture to assist in testing options_block
	struct options_block_test_suite_fixture : protected options_block_tester {
	protected:
		~options_block_test_suite_fixture() noexcept = default;

	public:
		template <typename OB>
		OB construct_options_block_for_testing();
	};


	/// \brief Template method for constructing options_blocks
	///
	/// This is here to allow template specialisation for any options_block class that cannot be default constructed
	template <typename OB>
	OB options_block_test_suite_fixture::construct_options_block_for_testing() {
		return OB();
	}

	/// \brief Specialisation of template method to construct a detail_help_options_block
	///
	/// This is needed because detail_help_options_block cannot be default constructed
	template <>
	detail_help_options_block options_block_test_suite_fixture::construct_options_block_for_testing<detail_help_options_block>() {
		return detail_help_options_block( TEST_DESC_AND_HELP_OF_OPTION_NAME() );
	}
} // namespace


BOOST_FIXTURE_TEST_SUITE(options_block_test_suite, options_block_test_suite_fixture)

/// \brief Check that each type of options_block can be default constructed without throwing
BOOST_AUTO_TEST_CASE_TEMPLATE(ctor_does_not_throw, options_block_type, all_options_block_types) {
	BOOST_CHECK_NO_THROW_DIAG(options_block_type a((construct_options_block_for_testing<options_block_type>())));
}

/// \brief Check that each type of options_block can be successfully cloned to the correct dynamic type
BOOST_AUTO_TEST_CASE_TEMPLATE(clone_works, options_block_type, all_options_block_types) {
	const options_block_type the_options_block(construct_options_block_for_testing<options_block_type>());
	const unique_ptr<options_block> options_block_new_clone_ptr = the_options_block.clone();
	const auto &options_block_ref = *options_block_new_clone_ptr;
	BOOST_CHECK_EQUAL(           typeid( options_block_ref ),                    typeid( the_options_block )          );
	BOOST_CHECK_EQUAL( demangle( typeid( options_block_ref ).name() ), demangle( typeid( the_options_block ).name() ) );
}

/// \brief Check that each type of options_block will return an options_description from get_options_description()
BOOST_AUTO_TEST_CASE_TEMPLATE(get_options_description_works, options_block_type, all_options_block_types) {
	options_block_type the_options_block(construct_options_block_for_testing<options_block_type>());
	BOOST_CHECK_NO_THROW_DIAG( const options_description desc = the_options_block.get_visible_options_description( 100 ) );
	BOOST_CHECK_NO_THROW_DIAG( const options_description desc = the_options_block.get_hidden_options_description ( 100 ) );
}

/// \brief Check that each type of options_block will return an options_description from get_options_description()
BOOST_AUTO_TEST_CASE_TEMPLATE(invalid_string_works, options_block_type, all_options_block_types) {
	const options_block_type the_options_block(construct_options_block_for_testing<options_block_type>());
	variables_map vm;
	BOOST_CHECK_NO_THROW_DIAG( const str_opt desc = the_options_block.invalid_string( vm ) );
}

/// \brief Check that the parsing throws on an unrecognised option
BOOST_AUTO_TEST_CASE_TEMPLATE(throws_on_unrecognised_option, options_block_type, all_options_block_types) {
	BOOST_CHECK_THROW(
		parse_into_options_block_copy( construct_options_block_for_testing<options_block_type>(), { string( UNKNOWN_OPT ) } ),
		unknown_option
	);
}

BOOST_AUTO_TEST_SUITE_END()

