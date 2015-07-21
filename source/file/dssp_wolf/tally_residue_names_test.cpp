/// \file
/// \brief The tally_residue_names test suite

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/boost_check_no_throw_diag.h"
#include "common/pair_insertion_operator.h"
#include "exception/invalid_argument_exception.h"
#include "file/dssp_wolf/tally_residue_names.h"
#include "structure/residue_name.h"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

//namespace std {
//	/// \brief Naughty addition of an insertion operator into std:: to get Boost.Test to output size_size_pairs
//	ostream & operator<<(ostream              &arg_os,            ///< ostream to which to output the size_size_pair
//	                     const size_size_pair &arg_size_size_pair ///< size_size_pair to output
//	                     ) {
//		arg_os << "pair[" << arg_size_size_pair.first  << ", ";
//		arg_os            << arg_size_size_pair.second << "]";
//		return arg_os;
//	}
//}

namespace cath {
	namespace test {

		/// \brief The tally_residue_names_test_suite_fixture to assist in testing tally_residue_names
		struct tally_residue_names_test_suite_fixture {
		protected:
			~tally_residue_names_test_suite_fixture() noexcept = default;

		public:
			const residue_name_vec   STANDARD_RES_NAMES = { residue_name( 2 ),
			                                                residue_name( 3 ),
			                                                residue_name( 4 ) };
			const residue_name_vec   MISMATCH_RES_NAMES = { residue_name( 5 ),
			                                                residue_name( 6 ),
			                                                residue_name( 7 ) };

			const size_size_pair_vec STANDARD_MATCH     = { { 0, 0 },
			                                                { 1, 1 },
			                                                { 2, 2 } };
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(tally_residue_names_test_suite, cath::test::tally_residue_names_test_suite_fixture)


/// \brief Check that tally_residue_names() works as expected on simple, completely matching lists
BOOST_AUTO_TEST_CASE(works_on_simple_match) {
	const size_size_pair_vec got = tally_residue_names( STANDARD_RES_NAMES, STANDARD_RES_NAMES, false );
	BOOST_CHECK_EQUAL_RANGES( STANDARD_MATCH, got );
}



/// \brief Check that tally_residue_names() throws on a PDB list that contains an empty string
BOOST_AUTO_TEST_CASE(throws_on_pdb_list_contains_empty) {
	BOOST_CHECK_THROW(tally_residue_names(
		{ residue_name( 2 ),
		  residue_name(   ),
		  residue_name( 4 ) },
		STANDARD_RES_NAMES,
		false
	), invalid_argument_exception);
}

/// \brief Check that tally_residue_names() throws on a PDB list that contains duplicates
BOOST_AUTO_TEST_CASE(throws_on_pdb_list_contains_duplicates) {
	BOOST_CHECK_THROW(tally_residue_names(
		{ residue_name( 2 ),
		  residue_name( 3 ),
		  residue_name( 2 ),
		  residue_name( 4 ) },
		STANDARD_RES_NAMES,
		false
	), invalid_argument_exception);
}

/// \brief Check that tally_residue_names() throws on a DSSP list that contains consecutive duplicate empty strings
BOOST_AUTO_TEST_CASE(throws_on_DSSP_list_contains_consecutive_duplicate_empty) {
	BOOST_CHECK_THROW(tally_residue_names(
		STANDARD_RES_NAMES,
		{ residue_name( 2 ),
		  residue_name(   ),
		  residue_name(   ),
		  residue_name( 4 ) },
		false
	), invalid_argument_exception);
}



/// \brief Check that tally_residue_names() throws on completely mismatching lists, with arg_permit_breaks_without_null_residues off or on
BOOST_AUTO_TEST_CASE(throws_on_bad_mismatch) {
	for (const bool &arg_permit : { false, true }) {
		BOOST_CHECK_THROW(tally_residue_names(STANDARD_RES_NAMES, MISMATCH_RES_NAMES, arg_permit), invalid_argument_exception);
	}
}



/// \brief Check that tally_residue_names() with arg_permit_breaks_without_null_residues off throws on a DSSP gap at the start
BOOST_AUTO_TEST_CASE(throws_if_not_permit_dssp_gap_at_start) {
	BOOST_CHECK_THROW(tally_residue_names(
		STANDARD_RES_NAMES,
		{ residue_name( 3 ),
		  residue_name( 4 ) },
		false,
		false
	), invalid_argument_exception);

	BOOST_CHECK_NO_THROW_DIAG(tally_residue_names(
		STANDARD_RES_NAMES,
		{ residue_name( 3 ),
		  residue_name( 4 ) },
		false,
		true
	));
}

/// \brief Check that tally_residue_names() with arg_permit_breaks_without_null_residues off throws on a DSSP gap at the middle
BOOST_AUTO_TEST_CASE(throw_if_not_permit_dssp_gap_at_middle) {
	BOOST_CHECK_THROW(tally_residue_names(
		STANDARD_RES_NAMES,
		{ residue_name( 2 ),
		  residue_name( 4 ) },
		false
	), invalid_argument_exception);
}

/// \brief Check that tally_residue_names() with arg_permit_tail_break_without_null_residue off throws on a DSSP gap at the end
BOOST_AUTO_TEST_CASE(throw_if_not_permit_dssp_gap_at_end) {
	BOOST_CHECK_THROW(tally_residue_names(
		STANDARD_RES_NAMES,
		{ residue_name( 2 ),
		  residue_name( 3 ) },
		false,
		false
	), invalid_argument_exception);
}



/// \brief Check that tally_residue_names() with arg_permit_breaks_without_null_residues on allows a gap at the start
BOOST_AUTO_TEST_CASE(works_if_permit_dssp_gap_at_start) {
	const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, {      "3", "4" }, true);
	const size_size_pair_vec expected = { { 1, 0 },
	                                      { 2, 1 } };
	BOOST_CHECK_EQUAL_RANGES( expected, got );
}

/// \brief Check that tally_residue_names() with arg_permit_breaks_without_null_residues on allows a gap at the middle
BOOST_AUTO_TEST_CASE(works_if_permit_dssp_gap_at_middle) {
	const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, { "2",      "4" }, true);
	const size_size_pair_vec expected = { { 0, 0 },
	                                      { 2, 1 } };
	BOOST_CHECK_EQUAL_RANGES( expected, got );
}

/// \brief Check that tally_residue_names() with arg_permit_breaks_without_null_residues on allows a gap at the end
BOOST_AUTO_TEST_CASE(works_if_permit_dssp_gap_at_end) {
	const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, { "2", "3"      }, true);
	const size_size_pair_vec expected = { { 0, 0 },
	                                      { 1, 1 } };
	BOOST_CHECK_EQUAL_RANGES( expected, got );
}



/// \brief Check that tally_residue_names() works with a DSSP null matching a residue at the start, with arg_permit_breaks_without_null_residues off or on
BOOST_AUTO_TEST_CASE(works_with_null_match_at_start) {
	for (const bool &arg_permit : { false, true }) {
		const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, {  "", "3", "4" }, arg_permit);
		const size_size_pair_vec expected = { { 1, 1 },
		                                      { 2, 2 } };
		BOOST_CHECK_EQUAL_RANGES( expected, got );
	}
}

/// \brief Check that tally_residue_names() works with a DSSP null matching a residue at the middle, with arg_permit_breaks_without_null_residues off or on
BOOST_AUTO_TEST_CASE(works_with_null_match_at_middle) {
	for (const bool &arg_permit : { false, true }) {
		const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, { "2",  "", "4" }, arg_permit);
		const size_size_pair_vec expected = { { 0, 0 },
		                                      { 2, 2 } };
		BOOST_CHECK_EQUAL_RANGES( expected, got );
	}
}

/// \brief Check that tally_residue_names() works with a DSSP null matching a residue at the end, with arg_permit_breaks_without_null_residues off or on
BOOST_AUTO_TEST_CASE(works_with_null_match_at_end) {
	for (const bool &arg_permit : { false, true }) {
		const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, { "2", "3",  "" }, arg_permit);
		const size_size_pair_vec expected = { { 0, 0 },
		                                      { 1, 1 } };
		BOOST_CHECK_EQUAL_RANGES( expected, got );
	}
}



/// \brief Check that tally_residue_names() works with a DSSP null inserted at the start, with arg_permit_breaks_without_null_residues off or on
BOOST_AUTO_TEST_CASE(works_with_null_insert_at_start) {
	for (const bool &arg_permit : { false, true }) {
		const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, {  "", "2", "3", "4" }, arg_permit);
		const size_size_pair_vec expected = { { 0, 1 },
		                                      { 1, 2 },
		                                      { 2, 3 } };
		BOOST_CHECK_EQUAL_RANGES( expected, got );
	}
}

/// \brief Check that tally_residue_names() works with a DSSP null inserted at the middle, with arg_permit_breaks_without_null_residues off or on
BOOST_AUTO_TEST_CASE(works_with_null_insert_at_middle) {
	for (const bool &arg_permit : { false, true }) {
		const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, { "2",  "", "3", "4" }, arg_permit);
		const size_size_pair_vec expected = { { 0, 0 },
		                                      { 1, 2 },
		                                      { 2, 3 } };
		BOOST_CHECK_EQUAL_RANGES( expected, got );
	}
}

/// \brief Check that tally_residue_names() works with a DSSP null inserted at the end, with arg_permit_breaks_without_null_residues off or on
BOOST_AUTO_TEST_CASE(works_with_null_insert_at_end) {
	for (const bool &arg_permit : { false, true }) {
		const size_size_pair_vec got      = tally_residue_names_str(STANDARD_RES_NAMES, { "2", "3", "4",  "" }, arg_permit);
		const size_size_pair_vec expected = { { 0, 0 },
		                                      { 1, 1 },
		                                      { 2, 2 } };
		BOOST_CHECK_EQUAL_RANGES( expected, got );
	}
}



BOOST_AUTO_TEST_SUITE_END()
