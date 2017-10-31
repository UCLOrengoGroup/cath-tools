/// \file
/// \brief The alignment_split_mapping test suite

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

#include "alignment/align_type_aliases.hpp"
#include "alignment/refiner/detail/alignment_split_mapping.hpp"
#include "common/size_t_literal.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "test/boost_addenda/boost_check_no_throw_diag.hpp"

#include <iostream>

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;
using namespace std;

using boost::none;

namespace cath {
	namespace test {

		/// \brief The alignment_split_mapping_test_suite_fixture to assist in testing alignment_split_mapping
		struct alignment_split_mapping_test_suite_fixture {
		protected:
			~alignment_split_mapping_test_suite_fixture() noexcept = default;

		public:
			/// Source alignment:
			///     -1245-
			///     01--24
			///     012345
			///     01--34
			const aln_posn_opt_vec_vec SRC_ALN_DATA = {
				{ none, 1,    2,    4,    5,    none },
				{ 0_z,  1,    none, none, 2,    4    },
				{ 0_z,  1,    2,    3,    4,    5    },
				{ 0_z,  1,    none, none, 3,    4    }
			};
			const alignment SRC_ALN{ SRC_ALN_DATA };

			/// \brief TODOCUMENT
			const size_vec SRC_CORR_LENGTHS = { 6, 5, 7, 5 };

			/// Mapping a (from 0 and 2):
			///     -012345--
			///     0-12-3456
			const size_set SPLIT_SET_A = { 0, 2 };

			/// Mapping b (from 1 and 3):
			///     01-234
			///     0123-4
			const size_set SPLIT_SET_B = { 1, 3 };

			/// Source inter-mapping alignment:
			///     -01.2345--
			///     0-1.2-3456
			///     ||  a   ||
			///     vv      vv
			///     012-345678
			///     -012345---
			///     ^^      ^^
			///     ||  b   ||
			///     .01-234...
			///     .0123-4...
			const aln_posn_opt_vec_vec INTER_ALN_DATA = {
				{ 0_z,  1,    2,    none, 3,    4,    5,    6,    7,    8    },
				{ none, 0_z,  1,    2,    3,    4,    5,    none, none, none }
			};
			const alignment INTER_ALN{ INTER_ALN_DATA };

			/// Final alignment:
			///     -01-2345-
			///     -01-234--
			///     0-1-2-345
			///     -0123-4--
			const aln_posn_opt_vec_vec FINAL_ALN_DATA = {
				{ none, 0_z,  1, none, 2,    3,    4,    5,    none, none       },
				{ none, 0_z,  1, none, 2,    3,    4,    none, none, none       },
				{ 0_z,  none, 1, none, 2,    none, 3,    4,    5,    6          },
				{ none, 0_z,  1,       2,    3,    none, 4,    none, none, none }
			};
			const alignment FINAL_ALN{ FINAL_ALN_DATA };

		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(alignment_split_mapping_test_suite, cath::test::alignment_split_mapping_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ctor_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( alignment_split_mapping( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS ) );
	BOOST_CHECK_NO_THROW_DIAG( alignment_split_mapping( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(orig_aln_length_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );
	BOOST_CHECK_EQUAL( test_mapping_a.orig_aln_length(), 6_z );
	BOOST_CHECK_EQUAL( test_mapping_b.orig_aln_length(), 6_z );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(orig_aln_num_entries_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );
	BOOST_CHECK_EQUAL( test_mapping_a.orig_aln_num_entries(), 4_z );
	BOOST_CHECK_EQUAL( test_mapping_b.orig_aln_num_entries(), 4_z );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(inserted_entries_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );
	BOOST_CHECK_EQUAL( test_mapping_a.inserted_entries(), true );
	BOOST_CHECK_EQUAL( test_mapping_b.inserted_entries(), true );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(length_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );
	BOOST_CHECK_EQUAL( test_mapping_a.length(), 9_z );
	BOOST_CHECK_EQUAL( test_mapping_b.length(), 6_z );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(num_entries_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );
	BOOST_CHECK_EQUAL( test_mapping_a.num_entries(), 2_z );
	BOOST_CHECK_EQUAL( test_mapping_b.num_entries(), 2_z );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(index_of_orig_aln_index_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );

	BOOST_CHECK_EQUAL( test_mapping_a.index_of_orig_aln_index( 0 ), size_opt( 0_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_orig_aln_index( 1 ), size_opt( 2_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_orig_aln_index( 2 ), size_opt( 3_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_orig_aln_index( 3 ), size_opt( 5_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_orig_aln_index( 4 ), size_opt( 6_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_orig_aln_index( 5 ), size_opt( 7_z   ) );

	BOOST_CHECK_EQUAL( test_mapping_b.index_of_orig_aln_index( 0 ), size_opt( 0_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_orig_aln_index( 1 ), size_opt( 1_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_orig_aln_index( 2 ), size_opt( none ) );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_orig_aln_index( 3 ), size_opt( none ) );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_orig_aln_index( 4 ), size_opt( 3_z   ) );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_orig_aln_index( 5 ), size_opt( 5_z   ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(entry_of_orig_aln_entry_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );

	BOOST_CHECK_EQUAL( test_mapping_a.entry_of_orig_aln_entry( 0 ), size_opt( 0_z   )  );
	BOOST_CHECK_EQUAL( test_mapping_a.entry_of_orig_aln_entry( 1 ), size_opt( none )  );
	BOOST_CHECK_EQUAL( test_mapping_a.entry_of_orig_aln_entry( 2 ), size_opt( 1_z   )  );
	BOOST_CHECK_EQUAL( test_mapping_a.entry_of_orig_aln_entry( 3 ), size_opt( none )  );

	BOOST_CHECK_EQUAL( test_mapping_b.entry_of_orig_aln_entry( 0 ), size_opt( none )  );
	BOOST_CHECK_EQUAL( test_mapping_b.entry_of_orig_aln_entry( 1 ), size_opt( 0_z   )  );
	BOOST_CHECK_EQUAL( test_mapping_b.entry_of_orig_aln_entry( 2 ), size_opt( none )  );
	BOOST_CHECK_EQUAL( test_mapping_b.entry_of_orig_aln_entry( 3 ), size_opt( 1_z   )  );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(orig_aln_entry_of_entry_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );

	BOOST_CHECK_EQUAL( test_mapping_a.orig_aln_entry_of_entry( 0 ), 0_z );
	BOOST_CHECK_EQUAL( test_mapping_a.orig_aln_entry_of_entry( 1 ), 2_z );

	BOOST_CHECK_EQUAL( test_mapping_b.orig_aln_entry_of_entry( 0 ), 1_z );
	BOOST_CHECK_EQUAL( test_mapping_b.orig_aln_entry_of_entry( 1 ), 3_z );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(position_of_entry_of_index_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );

	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 0 ), none               );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 1 ), aln_posn_type( 0 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 2 ), aln_posn_type( 1 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 3 ), aln_posn_type( 2 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 4 ), aln_posn_type( 3 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 5 ), aln_posn_type( 4 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 6 ), aln_posn_type( 5 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 7 ), none               );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 0, 8 ), none               );

	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 0 ), aln_posn_type( 0 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 1 ), none               );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 2 ), aln_posn_type( 1 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 3 ), aln_posn_type( 2 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 4 ), none               );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 5 ), aln_posn_type( 3 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 6 ), aln_posn_type( 4 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 7 ), aln_posn_type( 5 ) );
	BOOST_CHECK_EQUAL( test_mapping_a.position_of_entry_of_index( 1, 8 ), aln_posn_type( 6 ) );

	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 0, 0 ), aln_posn_type( 0 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 0, 1 ), aln_posn_type( 1 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 0, 2 ), none               );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 0, 3 ), aln_posn_type( 2 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 0, 4 ), aln_posn_type( 3 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 0, 5 ), aln_posn_type( 4 ) );

	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 1, 0 ), aln_posn_type( 0 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 1, 1 ), aln_posn_type( 1 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 1, 2 ), aln_posn_type( 2 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 1, 3 ), aln_posn_type( 3 ) );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 1, 4 ), none               );
	BOOST_CHECK_EQUAL( test_mapping_b.position_of_entry_of_index( 1, 5 ), aln_posn_type( 4 ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(index_of_protein_index_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );

	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 0, 0 ), 1_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 0, 1 ), 2_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 0, 2 ), 3_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 0, 3 ), 4_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 0, 4 ), 5_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 0, 5 ), 6_z );

	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 1, 0 ), 0_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 1, 1 ), 2_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 1, 2 ), 3_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 1, 3 ), 5_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 1, 4 ), 6_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 1, 5 ), 7_z );
	BOOST_CHECK_EQUAL( test_mapping_a.index_of_protein_index( 1, 6 ), 8_z );

	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 0, 0 ), 0_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 0, 1 ), 1_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 0, 2 ), 3_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 0, 3 ), 4_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 0, 4 ), 5_z );

	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 1, 0 ), 0_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 1, 1 ), 1_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 1, 2 ), 2_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 1, 3 ), 3_z );
	BOOST_CHECK_EQUAL( test_mapping_b.index_of_protein_index( 1, 4 ), 5_z );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(build_from_mappings_works) {
	const alignment_split_mapping test_mapping_a( SRC_ALN, SPLIT_SET_A, SRC_CORR_LENGTHS );
	const alignment_split_mapping test_mapping_b( SRC_ALN, SPLIT_SET_B, SRC_CORR_LENGTHS );
	BOOST_CHECK_EQUAL( build_alignment( INTER_ALN, test_mapping_a, test_mapping_b ), FINAL_ALN );
}

BOOST_AUTO_TEST_SUITE_END()
