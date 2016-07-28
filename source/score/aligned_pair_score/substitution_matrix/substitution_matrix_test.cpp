/// \file
/// \brief The substitution_matrix test suite

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

#include "score/aligned_pair_score/substitution_matrix/blosum62_substitution_matrix.h"
#include "score/aligned_pair_score/substitution_matrix/identity_substitution_matrix.h"
#include "score/aligned_pair_score/substitution_matrix/substitution_matrix.h"
#include "structure/protein/amino_acid.h"
#include "structure/structure_type_aliases.h"

using namespace cath;
using namespace cath::score;

namespace cath {
	namespace test {

		/// \brief The substitution_matrix_test_suite_fixture to assist in testing substitution_matrix
		struct substitution_matrix_test_suite_fixture {
		protected:
			~substitution_matrix_test_suite_fixture() noexcept = default;

		public:
			/// \brief A vector of the normal amino acids (including extras like 'B', but not 'X')
			const amino_acid_vec normal_amino_acids{ make_amino_acids_of_chars(
				{ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'Z' }
			) };

			/// \brief The amino acid entry for letter 'X'
			const amino_acid     amino_acid_x{ amino_acid( 'X' ) };
		};

	}
}

/// \brief Suite of tests to unit test substitution_matrix
BOOST_FIXTURE_TEST_SUITE(substitution_matrix_test_suite, cath::test::substitution_matrix_test_suite_fixture)

/// \brief Check that the identity matrix behaves as expected
BOOST_AUTO_TEST_CASE(identity) {
	const substitution_matrix identity_matrix = make_subs_matrix_identity();
	BOOST_CHECK_EQUAL( identity_matrix.get_name(),                              "identity" );
	BOOST_CHECK_EQUAL( identity_matrix.get_highest_score(),                     1          );
	BOOST_CHECK_EQUAL( identity_matrix.get_score( amino_acid_x, amino_acid_x ), 0          );
	for (const amino_acid &amino_acid_a : normal_amino_acids) {
		BOOST_CHECK_EQUAL( identity_matrix.get_score( amino_acid_a, amino_acid_x ), 0 );
		BOOST_CHECK_EQUAL( identity_matrix.get_score( amino_acid_x, amino_acid_a ), 0 );
		for (const amino_acid &amino_acid_b : normal_amino_acids) {
			const score_type correct_score = ( amino_acid_a == amino_acid_b ) ? 1 : 0;
			BOOST_CHECK_EQUAL( identity_matrix.get_score( amino_acid_a, amino_acid_b ), correct_score );
			BOOST_CHECK_EQUAL( identity_matrix.get_score( amino_acid_b, amino_acid_a ), correct_score );
		}
	}
}

/// \brief Perform some spot checks on the BLOSUM62 matrix
BOOST_AUTO_TEST_CASE(blosum62) {
	const substitution_matrix blosum62_matrix = make_subs_matrix_blosum62();
	BOOST_CHECK_EQUAL( blosum62_matrix.get_name(),                                        "blosum62" );
	BOOST_CHECK_EQUAL( blosum62_matrix.get_highest_score(),                               11         );
	BOOST_CHECK_EQUAL( blosum62_matrix.get_score( amino_acid( 'A' ), amino_acid( 'S' ) ),  1         );
	BOOST_CHECK_EQUAL( blosum62_matrix.get_score( amino_acid( 'S' ), amino_acid( 'A' ) ),  1         );
	BOOST_CHECK_EQUAL( blosum62_matrix.get_score( amino_acid( 'H' ), amino_acid( 'V' ) ), -3         );
	BOOST_CHECK_EQUAL( blosum62_matrix.get_score( amino_acid( 'V' ), amino_acid( 'H' ) ), -3         );
}

BOOST_AUTO_TEST_SUITE_END()
