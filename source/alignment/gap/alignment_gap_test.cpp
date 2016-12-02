/// \file
/// \brief The alignment_gap test suite

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

#include <boost/algorithm/string/erase.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "alignment/align_type_aliases.hpp"
#include "alignment/alignment.hpp"
#include "alignment/gap/alignment_gap.hpp"
#include "alignment/gap/gap_penalty.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace boost::algorithm;
using namespace cath;
using namespace cath::align;
using namespace cath::align::gap;
using namespace cath::common;
using namespace std;

using boost::none;
using boost::numeric_cast;

namespace cath {
	namespace test {

		/// \brief The alignment_gap_test_suite_fixture to assist in testing alignment_gap
		struct alignment_gap_test_suite_fixture {
		protected:
			~alignment_gap_test_suite_fixture() noexcept = default;

		public:
			alignment make_gap_alignment_of_strings(const str_vec &);

			score_type  open_penalty    = { 100 };
			score_type  extend_penalty  = {   1 };
			gap_penalty the_gap_penalty = { open_penalty, extend_penalty };
		};

	}
}

/// \brief Build an alignment from strings of 'x's for positions and '-'s for gaps (and ignored spaces for formatting)
alignment cath::test::alignment_gap_test_suite_fixture::make_gap_alignment_of_strings(const str_vec &arg_strings ///< The strings from which to make the alignment entries
                                                                                      ) {
	aln_posn_opt_vec_vec entries;
	entries.reserve( arg_strings.size() );
	for (const string &the_string : arg_strings) {
		const string stripped_string = erase_all_copy( the_string, " " );
		aln_posn_opt_vec positions;
		positions.reserve( stripped_string.length() );
		size_t counter = 0;
		for (const char &character : stripped_string) {
			if ( character == 'x' ) {
				positions.push_back( counter );
				++counter;
			}
			else if ( character == '-' ) {
				positions.push_back( none );
			}
			else {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(
					"Unable to recognise character "
					+ string( 1, character )
					+ " in gap alignment string"
				));
			}
		}
		entries.push_back( positions );
	}
	return alignment( entries );
}

/// \brief Test the code the calculate gap penalties for alignments
BOOST_FIXTURE_TEST_SUITE(alignment_gap_test_suite, cath::test::alignment_gap_test_suite_fixture)

/// \brief Check that gap_penalty_value_of_alignment() works for a simple alignment with a mix of issues
BOOST_AUTO_TEST_CASE(pair_bounce) {
	BOOST_CHECK_EQUAL(
		  2 * numeric_cast<float_score_type>( open_penalty   )
		+ 7 * numeric_cast<float_score_type>( extend_penalty ),
		gap_penalty_value_of_alignment(
			make_gap_alignment_of_strings( {
				"x x - x - x - x - x x",
				"x - x - x - x - x - x"
			} ),
			the_gap_penalty
		)
	);
}

/// \brief Check that gap_penalty_value_of_alignment() works for a simple alignment standard holes
BOOST_AUTO_TEST_CASE(standard) {
	const alignment the_alignment = make_gap_alignment_of_strings( {
		"x x x x x x x x x",
		"x x x x x x x - x",
		"x x x - x - x x x",
		"x x x x x x x x x",
		"x - x x x x x x x"
	} );
	BOOST_CHECK_EQUAL(
		4.0 * numeric_cast<float_score_type>( open_penalty ),
		gap_penalty_value_of_alignment( the_alignment, the_gap_penalty )
	);
	BOOST_CHECK_EQUAL(
		4.0,
		gap_count_of_alignment( the_alignment )
	);
	BOOST_CHECK_EQUAL(
		4_z,
		get_naive_num_gaps( the_alignment  )
	);
}

/// \brief Check that gap_penalty_value_of_alignment() works for a hole of width 1 in x out of 4 sequences
BOOST_AUTO_TEST_CASE(hole_in_four) {
	BOOST_CHECK_EQUAL(
		0.0 * numeric_cast<float_score_type>( open_penalty ),
		gap_penalty_value_of_alignment(
			make_gap_alignment_of_strings( {
				"x x x",
				"x x x",
				"x x x",
				"x x x",
			} ),
			the_gap_penalty
		)
	);
	BOOST_CHECK_EQUAL(
		1.0 * numeric_cast<float_score_type>( open_penalty ),
		gap_penalty_value_of_alignment(
			make_gap_alignment_of_strings( {
				"x x x",
				"x x x",
				"x - x",
				"x x x"
			} ),
			the_gap_penalty
		)
	);
	BOOST_CHECK_EQUAL(
		( 4.0 / 3.0 ) * numeric_cast<float_score_type>( open_penalty ),
		gap_penalty_value_of_alignment(
			make_gap_alignment_of_strings( {
				"x x x",
				"x - x",
				"x - x",
				"x x x"
			} ),
			the_gap_penalty
		)
	);
	BOOST_CHECK_EQUAL(
		1.0 * numeric_cast<float_score_type>( open_penalty ),
		gap_penalty_value_of_alignment(
			make_gap_alignment_of_strings( {
				"x - x",
				"x - x",
				"x - x",
				"x x x"
			} ),
			the_gap_penalty
		)
	);
}

/// \brief Check that gap_penalty_value_of_alignment() works for a simple alignment with a mix of issues
BOOST_AUTO_TEST_CASE(mix) {
	BOOST_CHECK_EQUAL(
		3.5 * numeric_cast<float_score_type>( open_penalty ),
		gap_penalty_value_of_alignment(
			make_gap_alignment_of_strings( {
				"x x - - x",
				"x - x - x",
				"x - - x x",
				"x - - x x",
				"x - - x x"
			} ),
			the_gap_penalty
		)
	);
}

/// \brief Check that gap_penalty_value_of_alignment() works for a complex alignment with a mix of issues
BOOST_AUTO_TEST_CASE(complex_mix) {
	BOOST_CHECK_EQUAL(
		  ( 47.0 / 11.0 ) * numeric_cast<float_score_type>( open_penalty   )
		+ ( 50.0 / 11.0 ) * numeric_cast<float_score_type>( extend_penalty ),
		gap_penalty_value_of_alignment(
			make_gap_alignment_of_strings( {
				"x x x x x",
				"x x x x x",
				"x x x x x",

				"x x x - x",
				"x x x - x",
				"x x x - x",
				"x x x - x",

				"x - - - x",
				"x - - - x",
				"x - - - x",
				"x - - - x",
				"x - - - x"
			} ),
			the_gap_penalty
		)
	);
}

BOOST_AUTO_TEST_SUITE_END()
