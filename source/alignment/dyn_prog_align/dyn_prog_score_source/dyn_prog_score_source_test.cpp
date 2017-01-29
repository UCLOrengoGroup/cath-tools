/// \file
/// \brief The dyn_prog_score_source test suite

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

#include "dyn_prog_score_source.hpp"

#include <boost/test/unit_test.hpp>

#include <boost/mpl/list.hpp>

#include "alignment/dyn_prog_align/dyn_prog_score_source/entry_querier_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/test/dyn_prog_score_source_fixture.hpp"
#include "common/size_t_literal.hpp"

 using namespace cath::align;
 using namespace cath::common;

/// \brief TODOCUMENT
///
/// \todo Re-add entry_querier_dyn_prog_score_source to the list (which will require
///       setting up a pair of example proteins for tests).
using all_dyn_prog_score_source_types = boost::mpl::list< /*entry_querier_dyn_prog_score_source,*/
                                                          mask_dyn_prog_score_source,
                                                          old_matrix_dyn_prog_score_source,
                                                          sequence_string_dyn_prog_score_source >;

BOOST_FIXTURE_TEST_SUITE(dyn_prog_score_source_test_suite, dyn_prog_score_source_fixture)

/// \brief Check that each dyn_prog_score_source returns two strictly positive lengths
///        and gets a score for indices 0, 0 without throwing
BOOST_AUTO_TEST_CASE_TEMPLATE(basic_indices_and_score, dyn_prog_score_source_type, all_dyn_prog_score_source_types) {
	const dyn_prog_score_source_type score_source = make_example_dyn_prog_score_source<dyn_prog_score_source_type>();
	BOOST_REQUIRE_GT(     score_source.get_length_a(), 0_z );
	BOOST_REQUIRE_GT(     score_source.get_length_b(), 0_z );
	BOOST_CHECK_NO_THROW( score_source.get_score(0, 0) );
}

BOOST_AUTO_TEST_SUITE_END()
