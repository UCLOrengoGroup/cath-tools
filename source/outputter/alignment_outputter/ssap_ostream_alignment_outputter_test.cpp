/// \file
/// \brief The ssap_ostream_alignment_outputter test suite

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

#include "alignment/alignment_context.hpp"
#include "common/test_predicate/files_equal.hpp"
#include "outputter/alignment_outputter/alignment_outputter_fixture.hpp"
#include "outputter/alignment_outputter/file_alignment_outputter.hpp"
#include "outputter/alignment_outputter/ssap_ostream_alignment_outputter.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::opts;
using namespace cath::test;

BOOST_FIXTURE_TEST_SUITE(ssap_ostream_alignment_outputter_test_suite, alignment_outputter_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	file_alignment_outputter{ out_file, ssap_ostream_alignment_outputter() }
		.output_alignment( get_example_alignment_context(), out_ss );

	BOOST_CHECK_FILES_EQUAL( out_file, global_test_constants::TEST_SOURCE_DATA_DIR() / "1c0pA01_1hdoA00.partial_aln.ssap" );
}

BOOST_AUTO_TEST_SUITE_END()
