/// \file
/// \brief The alignment_outputter_fixture class definitions

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

#include "alignment_outputter_fixture.hpp"

#include <boost/test/auto_unit_test.hpp>

#include "acquirer/alignment_acquirer/ssap_aln_file_alignment_acquirer.hpp"
#include "acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp"
#include "alignment/alignment_context.hpp"
#include "chopping/domain/domain.hpp"
#include "test/global_test_constants.hpp"

using namespace cath::align;
using namespace cath::chop;
using namespace cath::opts;
using namespace cath::test;

using std::istringstream;

/// \brief Check that the output stream is empty at the end  of the test
alignment_outputter_fixture::~alignment_outputter_fixture() {
	try {
		BOOST_TEST( out_ss.str().empty() );
	}
	catch (...) {
	}
}

/// \brief Get an example alignment_context
///
/// This only covers the end of the two structures so it can be used for testing alignment_outputter's
/// handling of regions
alignment_context alignment_outputter_fixture::get_example_alignment_context() const {
	constexpr bool       remove_partial_residues = false;
	const     str_vec    ids                     = { "1c0pA01", "1hdoA00" };
	const     domain_vec domains                 = {
		domain{ { make_simple_region( 'A', 1326, 1361 ) }, "1c0pA01_end" },
		domain{ { make_simple_region( 'A',  105,  205 ) }, "1hdoA00_end" }
	};

	istringstream the_istream{};
	return get_alignment_context(
		ssap_aln_file_alignment_acquirer{ global_test_constants::TEST_SOURCE_DATA_DIR()  / "1c0pA01_1hdoA00.refined.sublist" },
		file_list_pdbs_acquirer{ {
			global_test_constants::TEST_SOURCE_DATA_DIR() / "1c0pA01",
			global_test_constants::TEST_SOURCE_DATA_DIR() / "1hdoA00",
		} },
		the_istream,
		remove_partial_residues,
		ids,
		domains
	);
}
