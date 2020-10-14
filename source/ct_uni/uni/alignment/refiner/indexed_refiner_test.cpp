/// \file
/// \brief The indexed_refiner test suite

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

#include <boost/range/adaptor/transformed.hpp>

#include "alignment/io/alignment_io.hpp"
#include "alignment/refiner/indexed_refiner.hpp"
#include "structure/protein/protein_source_file_set/protein_from_pdb.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::align;

using std::string;

namespace cath {
	namespace test {

		/// \brief The indexed_refiner_test_suite_fixture to assist in testing indexed_refiner
		struct indexed_refiner_test_suite_fixture {
		protected:
			~indexed_refiner_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(indexed_refiner_test_suite, cath::test::indexed_refiner_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	// const protein prot_1qbkB00 = read_protein_from_files( protein_from_pdb(), "refine_stuff/really_long", "1qbkB00" );
	// const protein prot_1u6gC00 = read_protein_from_files( protein_from_pdb(), "refine_stuff/really_long", "1u6gC00" );
	// const alignment the_aln = read_alignment_from_cath_ssap_legacy_format( "refine_stuff/really_long/1qbkB001u6gC00.list", prot_1qbkB00, prot_1u6gC00 );

	// do_some_gubbins( prot_1qbkB00, prot_1u6gC00, the_aln );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
