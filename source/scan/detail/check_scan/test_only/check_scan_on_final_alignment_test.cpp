/// \file
/// \brief The check_scan_on_final_alignment test suite

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

#include <boost/filesystem/path.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "alignment/alignment.hpp"
#include "alignment/residue_score/residue_scorer.hpp"
#include "scan/detail/check_scan/test_only/alignment_scan_comparison.hpp"
#include "scan/detail/check_scan/test_only/check_scan_on_final_alignment.hpp"
#include "scan/detail/scan_role.hpp"
#include "scan/detail/stride/roled_scan_stride.hpp"
#include "scan/quad_criteria.hpp"
#include "scan/scan_stride.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/protein_source_file_set/protein_source_file_set.hpp"
#include "structure/protein/protein_source_file_set/protein_source_from_pdb.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/global_test_constants.hpp"

#include <iostream>

using namespace cath;
using namespace cath::align;
using namespace cath::scan;
using namespace cath::scan::detail;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The check_scan_on_final_alignment_test_suite_fixture to assist in testing check_scan_on_final_alignment
		struct check_scan_on_final_alignment_test_suite_fixture : protected global_test_constants {
		protected:
			~check_scan_on_final_alignment_test_suite_fixture() noexcept = default;
		};

	}
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(check_scan_on_final_alignment_test_suite, cath::test::check_scan_on_final_alignment_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
//	//SNAP: -100.00; SSAP: 1a1hA01 2j7jA03 30 27 88.49 26 86 18  1.59 very similar; short
//	//SNAP: -100.00; SSAP: 2j7jA03 1a1hA01 27 30 88.49 26 86 18  1.59 very similar; short
//	//SNAP:   23.47; SSAP: 1a04A02 1au7A02 80 57 77.50 52 65 10  7.62
//	//SNAP: -100.00; SSAP: 1fseB00 1au7A02 70 57 80.82 52 74 12  7.46 not too bad
//	//SNAP:    9.79; SSAP: 1fseB00 1ufmA00 70 84 76.72 66 78 11  9.68 pretty good for much but different at end
//	//SNAP:   15.76; SSAP: 1au7A02 1fseB00 57 70 80.82 52 74 12  7.46
//	//SNAP: -100.00; SSAP: 1au7A02 1rr7A02 57 48 81.49 44 77  6  4.19
//	//SNAP: -100.00; SSAP: 1au7A02 1ufmA00 57 84 78.07 54 64  3  3.93
//	//SNAP: -100.00; SSAP: 1rr7A02 1fseB00 48 70 81.24 48 68 10  4.24 not bad
//	//SNAP: -100.00; SSAP: 1ufmA00 1fseB00 84 70 76.68 68 80 11 10.06
//	//SNAP: -100.00; SSAP: 1ufmA00 1cf7B00 84 82 79.64 63 75  2  3.77 pretty good for much but bit different in middle
//	ostringstream parse_ss;
//	const auto the_scofa       = check_scan_on_final_alignment{};
////	const auto the_scan_stride = scan_stride{ 4, 4, 2, 2 };
//	const auto the_scan_stride = scan_stride{ 2, 2, 1, 1 };
////	const auto the_scan_stride = scan_stride{ 1, 1, 0, 0 };
////	const auto the_scan_stride = scan_stride{ 0, 0, 0, 0 };
//	//                                             0          1          2          3          4          5          6          7
//	const auto ids                  = str_vec{ "1a04A02", "1a1hA01", "1au7A02", "1cf7B00", "1fseB00", "1rr7A02", "1ufmA00", "2j7jA03" };
//	const auto proteins             = read_proteins_from_files( protein_source_from_pdb(), TEST_EXAMPLE_PDBS_DATA_DIR(), ids );
//	const auto prot_1ufmA00         = proteins[ 6 ];
//	const auto prot_1cf7B00         = proteins[ 3 ];
//	const auto prot_1ufmA00_1cf7B00 = make_protein_list( { prot_1ufmA00, prot_1cf7B00 } );
//	const auto aln_1c55A_1c55A      = read_and_rescore_fasta_alignment( TEST_EXAMPLE_PDBS_DATA_DIR() / "1ufmA00_1cf7B00.fa", prot_1ufmA00_1cf7B00, residue_scorer(), parse_ss );
//
//	const auto the_comparison = the_scofa.do_check(
//		aln_1c55A_1c55A,
//		prot_1ufmA00,
//		prot_1cf7B00,
//		make_default_quad_criteria(),
//		the_scan_stride
//	);
//	cerr << "the_comparison is : " << the_comparison << "\n";
//
//	print_highlight_rep_pymol_commands( cerr, the_scofa, prot_1ufmA00, roled_scan_stride{ scan_role::QUERY, the_scan_stride } );
//	print_highlight_rep_pymol_commands( cerr, the_scofa, prot_1cf7B00, roled_scan_stride{ scan_role::INDEX, the_scan_stride } );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
