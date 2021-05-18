/// \file
/// \brief The view_cache_index test suite

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

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/io/alignment_io.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/ssap/ssap.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"
#include "cath/structure/view_cache/index/detail/vcie_match_criteria.hpp"
#include "cath/structure/view_cache/index/quad_find_action.hpp"
#include "cath/structure/view_cache/index/quad_find_action_check.hpp"
#include "cath/structure/view_cache/index/view_cache_index.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::index;
using namespace ::cath::index::detail;
using namespace ::cath::geom;
using namespace ::std;

namespace cath {
	namespace test {

		/// \brief The view_cache_index_test_suite_fixture to assist in testing view_cache_index
		struct view_cache_index_test_suite_fixture : protected global_test_constants {
		protected:
			~view_cache_index_test_suite_fixture() noexcept = default;
		};

	} // namespace test
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(view_cache_index_test_suite, cath::test::view_cache_index_test_suite_fixture)

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(basic) {
//	// const double cell_width = 0.5 * sqrt( 40.0 );
//	const double cell_width =       sqrt( 40.0 );
////	const double cell_width = 2.0 * sqrt( 40.0 );
//
////	const protein protein_a = read_protein_from_files( protein_from_pdb(), TEST_SOURCE_DATA_DIR(), "1a04A02" );
////	const protein protein_b = read_protein_from_files( protein_from_pdb(), TEST_SOURCE_DATA_DIR(), "1fseB00" );
//////    const alignment correct_alignment = read_alignment_from_fasta_file(
//////    	TEST_SOURCE_DATA_DIR() / "1a04A02_1fseB00.fa",
//////		make_pdb_list( { read_pdb_file( TEST_SOURCE_DATA_DIR() / "1a04A02" ),
//////		                 read_pdb_file( TEST_SOURCE_DATA_DIR() / "1fseB00" ) } ),
//////		cerr
//////	);
//
////	//	1n3lA01  1h3fA01  209  195  84.91  186   88   26   2.68
////	const protein protein_a = read_protein_from_files( protein_from_pdb(), TEST_SOURCE_DATA_DIR(), "1n3lA01" );
////	const protein protein_b = read_protein_from_files( protein_from_pdb(), TEST_SOURCE_DATA_DIR(), "1h3fA01" );
//////    const alignment correct_alignment = read_alignment_from_fasta_file(
//////    	TEST_SOURCE_DATA_DIR() / "1n3lA01_1h3fA01.fa",
//////		make_pdb_list( { read_pdb_file( TEST_SOURCE_DATA_DIR() / "1n3lA01" ),
//////		                 read_pdb_file( TEST_SOURCE_DATA_DIR() / "1h3fA01" ) } ),
//////		cerr
//////	);
//
//	//	1n3lA01  1r6xA02  209  213  69.28  129   60    5   6.02
//	const protein protein_a = read_protein_from_files( protein_from_pdb(), TEST_SOURCE_DATA_DIR(), "1n3lA01" );
//	const protein protein_b = read_protein_from_files( protein_from_pdb(), TEST_SOURCE_DATA_DIR(), "1r6xA02" );
////    const alignment correct_alignment = read_alignment_from_fasta_file(
////    	TEST_SOURCE_DATA_DIR() / "1n3lA01_1r6xA02.fa",
////		make_pdb_list( { read_pdb_file( TEST_SOURCE_DATA_DIR() / "1n3lA01" ),
////		                 read_pdb_file( TEST_SOURCE_DATA_DIR() / "1r6xA02" ) } ),
////		cerr
////	);
//
////    34 and 35 vs 18 to 19
//
////	view_cache_index new_view_cache_index( cell_width );
////	const view_cache_index view_cache_index_a = build_view_cache_index( cell_width, protein_a );
////	const view_cache_index view_cache_index_b = build_view_cache_index( cell_width, protein_b );
//
//	const auto the_criteria = make_default_vcie_match_criteria();
//
////	const size_t num_repeats =  500;
////	const size_t num_repeats = 20;
//	const size_t num_repeats = 1;
//	using durn_vec = vector<high_resolution_clock::duration>;
//	durn_vec durations;
//	for (const size_t &repeat_ctr : indices( num_repeats ) ) {
//		quad_find_action_check indexed_action( protein_a, protein_b, the_criteria );
//		durations.push_back( process_quads_indexed(
//			protein_a,
//			protein_b,
//			cell_width,
//			// make_angle_from_degrees<angle_base_type>( 67.5 / 2.0 ),
//			// make_angle_from_degrees<angle_base_type>( 67.5 / 2.0 ),
//			make_angle_from_degrees<angle_base_type>( 67.5 ),
//			make_angle_from_degrees<angle_base_type>( 67.5 ),
//			the_criteria,
//			indexed_action
//		) );
//		cerr << "The indexed score is    : " << indexed_action.get_total_score() << endl;
//		cerr << "The indexed scan takes  : " << numeric_cast<double>( duration_cast<nanoseconds>( durations.back() ).count() ) / 1000000.0 << "ms" << endl;
//	}
//
//	quad_find_action complete_action( protein_a, protein_b );
//	const auto complete_duration = process_quads_complete(
//		protein_a,
//		protein_b,
//		cell_width,
//		the_criteria,
//		complete_action
//	);
//	cerr << "The complete score is   : " << complete_action.get_total_score() << endl;
//	cerr << "The complete scan takes : " << numeric_cast<double>( duration_cast<nanoseconds>( complete_duration ).count() ) / 1000000.0 << "ms" << endl;
//
////    const duration total_duration = accumulate(
////    	durations,
////		duration( 0, 0, 0, 0 )
////    );
////    const duration average_duration   = total_duration / num_repeats;
////    const double        avg_num_per_second = ( 1000000.0 / numeric_cast<double>( duration_cast<microseconds>( average_duration.total_microseconds() ).count() ) );
////    cerr << "Durations : " << endl;
////    for (const duration &the_durn : durations) {
////    	cerr << "  duration : " << the_durn.total_microseconds() << " microseconds" << endl;
////    }
////    cerr << "Average duration was " << average_duration.total_microseconds() << " microseconds";
////    cerr << ", (which equates to "  << avg_num_per_second                    << " per second)";
////    cerr << endl;
//
////    protein_a
//
////    187 H 0 V   91  V 0 H   43
//
//    BOOST_CHECK( true );
//}

BOOST_AUTO_TEST_SUITE_END()

