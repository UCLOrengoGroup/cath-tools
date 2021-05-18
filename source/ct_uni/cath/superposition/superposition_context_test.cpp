/// \file
/// \brief The superposition_context test suite

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

#include <filesystem>

#include <boost/test/unit_test.hpp>

#include "cath/biocore/residue_id.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/boost_addenda/range/back.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/property_tree/from_json_string.hpp"
#include "cath/common/property_tree/to_json_string.hpp"
#include "cath/file/options/data_dirs_options_block.hpp"
#include "cath/superposition/superposition_content_spec.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"
#include "cath/test/superposition_fixture.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;
using namespace ::cath::sup;

using ::std::filesystem::path;

namespace cath {
	namespace test {

		/// \brief The superposition_context_test_suite_fixture to assist in testing superposition_context
		struct superposition_context_test_suite_fixture : protected superposition_fixture {
		protected:
			~superposition_context_test_suite_fixture() noexcept = default;
		};
	}
}

BOOST_FIXTURE_TEST_SUITE(superposition_context_test_suite, cath::test::superposition_context_test_suite_fixture)

BOOST_AUTO_TEST_CASE(load_pdbs_from_names_copy_leaves_orig_empty_pdbs) {
	const auto loaded_sup_con = load_pdbs_from_names_copy( the_sup_con, build_data_dirs_spec_of_dir( TEST_SOURCE_DATA_DIR() ) );
	BOOST_REQUIRE_EQUAL( get_pdbs( the_sup_con ).size(), 2 );
	BOOST_CHECK_EQUAL  ( get_pdbs( the_sup_con )[ 0 ].get_num_residues(), 0 );
	BOOST_CHECK_EQUAL  ( get_pdbs( the_sup_con )[ 1 ].get_num_residues(), 0 );
}

BOOST_AUTO_TEST_CASE(load_pdbs_from_names_copy_sets_two_pdbs_with_199_and_205_residues) {
	const auto loaded_sup_con = load_pdbs_from_names_copy( the_sup_con, build_data_dirs_spec_of_dir( TEST_SOURCE_DATA_DIR() ) );
	BOOST_REQUIRE_EQUAL( get_pdbs( loaded_sup_con ).size(), 2 );
	BOOST_CHECK_EQUAL  ( get_pdbs( loaded_sup_con )[ 0 ].get_num_residues(), 199 );
	BOOST_CHECK_EQUAL  ( get_pdbs( loaded_sup_con )[ 1 ].get_num_residues(), 205 );
}

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup_con) {
	BOOST_CHECK_EQUAL(
		to_json_string( the_sup_con, json_style::COMPACT ),
		sup_context_json_str
	);
}

BOOST_AUTO_TEST_CASE(from_json_string_works) {
	const auto from_json_str = from_json_string<superposition_context>( sup_context_json_str );
	BOOST_REQUIRE_EQUAL     ( get_pdbs ( from_json_str ).size(),                  2       );
	BOOST_CHECK_EQUAL       ( get_pdbs ( from_json_str )[ 0 ].get_num_residues(), 0       );
	BOOST_CHECK_EQUAL       ( get_pdbs ( from_json_str )[ 1 ].get_num_residues(), 0       );

	BOOST_CHECK_EQUAL_RANGES( get_name_sets( from_json_str ),                     names   );

	BOOST_CHECK_EQUAL       ( from_json_str.get_superposition(),                  the_sup );

	BOOST_CHECK_EQUAL       ( from_json_str.has_alignment(),                      false   );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_CASE(get_supn_content_pdb_works) {
	const path pdb_file = TEST_SOURCE_DATA_DIR() / "supn_content" / "1bdh";
	const auto the_pdb = read_pdb_file( pdb_file );

	const region_vec regions_1bdhA01 = { make_simple_region( 'A',   3,  59 )                                      };
	const region_vec regions_1bdhA02 = { make_simple_region( 'A',  60, 161 ), make_simple_region( 'A', 292, 323 ) };
	const region_vec regions_1bdhA03 = { make_simple_region( 'A', 162, 291 ), make_simple_region( 'A', 324, 340 ) };

	const pdb supn_content_pdb_1bdhA01 = get_supn_content_pdb( the_pdb, regions_1bdhA01, superposition_content_spec{} );
	const pdb supn_content_pdb_1bdhA02 = get_supn_content_pdb( the_pdb, regions_1bdhA02, superposition_content_spec{} );
	const pdb supn_content_pdb_1bdhA03 = get_supn_content_pdb( the_pdb, regions_1bdhA03, superposition_content_spec{} );

	BOOST_CHECK_EQUAL(        supn_content_pdb_1bdhA01.get_num_residues(),                          74                         );
	BOOST_CHECK_EQUAL(        supn_content_pdb_1bdhA01.get_post_ter_residues().size(),               8                         );
	BOOST_CHECK_EQUAL( front( supn_content_pdb_1bdhA01                         ).get_residue_id(), make_residue_id( 'B', 699 ) );
	BOOST_CHECK_EQUAL( back ( supn_content_pdb_1bdhA01                         ).get_residue_id(), make_residue_id( 'A',  59 ) );
	BOOST_CHECK_EQUAL( front( supn_content_pdb_1bdhA01.get_post_ter_residues() ).get_residue_id(), make_residue_id( 'B', 765 ) );
	BOOST_CHECK_EQUAL( back ( supn_content_pdb_1bdhA01.get_post_ter_residues() ).get_residue_id(), make_residue_id( 'A', 777 ) );

	BOOST_CHECK_EQUAL(        supn_content_pdb_1bdhA02.get_num_residues(),                         134                         );
	BOOST_CHECK_EQUAL(        supn_content_pdb_1bdhA02.get_post_ter_residues().size(),              32                         );
	BOOST_CHECK_EQUAL( front( supn_content_pdb_1bdhA02                         ).get_residue_id(), make_residue_id( 'A',  60 ) );
	BOOST_CHECK_EQUAL( back ( supn_content_pdb_1bdhA02                         ).get_residue_id(), make_residue_id( 'A', 323 ) );
	BOOST_CHECK_EQUAL( front( supn_content_pdb_1bdhA02.get_post_ter_residues() ).get_residue_id(), make_residue_id( 'A', 599 ) );
	BOOST_CHECK_EQUAL( back ( supn_content_pdb_1bdhA02.get_post_ter_residues() ).get_residue_id(), make_residue_id( 'A', 788 ) );

	BOOST_CHECK_EQUAL(        supn_content_pdb_1bdhA03.get_num_residues(),                         147                         );
	BOOST_CHECK_EQUAL(        supn_content_pdb_1bdhA03.get_post_ter_residues().size(),              36                         );
	BOOST_CHECK_EQUAL( front( supn_content_pdb_1bdhA03                         ).get_residue_id(), make_residue_id( 'A', 162 ) );
	BOOST_CHECK_EQUAL( back ( supn_content_pdb_1bdhA03                         ).get_residue_id(), make_residue_id( 'A', 340 ) );
	BOOST_CHECK_EQUAL( front( supn_content_pdb_1bdhA03.get_post_ter_residues() ).get_residue_id(), make_residue_id( 'A', 599 ) );
	BOOST_CHECK_EQUAL( back ( supn_content_pdb_1bdhA03.get_post_ter_residues() ).get_residue_id(), make_residue_id( 'A', 788 ) );
}



//  To calculate the figures on the furthest distances from CA atoms, scanned a large number of PDB files and recorded the
//  largest distance. Was going to take 4th largest (or 1st for apart from ASX and GLX), multiply by 1.1 and round to 3sf
//  and result of that is second column.
//  ALA      3.62     ( 6.63459 6.08628 4.25559  3.2929 3.23273 3.22772 3.20309 3.07452 2.92904 2.90389 )
//  ARG      9.22     ( 11.3549 8.38962 8.38332 8.38293  8.3727 8.36276 8.34231 8.33822 8.33058 8.32651 )
//  ASN      5.21     ( 6.81096 5.06171 4.73937 4.73857 4.72852 4.72647  4.7062 4.69577 4.68606 4.67156 )
//  ASP      4.97     ( 8.54586 5.07719 4.89187 4.51591 4.49497 4.49368 4.49313  4.4878 4.48468 4.47677 )
//  ASX      3.75     ( 3.41174 3.33253                                                                 )
//  CYS      4.32     ( 4.12741 3.96795 3.93771 3.92844 3.92687 3.91458 3.91064 3.90452 3.89417 3.88587 )
//  GLN      6.65     ( 6.08122 6.05941 6.05692  6.0458 6.02797 5.99374  5.9786  5.9725 5.96889  5.9618 )
//  GLU      6.78     (  36.059 6.16626 6.16616 6.16609 6.16595 6.16593 6.16582 6.16578 6.16554 6.15341 )
//  GLX      4.32     (  3.9289                                                                         )
//  GLY      5.21     ( 9.38445 6.74371 4.97226 4.73755 3.65506 3.50463 3.34141 3.29895 3.18099 3.00473 )
//  HIS      6.51     ( 5.93011  5.9212 5.91549 5.91509 5.90462 5.90172 5.89601  5.8907  5.8819 5.87362 )
//  ILE      5.36     ( 4.98606 4.91932  4.9119 4.87389 4.86807 4.86785 4.85774 4.85208  4.8416 4.84157 )
//  LEU      5.49     ( 18.2554 5.15672 5.00275 4.99124 4.98429 4.97682 4.96843 4.96836 4.96429 4.96271 )
//  LYS      8.90     ( 85.6703 84.4392 10.9482  8.0909 7.77396 7.61505 7.54366  7.5356 7.48096 7.37498 )
//  MET      7.03     ( 6.51569 6.46256 6.40321 6.38647 6.37399 6.36878 6.36691 6.36144 6.35938 6.35769 )
//  PHE      7.19     ( 6.59313 6.56639 6.54296 6.53252  6.4929 6.48815 6.47493  6.4728 6.46114 6.45432 )
//  PRO      3.83     ( 3.74101 3.54677 3.53948 3.48444  3.4463 3.44367 3.43948  3.4392 3.43858 3.43731 )
//  SEC      3.16     ( 2.92998 2.90802  2.8932 2.87171 2.85696 2.80445                                 )
//  SER      3.76     ( 6.40928 3.47836 3.43206 3.42134 3.41762 3.41595 3.41149 3.40422 3.40274 3.40212 )
//  THR      4.10     ( 4.07216 3.75104 3.73069 3.72622 3.72524 3.72437 3.67135  3.6619 3.63288  3.6324 )
//  TRP      8.73     ( 7.95336 7.94141 7.93346 7.93188 7.92058 7.91838 7.91522 7.91332 7.91034 7.90809 )
//  TYR      8.17     ( 7.60994 7.45641 7.42973  7.4281 7.36718 7.36226 7.31645  7.3126 7.30323 7.30236 )
//  UNK      5.22     ( 5.15753 5.15141 4.95181 4.74281 4.68458 4.68211 4.39289 4.36692 4.32534 4.22927 )
//  VAL      4.14     ( 4.45404 3.93986  3.8088 3.76094 3.66825 3.65856 3.65216 3.64956 3.64491 3.64453 )
// All others would have been taken to be 14.000
//
//
//
//
// #include <boost/range/algorithm/sort.hpp> // ***** TEMPORARY *****
// #include <boost/range/algorithm/lower_bound.hpp> // ***** TEMPORARY *****
// #include "cath/common/boost_addenda/range/max_proj_element.hpp" // ***** TEMPORARY *****
//
// using ::boost::range::sort;
// using ::cath::file::pdb_atom;
// using ::cath::file::pdb_residue;
// using ::cath::file::read_pdb_file;
// using ::cath::path_vec;
// using ::cath::str_doub_map;
//
// BOOST_AUTO_TEST_CASE(scan_all_pdbs_for_max_intra_residue_to_ca_dists) {
// 	constexpr size_t NUM_MAX_DISTS = 10;
// 	const path pdb_dir = "/cath/data/v4_1_0/wholepdb";
//
// 	BOOST_REQUIRE( is_directory( pdb_dir  ) );
//
// 	const path_vec sorted_pdb_files = [&] {
// 		path_vec results;
// 		for (const directory_entry &pdb_entry : directory_iterator( pdb_dir ) ) {
// 			const path &pdb_file = pdb_entry.path();
// 			if ( is_regular_file( pdb_file ) || is_symlink( pdb_file ) ) {
// 				results.push_back( pdb_file );
// 			}
// 		}
// 		boost::range::sort( results );
// 		return results;
// 	} ();
//
// 	using str_doub_vec_map = map<string, cath::doub_vec>;
// 	str_doub_vec_map max_ca_to_atom_dists_by_aa;
//
// 	for (const path &pdb_file : sorted_pdb_files ) {
// 		const auto parsed_pdb = read_pdb_file( pdb_file );
//
// 		for (const pdb_residue &res : parsed_pdb) {
// 			const double max_dist = max_proj(
// 				res,
// 				std::less<>{},
// 				[&] (const pdb_atom &x) {
// 					return ( res.has_carbon_alpha() && res.get_carbon_alpha().get_alt_locn() == x.get_alt_locn() )
// 						? distance_between_points( res.get_carbon_alpha().get_coord(), x.get_coord() )
// 						: 0.0;
// 				}
// 			);
// 			// std::cerr << res.get_amino_acid() << "\t" << max_dist << "\n";
// 			const auto aa_code = res.get_amino_acid().get_code();
// 			// \todo Come C++17, do this with insert_or_assign()
// 			if ( ! contains( max_ca_to_atom_dists_by_aa, aa_code ) ) {
// 				max_ca_to_atom_dists_by_aa.emplace( aa_code, cath::doub_vec{ max_dist } );
// 			}
// 			else {
// 				auto &dists = max_ca_to_atom_dists_by_aa.at( aa_code );
// 				if ( dists.size() < NUM_MAX_DISTS || dists.back() < max_dist ) {
// 					const auto posn_itr = boost::range::lower_bound( dists, max_dist, greater<>{} );
// 					if ( posn_itr == cend( dists ) || *posn_itr != max_dist ) {
// 						dists.insert(
// 							posn_itr,
// 							max_dist
// 						);
// 						if ( dists.size() > NUM_MAX_DISTS ) {
// 							dists.pop_back();
// 						}
// 					}
// 				}
// 			}
// 		}

// 		// break;
// 		for (const auto &x : max_ca_to_atom_dists_by_aa) {
// 			std::cerr << x.first;
// 			for (const double &y : x.second) {
// 				std::cerr << std::right << std::setw( 8 ) << y;
// 			}
// 			std::cerr << "\n";
// 		}
// 		std::cerr << "\n";
// 	}

// }

BOOST_AUTO_TEST_SUITE_END()

