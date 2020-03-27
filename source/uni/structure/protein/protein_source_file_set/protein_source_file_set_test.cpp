/// \file
/// \brief The protein_source_file_set test suite

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
#include <boost/test/unit_test.hpp>

#include "common/file/simple_file_read_write.hpp"
#include "file/options/data_dirs_spec.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_source_file_set/protein_from_pdb_dssp_and_sec.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

#include <iostream>
#include <random>

using namespace cath::common;
using namespace cath::opts;
using namespace cath;

using std::string;

BOOST_AUTO_TEST_SUITE(protein_source_file_set_test_suite)

BOOST_AUTO_TEST_CASE(file_set_contains_primary_file) {
	const auto all_source_file_sets = get_all_protein_source_file_sets();
	for (const auto &x : all_source_file_sets) {
		BOOST_TEST( contains( x.get_file_set(), x.get_primary_file() ) );
	}
}

// Commented-out test to load a large list of proteins via their PDB+DSSP+SEC files to check for any breaking errors
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );

//	const path    list_file     = "/cath/homes2/ucbctnl/sec_dssp_pdb_ids.txt";
//	const path    root_dir      = "/cath/mothra-data1/people/ucbctnl/temp_tony_20160601";
//	const str_vec ids           = read_file<string>( list_file );
//	const auto    the_data_dirs = data_dirs_spec{}.set_cath_root_dir( root_dir );
//	for (const string id : ids) {
////		std::mt19937 rng{ std::random_device{}() };
////		std::shuffle( ::std::begin( ids ), ::std::end( ids ), rng );
//		cerr << "Attempting to load PDB, DSSP and SEC for " << id << "\n";
//		protein_from_pdb_dssp_and_sec{}.read_files(
//			the_data_dirs,
//			id,
//			cerr
//		);
//	}
}

BOOST_AUTO_TEST_SUITE_END()
