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
#include <boost/test/auto_unit_test.hpp>

#include "common/file/simple_file_read_write.h"
#include "options/options_block/data_dirs_spec.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_source_file_set/protein_source_from_pdb_dssp_and_sec.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

#include <iostream>
#include <random>

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::opts;
using namespace cath;

using std::cerr;
using std::string;

BOOST_AUTO_TEST_SUITE(protein_source_file_set_test_suite)

// Commented-out test to load a large list of proteins via their PDB+DSSP+SEC files to check for any breaking errors
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );

//	const path    list_file     = "/cath/homes2/ucbctnl/sec_dssp_pdb_ids.txt";
//	const path    root_dir      = "/cath/mothra-data1/people/ucbctnl/temp_tony_20160601";
//	const str_vec ids           = read_file<string>( list_file );
//	const auto    the_data_dirs = data_dirs_spec{}.set_cath_root_dir( root_dir );
//	for (const string id : ids) {
//		if ( id == "2br0A" || id == "2c28A" || id == "3nirA" || id == "4c9aB" || id == "4c9aD" ) {
//			continue;
//		}
////		std::mt19937 rng{ std::random_device{}() };
////		std::shuffle( ::std::begin( ids ), ::std::end( ids ), rng );
//		cerr << "Attempting to load PDB, DSSP and SEC for " << id << "\n";
//		protein_source_from_pdb_dssp_and_sec{}.read_files(
//			the_data_dirs,
//			id,
//			cerr
//		);
//	}
}

BOOST_AUTO_TEST_SUITE_END()
