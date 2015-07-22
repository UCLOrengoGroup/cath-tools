/// \file
/// \brief The snap_judgement program_exception_wrapper

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include <boost/units/quantity.hpp>

#include "exception/program_exception_wrapper.h"
#include "scan/scan_tools/all_vs_all.h"
#include "scan/scan_tools/load_and_scan.h"
#include "scan/scan_tools/load_and_scan_metrics.h"
#include "scan/scan_tools/single_pair.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_source_file_set/protein_source_from_pdb.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "test/global_test_constants.h"

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::scan;
using namespace std;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_align_refiner::refine()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class snap_judgement_program_exception_wrapper final : public program_exception_wrapper {
		virtual string do_get_program_name() const override final {
			return "snap-judgement";
		}

		/// \brief Parse the options and then pass them to snap_judgementr::superpose()
		virtual void do_run_program(int /*argc*/, char * /*argv*/[]) override final {
			cerr << "Running snap-judgement\n";

			// The details of the pair to use
			//
			/// \todo Consider whether there are better ways to get the dir/ids here
			///
			/// These are just EXAMPLE_A_PDB_STEMNAME(), EXAMPLE_B_PDB_STEMNAME() and TEST_SOURCE_DATA_DIR()
			/// from global_test_constants 
			const path   &the_dir = "build-test-data";
//			const string &name_a  = "1c0pA01";
//			const string &name_b  = "1hdoA00";
			const string &name_a  = "1n3lA01";
			const string &name_b  = "1r6xA02";

			const auto single_pair_lasm = load_and_scan{
				protein_list_loader{ protein_source_from_pdb(), the_dir, { name_a } },
				protein_list_loader{ protein_source_from_pdb(), the_dir, { name_b } },
				single_pair{}
			}.get_load_and_scan_metrics();

			const auto all_vs_all_ids = str_vec{ "1my7A00", "1my5A00", "2qjyB02", "2qjpB02", "2pw9A02", "2pw9C02", "2c4jA01", "1b4pA01", "2fmpA04", "2vanA03", "1okiA01", "1ytqA01", "1b06A01", "1ma1B01", "1a7sA02", "2xw9A02", "1avyB00", "1avyA00", "1m2tA02", "1hwmA02", "1d0cA01", "1m7vA01", "1a1hA01", "2j7jA03", "1a04A02", "1fseB00", "1fcyA00", "1pzlA00", "1avcA07", "1dk5B01", "1bd8A00", "1s70B01", "1atgA01", "1pc3A01", "1a2oA01", "2ayzA00", "1au7A02", "1rr7A02", "1arbA01", "1si5H01", "1ufmA00", "1a9xB02", "2nv0A00", "1aepA00", "1h6gA02", "1a4iB01", "1sc6A01", "2y1eA01", "1cf7B00", "1a32A00", "1go3F02", "3broD00", "1tnsA00", "2xblD00", "1a3qA01", "1g4mA01", "1a04A01", "2wjwA01", "1a02F00", "1mslA02" };
			const auto all_vs_all_lasm = load_and_scan{
				protein_list_loader{ protein_source_from_pdb(), the_dir, all_vs_all_ids },
				protein_list_loader{ protein_source_from_pdb(), the_dir, all_vs_all_ids },
				all_vs_all{}
			}.get_load_and_scan_metrics();

			cerr <<
R"(SNAP Judgement
===============

Single Pair
-----------

single pair comparison of largish structures 1n3lA01 (209 residues) and 1r6xA02 (213 residues).
SSAP takes about 2.7 seconds and gets a score of 69.28.

)" << to_markdown_string( single_pair_lasm ) << R"(

All versus all within a set of 60 domains
-----------------------------------------

single pair comparison of largish structures 1n3lA01 (209 residues) and 1r6xA02 (213 residues).
SSAPs typically take around a second so 60x60=3600 would normally take around an hour.

)" << to_markdown_string( single_pair_lasm ) << R"(

Build details
-------------

| Platform | Compiler | Library | Boost version |
|----------|----------|---------|---------------|
| )" << BOOST_PLATFORM << " | " << BOOST_COMPILER << " | " << BOOST_STDLIB << " | " << BOOST_LIB_VERSION  << " |" << "\n";
		}
	};
}

/// \brief A main function for snap_judgement that just calls run_program() on a snap_judgement_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::snap_judgement_program_exception_wrapper().run_program( argc, argv );
}
