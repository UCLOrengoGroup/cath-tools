/// \file


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

#include "common/boost_check_no_throw_diag.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "file/pdb/pdb_atom.hpp"

#include <iostream>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The pdb_atom_test_suite_fixture to assist in testing pdb_atom
		struct pdb_atom_test_suite_fixture {
		protected:
			~pdb_atom_test_suite_fixture() noexcept = default;

		public:
			const string ATOM_RECORD_SIMPLE           = { "ATOM      1  N   LEU A 999       0.041 148.800  54.967  1.00 35.61"             };
			const string ATOM_RECORD_WITH_SPACE_IN_AA = { "ATOM  31067  O5'  DT T  10      91.625  85.910 -23.177  0.00119.41"             };
			const string ATOM_RECORD_NEG_RES_NUM      = { "ATOM   2041  N   ASN B  -1     -27.445  -1.104  16.047  1.00 65.16"             };
			const string ATOM_RECORD_RALIGNED_RES_NUM = { "ATOM     89  N   GLU   200     -60.412  54.016  32.157  1.00  0.00           N" };
			const string ATOM_RECORD_LALIGNED_RES_NUM = { "ATOM    182  N   GLU  200       -0.292   3.837   4.911  1.00226.06           N" };

			/// \brief Check that parsing the string into a pdb_atom and then
			///        writing it back to a string produces identical results
			void check_parse_and_write_pdb_atom(const string &arg_input_string
		// 	                                    const size_t &arg_atom_index
			                                    ) {
				const chain_resname_atom_tuple parsed_details( parse_pdb_atom_record( arg_input_string ) );
				ostringstream pdb_atom_ss;
				write_pdb_file_entry(
					pdb_atom_ss,
					get<0>( parsed_details ),
					get<1>( parsed_details ),
					get<2>( parsed_details )
				);
				if (arg_input_string != pdb_atom_ss.str()) {
					cerr << endl;
					cerr << "Expected : \"" << arg_input_string  << "\"" << endl;
					cerr << "Got      : \"" << pdb_atom_ss.str() << "\"" << endl;
				}
				BOOST_CHECK_EQUAL(arg_input_string, pdb_atom_ss.str());
			}
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(pdb_atom_test_suite, cath::test::pdb_atom_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(parse_from_simple_line) {
	const chain_resname_atom_tuple parsed_details( parse_pdb_atom_record(
		ATOM_RECORD_SIMPLE
	) );
	const pdb_atom &my_atom( get<2>( parsed_details ) );

	BOOST_CHECK_EQUAL(   0.041, my_atom.get_coord().get_x() );
	BOOST_CHECK_EQUAL( 148.800, my_atom.get_coord().get_y() );
	BOOST_CHECK_EQUAL(  54.967, my_atom.get_coord().get_z() );
}

/// \brief Check that reading/writing atom records handles a standard case
BOOST_AUTO_TEST_CASE(parse_and_write_simple) {
	check_parse_and_write_pdb_atom(ATOM_RECORD_SIMPLE);
}

/// \brief Check that reading/writing atom records handles a space in the amino acid name
BOOST_AUTO_TEST_CASE(parse_and_write_with_space_in_aa) {
	BOOST_CHECK_THROW( parse_pdb_atom_record( ATOM_RECORD_WITH_SPACE_IN_AA ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_pdb_atom_record( ATOM_RECORD_WITH_SPACE_IN_AA ), invalid_argument_exception );
}

/// \brief Check that reading/writing atom records handles a negative residue number
BOOST_AUTO_TEST_CASE(parse_and_write_neg_res_num) {
	check_parse_and_write_pdb_atom(ATOM_RECORD_NEG_RES_NUM);
}

/// \brief Check that ATOM parsing doesn't permit a right aligned residue number if that's specified
BOOST_AUTO_TEST_CASE(requiring_right_aligned_res_num_throws_if_not_right_aligned) {
	BOOST_CHECK_NO_THROW_DIAG( parse_pdb_atom_record( ATOM_RECORD_RALIGNED_RES_NUM )                             );
	BOOST_CHECK_THROW(         parse_pdb_atom_record( ATOM_RECORD_LALIGNED_RES_NUM ), invalid_argument_exception );
}

BOOST_AUTO_TEST_SUITE_END()