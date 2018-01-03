/// \file
/// \brief The pdb test suite

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

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "chopping/region/region.hpp"
#include "common/boost_addenda/log/log_to_ostream_guard.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/file/temp_file.hpp"
#include "common/size_t_literal.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "test/boost_addenda/boost_check_equal_ranges.hpp"
#include "test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "test/global_test_constants.hpp"
#include "test/predicate/files_equal.hpp"

#include <regex>
#include <vector>

using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;

using boost::algorithm::icontains;
using boost::algorithm::join;
using boost::filesystem::path;
using boost::irange;
using std::istringstream;
using std::ostringstream;
using std::regex;
using std::string;
using std::stringstream;

namespace cath {
	namespace test {

		/// \brief The pdb_test_suite_fixture to assist in testing pdb
		struct pdb_test_suite_fixture : protected global_test_constants {
		protected:
			~pdb_test_suite_fixture() noexcept = default;

			/// \brief Check that the number of atoms in a vector of PDBs matches the expected numbers
			void check_nums_of_atoms(const pdb_list &,
			                         const size_vec &) const;
		};

	}  // namespace test
}  // namespace cath

/// \brief Check that the number of atoms in a vector of PDBs matches the expected numbers
void cath::test::pdb_test_suite_fixture::check_nums_of_atoms(const pdb_list &arg_pdbs,                  ///< TODOCUMENT
                                                             const size_vec &arg_expected_nums_of_atoms ///< TODOCUMENT
                                                             ) const {
	// Grab the numbers of atoms
	size_vec got_nums_of_atoms;
	got_nums_of_atoms.reserve(arg_pdbs.size());
	for (const pdb &arg_pdb : arg_pdbs) {
		got_nums_of_atoms.push_back( arg_pdb.get_num_atoms() );
	}

	// Check that they match
	BOOST_CHECK_EQUAL_RANGES( arg_expected_nums_of_atoms, got_nums_of_atoms );
}

BOOST_FIXTURE_TEST_SUITE(pdb_test_suite, cath::test::pdb_test_suite_fixture)

/// \brief Test the parsing of END separated pdbs by sticking every combination of 0-3 ENDS before, between and after two 1-atom PDBs
///        (except there must be at least one END between)
BOOST_AUTO_TEST_CASE(check_parsing_of_end_separators) {
	constexpr size_t MAX_NUM_CONSECUTIVE_ENDS =3;
	for (const size_t &num_before_ends : indices( MAX_NUM_CONSECUTIVE_ENDS ) ) {
		for (const size_t &num_between_ends : irange( 1_z, MAX_NUM_CONSECUTIVE_ENDS ) ) {
			for (const size_t &num_after_ends : indices( MAX_NUM_CONSECUTIVE_ENDS ) ) {
				stringstream end_separated_stream;
				end_separated_stream << join(str_vec(num_before_ends,  "END"), "\n") << "\n";
				end_separated_stream << "ATOM   2952  OXT ALA   385      70.681 -13.748  36.367  1.00 26.84           O\n";
				end_separated_stream << join(str_vec(num_between_ends, "END"), "\n") << "\n";
				end_separated_stream << "ATOM   2968  OXT ARG   387      20.593  77.271 -19.667  1.00  0.00           O\n";
				end_separated_stream << join(str_vec(num_after_ends,   "END"), "\n") << "\n";
				check_nums_of_atoms(
					read_end_separated_pdb_files(end_separated_stream),
					{ 1_z, 1_z }
				);
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(parses_mse_resiude_as_unk_x) {
	// Example is from chain B of 4c9a
	const string input_string = R"(ATOM   1861  N   THR B 138     -33.417  42.721 103.639  1.00142.96           N  
ATOM   1862  CA  THR B 138     -33.726  41.382 103.094  1.00124.05           C  
ATOM   1863  C   THR B 138     -34.247  41.593 101.682  1.00122.85           C  
ATOM   1864  O   THR B 138     -34.642  42.709 101.356  1.00126.34           O  
ATOM   1865  CB  THR B 138     -34.816  40.679 103.927  1.00108.88           C  
ATOM   1866  N   MSE B 139     -34.243  40.568 100.830  1.00114.60           N  
ATOM   1867  CA  MSE B 139     -34.924  40.699  99.523  1.00114.24           C  
ATOM   1868  C   MSE B 139     -36.170  39.811  99.512  1.00107.49           C  
ATOM   1869  O   MSE B 139     -36.514  39.225  98.475  1.00116.82           O  
ATOM   1870  CB  MSE B 139     -33.991  40.377  98.339  1.00105.83           C  
ATOM   1871  N   VAL B 140     -36.867  39.765 100.654  1.00 94.64           N  
ATOM   1872  CA  VAL B 140     -37.904  38.745 100.917  1.00101.80           C  
ATOM   1873  C   VAL B 140     -39.124  39.279 101.672  1.00100.30           C  
ATOM   1874  O   VAL B 140     -39.041  40.306 102.330  1.00 97.38           O  
ATOM   1875  CB  VAL B 140     -37.331  37.579 101.757  1.00 93.10           C)";
	istringstream input_ss{ input_string };
	const pdb the_pdb = read_pdb_file( input_ss );
	BOOST_REQUIRE_EQUAL( the_pdb.get_num_residues(), 3 );
	BOOST_CHECK_EQUAL( get_code_string( the_pdb.get_residue_of_index__backbone_unchecked( 1 ).get_amino_acid() ), "UNK");
	// \todo: Create a better way to refer to check for UNK in code/tests without a string literal
}


BOOST_AUTO_TEST_CASE(skips_and_warns_on_silly_amino_acid) {
	// Example is a modified version of chain B of 4c9a
	const string input_string = R"(ATOM   1861  N   THR B 138     -33.417  42.721 103.639  1.00142.96           N  
ATOM   1862  CA  THR B 138     -33.726  41.382 103.094  1.00124.05           C  
ATOM   1863  C   THR B 138     -34.247  41.593 101.682  1.00122.85           C  
ATOM   1864  O   THR B 138     -34.642  42.709 101.356  1.00126.34           O  
ATOM   1865  CB  THR B 138     -34.816  40.679 103.927  1.00108.88           C  
ATOM   1866  N   FOO B 139     -34.243  40.568 100.830  1.00114.60           N  
ATOM   1867  CA  FOO B 139     -34.924  40.699  99.523  1.00114.24           C  
ATOM   1868  C   FOO B 139     -36.170  39.811  99.512  1.00107.49           C  
ATOM   1869  O   FOO B 139     -36.514  39.225  98.475  1.00116.82           O  
ATOM   1870  CB  FOO B 139     -33.991  40.377  98.339  1.00105.83           C  
ATOM   1871  N   VAL B 140     -36.867  39.765 100.654  1.00 94.64           N  
ATOM   1872  CA  VAL B 140     -37.904  38.745 100.917  1.00101.80           C  
ATOM   1873  C   VAL B 140     -39.124  39.279 101.672  1.00100.30           C  
ATOM   1874  O   VAL B 140     -39.041  40.306 102.330  1.00 97.38           O  
ATOM   1875  CB  VAL B 140     -37.331  37.579 101.757  1.00 93.10           C)";
	ostringstream test_ss;
	istringstream input_ss{ input_string };
	const log_to_ostream_guard the_guard{ test_ss };

	const pdb the_pdb = read_pdb_file( input_ss );

	BOOST_REQUIRE_EQUAL( the_pdb.get_num_residues(), 2 );
	BOOST_CHECK( icontains( test_ss.str(), "skip" ) );
}

// Note: this test was originally written to check that such residues would be dropped
//       but we switched to a better system instead...
//
//       Now: If some residue has all alternative locations other than space or 'a', that increases tolerance
//       of that residue being absent in a corresponding DSSP file
BOOST_AUTO_TEST_CASE(allows_alt_loc_other_than_space_or_a) {
	// Example is from chain A of 3nir
	const string input_string = R"(ATOM    421  N  APRO A  22      -3.223 -12.679   4.864  0.70  1.90           N  
ATOM    422  CA APRO A  22      -2.719 -13.510   5.958  0.66  1.95           C  
ATOM    423  C  APRO A  22      -1.267 -13.212   6.318  0.63  1.82           C  
ATOM    424  O  APRO A  22      -0.409 -13.058   5.429  0.70  2.27           O  
ATOM    425  CB APRO A  22      -2.869 -14.951   5.453  0.62  3.01           C  
ATOM    426  CG APRO A  22      -4.087 -14.852   4.583  0.71  3.52           C  
ATOM    427  CD APRO A  22      -3.965 -13.515   3.897  0.60  3.10           C  
ATOM    428  HA APRO A  22      -3.484 -13.237   6.939  0.65  2.12           H  
ATOM    429  HB2APRO A  22      -2.005 -15.242   4.842  0.65  3.12           H  
ATOM    430  HB3APRO A  22      -3.029 -15.649   6.272  0.65  3.12           H  
ATOM    431  HG2APRO A  22      -4.115 -15.658   3.853  0.65  3.84           H  
ATOM    432  HG3APRO A  22      -4.996 -14.869   5.207  0.65  3.84           H  
ATOM    433  HD2APRO A  22      -3.389 -13.631   2.975  0.65  3.36           H  
ATOM    434  HD3APRO A  22      -4.941 -13.091   3.702  0.65  3.36           H  
ATOM    435  N  BSER A1022      -3.018 -12.669   4.880  0.36  2.28           N  
ATOM    436  CA BSER A1022      -2.587 -13.339   6.108  0.35  2.17           C  
ATOM    437  C  BSER A1022      -1.146 -12.973   6.446  0.30  2.37           C  
ATOM    438  O  BSER A1022      -0.307 -12.731   5.575  0.29  2.89           O  
ATOM    439  CB BSER A1022      -2.704 -14.834   5.765  0.28  2.84           C  
ATOM    440  OG BSER A1022      -1.766 -15.203   4.767  0.34  4.61           O  
ATOM    441  H  BSER A1022      -3.534 -13.194   4.157  0.35  2.20           H  
ATOM    442  HA BSER A1022      -3.031 -12.985   7.184  0.35  2.64           H  
ATOM    443  HB2BSER A1022      -2.478 -15.439   6.624  0.35  3.56           H  
ATOM    444  HB3BSER A1022      -3.676 -15.081   5.362  0.35  3.56           H  
ATOM    445  HG BSER A1022      -1.905 -14.683   3.926  0.35  6.84           H  
ATOM    446  N   GLU A  23      -0.939 -13.207   7.642  1.00  2.33           N  
ATOM    447  CA  GLU A  23       0.439 -13.004   8.080  1.00  2.21           C  
ATOM    448  C   GLU A  23       1.401 -13.972   7.382  1.00  2.08           C  
ATOM    449  O   GLU A  23       2.520 -13.584   7.062  1.00  2.46           O  
ATOM    450  CB  GLU A  23       0.597 -13.153   9.596  1.00  2.58           C  
ATOM    451  CG  GLU A  23      -0.085 -12.041  10.394  1.00  3.18           C  
ATOM    452  CD  GLU A  23       0.251 -12.089  11.902  1.00  4.41           C  
ATOM    453  OE1 GLU A  23       0.630 -13.127  12.424  1.00  5.47           O  
ATOM    454  OE2 GLU A  23       0.090 -11.025  12.579  1.00  6.76           O  
ATOM    455  H   GLU A  23      -1.608 -13.536   8.317  1.00  2.47           H  
ATOM    456  HA  GLU A  23       0.694 -12.123   7.844  1.00  2.32           H  
ATOM    457  HB2 GLU A  23       0.152 -14.118   9.897  1.00  2.72           H  
ATOM    458  HB3 GLU A  23       1.652 -13.179   9.840  1.00  2.72           H  
ATOM    459  HG2 GLU A  23       0.253 -11.081  10.015  1.00  3.48           H  
ATOM    460  HG3 GLU A  23      -1.163 -12.114  10.272  1.00  3.48           H  
ATOM    461  HE2 GLU A  23       1.108 -11.273  13.151  1.00  9.48           H)";
	ostringstream test_ss;
	istringstream input_ss{ input_string };
	const log_to_ostream_guard the_guard{ test_ss };

	const pdb the_pdb = read_pdb_file( input_ss );

	BOOST_CHECK_EQUAL( the_pdb.get_num_residues(), 3 );
}

BOOST_AUTO_TEST_CASE(accepts_different_aa_in_hetatm_records_without_warning) {
	const string input_string = R"(ATOM   1228  N   PHE A 165     -10.430  20.141  87.926  1.00 15.18           N  
ATOM   1229  CA  PHE A 165      -9.247  19.371  87.531  1.00 13.79           C  
ATOM   1230  C   PHE A 165      -9.033  18.252  88.537  1.00 12.81           C  
ATOM   1231  O   PHE A 165      -8.794  18.510  89.709  1.00 14.71           O  
HETATM 1239  O  AKPI A 166      -7.321  14.931  87.441  0.62 14.69           O  
HETATM 1240  N  AKPI A 166      -9.130  17.009  88.079  0.62 14.06           N  
HETATM 1241  CA AKPI A 166      -8.936  15.841  88.945  0.62 13.79           C  
HETATM 1242  CB AKPI A 166     -10.092  14.834  88.791  0.62 14.28           C  
HETATM 1243  CG AKPI A 166      -9.940  13.524  89.605  0.62 15.44           C  
HETATM 1244  CD AKPI A 166     -11.043  12.492  89.267  0.62 17.06           C  
HETATM 1245  CE AKPI A 166     -10.908  11.202  90.094  0.62 16.87           C  
HETATM 1246  NZ AKPI A 166     -11.829  10.079  89.676  0.62 17.67           N  
HETATM 1247  C  AKPI A 166      -7.608  15.179  88.604  0.62 14.38           C  
HETATM 1248  CX1AKPI A 166     -11.966   8.842  90.115  0.62 17.79           C  
HETATM 1249  C1 AKPI A 166     -11.040   8.345  91.186  0.62 21.96           C  
HETATM 1250  CX2AKPI A 166     -12.913   8.001  89.548  0.62 14.57           C  
HETATM 1251  O1 AKPI A 166     -13.275   6.876  90.123  0.62 16.54           O  
HETATM 1252  O2 AKPI A 166     -13.453   8.361  88.423  0.62 14.70           O  
ATOM   1253  N  BLYS A 166      -9.130  17.009  88.079  0.38 14.06           N  
ATOM   1254  CA BLYS A 166      -8.927  15.864  88.948  0.38 13.81           C  
ATOM   1255  C  BLYS A 166      -7.608  15.179  88.604  0.38 14.38           C  
ATOM   1256  O  BLYS A 166      -7.321  14.931  87.441  0.38 14.69           O  
ATOM   1257  CB BLYS A 166     -10.103  14.903  88.811  0.38 14.30           C  
ATOM   1258  CG BLYS A 166     -10.119  13.772  89.806  0.38 15.34           C  
ATOM   1259  CD BLYS A 166     -11.385  12.962  89.626  0.38 16.71           C  
ATOM   1260  CE BLYS A 166     -11.511  11.920  90.705  0.38 17.24           C  
ATOM   1261  NZ BLYS A 166     -11.704  12.533  92.029  0.38 16.94           N  
ATOM   1262  N   ASP A 167      -6.797  14.889  89.619  1.00 14.05           N  
ATOM   1263  CA  ASP A 167      -5.513  14.245  89.363  1.00 14.38           C  
ATOM   1264  C   ASP A 167      -5.469  12.788  89.749  1.00 14.39           C  
ATOM   1265  O   ASP A 167      -5.724  12.440  90.910  1.00 15.07           O  )";
	ostringstream test_ss;
	istringstream input_ss{ input_string };
	const log_to_ostream_guard the_guard{ test_ss };

	const pdb the_pdb = read_pdb_file( input_ss );

	BOOST_CHECK_EQUAL( the_pdb.get_num_residues(), 4  );
	BOOST_CHECK_EQUAL( test_ss.str(),              "" );
}

BOOST_AUTO_TEST_CASE(handle_chain_change_correctly) {
	const string input_string = R"(ATOM      1  N   LEU A   1      19.951  -0.078  26.341  1.00 74.36           N  
ATOM      2  CA  LEU A   1      20.671  -1.248  26.845  1.00 77.32           C  
ATOM      3  C   LEU A   1      20.314  -2.462  25.964  1.00 81.33           C  
ATOM      4  O   LEU A   1      19.481  -2.327  25.056  1.00 82.44           O  
TER       5      LEU A   1                                                      
HETATM    6  N   GLU H   2      41.172  36.192  51.348  1.00 17.41           N  
HETATM    7  CA  GLU H   2      40.491  35.926  50.091  1.00 15.21           C  
HETATM    8  C   GLU H   2      40.029  34.469  50.052  1.00 16.31           C  
HETATM    9  O   GLU H   2      40.630  33.666  50.805  1.00 17.29           O  
HETATM   10  C   ACT A   3      18.687   4.315  47.169  1.00 39.67           C  
END                                                                             
)";
	ostringstream test_ss;
	istringstream input_ss{ input_string };
	const log_to_ostream_guard the_guard{ test_ss };

	const pdb the_pdb = read_pdb_file( input_ss );

	BOOST_CHECK_EQUAL( the_pdb.get_num_residues(), 2  );
	BOOST_CHECK_EQUAL( test_ss.str(),              "" );
}

BOOST_AUTO_TEST_CASE(handles_clashing_residue_names) {
	ostringstream test_ss;
	const log_to_ostream_guard the_guard{ test_ss };

	// From 4tsw (as of ~April 2017 I think)
	// This has two different residues called 103 on chain B (and the same happens on C, D, E and F)
	const string pdb_data = R"(ATOM   2241  N   ARG B 103      30.036  20.481  30.582  1.00 67.41           N  
ATOM   2242  CA  ARG B 103      28.981  21.363  30.110  1.00 69.06           C  
ATOM   2243  C   ARG B 103      28.668  22.406  31.174  1.00 68.80           C  
ATOM   2244  O   ARG B 103      29.483  22.638  32.071  1.00 69.14           O  
ATOM   2245  CB  ARG B 103      27.725  20.557  29.768  1.00 71.39           C  
ATOM   2246  H   ARG B 103      29.950  20.114  31.491  1.00  0.00           H  
ATOM   2247  N   GLU B 103      27.492  23.017  31.036  1.00 69.00           N  
ATOM   2248  CA  GLU B 103      26.946  24.055  31.918  1.00 69.65           C  
ATOM   2249  C   GLU B 103      26.058  24.961  31.066  1.00 71.54           C  
ATOM   2250  O   GLU B 103      26.334  25.177  29.877  1.00 72.21           O  
ATOM   2251  CB  GLU B 103      28.051  24.882  32.590  1.00 67.34           C  
ATOM   2252  H   GLU B 103      26.896  22.772  30.305  1.00  0.00           H  
TER    2253      GLU B 103                                                      
END   
)";
	BOOST_TEST( pdb_file_to_string( read_pdb( pdb_data ) ) == pdb_data );
	BOOST_CHECK( regex_search( test_ss.str(), regex{ R"(found conflicting consecutive entries for residue "B:103")" } ) );
}

BOOST_AUTO_TEST_CASE(handles_very_large_temperature_factor_number_a) {
	// From 2pyi (as of ~April 2017 I think)
	// Atom 6358 has a very large temperature factor that can overflow the field
	const string pdb_data = R"(ATOM   6351  N   ILE A 805      15.633  14.745  13.895  1.00 14.96           N  
ATOM   6352  CA  ILE A 805      15.277  13.957  15.076  1.00 14.96           C  
ATOM   6353  C   ILE A 805      15.671  12.494  14.877  1.00 14.94           C  
ATOM   6354  O   ILE A 805      16.311  11.888  15.742  1.00 14.77           O  
ATOM   6355  CB  ILE A 805      13.774  14.071  15.424  1.00 14.94           C  
ATOM   6356  CG1 ILE A 805      13.408  15.522  15.764  1.00 15.09           C  
ATOM   6357  CG2 ILE A 805      13.423  13.123  16.576  1.00 15.41           C  
ATOM   6358  CD1 ILE A 805      11.913  15.784  15.849  1.001511.0           C  
TER    6359      ILE A 805                                                      
END   
)";
	BOOST_TEST( pdb_file_to_string( read_pdb( pdb_data ) ) == pdb_data );
}

BOOST_AUTO_TEST_CASE(handles_very_large_temperature_factor_number_b) {
	// From 2qlm (as of ~April 2017 I think)
	// Atom 850 has a very large temperature factor that can overflow the field
	const string pdb_data = R"(ATOM    850  N   TYR A 113      -2.998  15.281  42.352  1.003106.0           N  
ATOM    851  CA  TYR A 113      -4.266  15.757  41.788  1.00 31.85           C  
ATOM    852  C   TYR A 113      -5.440  15.471  42.728  1.00 31.76           C  
ATOM    853  O   TYR A 113      -6.257  16.360  42.988  1.00 31.50           O  
ATOM    854  CB  TYR A 113      -4.519  15.135  40.405  1.00 32.84           C  
ATOM    855  CG  TYR A 113      -5.862  15.488  39.790  1.00 34.00           C  
ATOM    856  CD1 TYR A 113      -6.090  16.749  39.235  1.00 35.44           C  
ATOM    857  CD2 TYR A 113      -6.903  14.559  39.763  1.00 35.50           C  
ATOM    858  CE1 TYR A 113      -7.322  17.078  38.670  1.00 35.89           C  
ATOM    859  CE2 TYR A 113      -8.143  14.876  39.200  1.00 36.07           C  
ATOM    860  CZ  TYR A 113      -8.342  16.137  38.654  1.00 35.90           C  
ATOM    861  OH  TYR A 113      -9.562  16.460  38.098  1.00 36.20           O  
TER     862      TYR A 113                                                      
END   
)";
	BOOST_TEST( pdb_file_to_string( read_pdb( pdb_data ) ) == pdb_data );
}

BOOST_AUTO_TEST_CASE(get_backbone_complete_residue_ids__removes_dupl_res_ids) {
	ostringstream test_ss;
	const log_to_ostream_guard the_guard{ test_ss };

	// From 4tsw (as of ~January 2018)
	//
	// Ensure that get_backbone_complete_residue_ids() removes duplicate residue_ids
	// so that it matches the behaviour of making proteins from PDBs.
	// If necessary, it'd be possible to make both try to handle them but
	// that will likely require work to fix arising problems and the issue
	// of duplicated residue is rare and is (mostly? completely?) restricted
	// to superseded PDBs
	const string pdb_data = R"(ATOM   2229  N   GLN B 102      31.966  19.126  31.868  1.00 57.91           N  
ATOM   2230  CA  GLN B 102      32.135  19.286  30.428  1.00 63.27           C  
ATOM   2231  C   GLN B 102      31.098  20.226  29.827  1.00 65.08           C  
ATOM   2232  O   GLN B 102      31.295  20.773  28.738  1.00 66.67           O  
ATOM   2241  N   ARG B 103      30.036  20.481  30.582  1.00 67.41           N  
ATOM   2242  CA  ARG B 103      28.981  21.363  30.110  1.00 69.06           C  
ATOM   2243  C   ARG B 103      28.668  22.406  31.174  1.00 68.80           C  
ATOM   2244  O   ARG B 103      29.483  22.638  32.071  1.00 69.14           O  
ATOM   2247  N   GLU B 103      27.492  23.017  31.036  1.00 69.00           N  
ATOM   2248  CA  GLU B 103      26.946  24.055  31.918  1.00 69.65           C  
ATOM   2249  C   GLU B 103      26.058  24.961  31.066  1.00 71.54           C  
ATOM   2250  O   GLU B 103      26.334  25.177  29.877  1.00 72.21           O  
ATOM   2253  N   THR B 105      24.903  25.407  31.312  1.00 71.90           N  
ATOM   2254  CA  THR B 105      24.124  26.292  30.435  1.00 69.49           C  
ATOM   2255  C   THR B 105      23.000  27.081  31.135  1.00 67.87           C  
ATOM   2256  O   THR B 105      22.911  27.101  32.373  1.00 63.94           O  
TER    2257      THR B 105                                                      
END   
)";
	const residue_id_vec expected = { {
		make_residue_id( 'B', 102 ),
		make_residue_id( 'B', 103 ),
		make_residue_id( 'B', 105 ),
	} };
	BOOST_TEST( get_backbone_complete_residue_ids( read_pdb( pdb_data ) ) == expected );
}

BOOST_AUTO_TEST_CASE(writes_partial_pdb_correctly) {
	const auto parsed_pdb = read_pdb_file( global_test_constants::EXAMPLE_A_PDB_FILENAME() );

	const temp_file temp_test{ ".pdb_test_output.%%%%-%%%%-%%%%-%%%%" };
	const path temp_test_file = get_filename( temp_test );

	write_pdb_file(
		temp_test_file,
		parsed_pdb,
		{ { make_simple_region( 'A', 1321, 1324 ) } }
	);

	const path expected{ TEST_EXAMPLE_PDBS_DATA_DIR() / "1c0pA01.1321-1324:A" };
	BOOST_CHECK_FILES_EQUAL( temp_test_file, expected );
}

BOOST_AUTO_TEST_SUITE_END()
