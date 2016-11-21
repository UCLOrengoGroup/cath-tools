/// \file
/// \brief The bifur_hbond_list class definitions

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

#include "bifur_hbond_list.h"

using namespace cath::sec;
// using namespace cath::file;

// double bifur_hbond_list::h_bond_energy_of_full_non_proline_residues(const pdb_residue     &arg_residue_i,          ///< TODOCUMENT
//                                                                const pdb_residue_opt &arg_residue_j_prev_opt, ///< TODOCUMENT
//                                                                const pdb_residue     &arg_residue_j           ///< TODOCUMENT
//                                                                ) {
// 	const auto i_c         = arg_residue_i.get_carbon_coord_of_residue();
// 	const auto i_o         = arg_residue_i.get_oxygen_coord_of_residue();

// 	const auto j_n         = arg_residue_j.get_nitrogen_coord_of_residue();
// 	const auto j_o         = j_narg_residue_j.get_oxygen_coord_of_residue();
// 	const bool do_h_offset = ( arg_residue_j_prev_opt && )
// 	if ( arg_residue_j_prev_opt && )
// 	return false;
// 	arg_residue_i;
// 	arg_residue_j_prev;
// 	arg_residue_j;
// }

// BOOST_AUTO_TEST_CASE(basic) {
// 	parse_dssp_for_calc_testing( "/cath-tools/build-test-data/1c0pA01.dssp" );


// 	geom::coord get_atom_of_id_of_residue(const pdb_residue &,
// 	                                      const std::string &);
// 	geom::coord get_nitrogen_coord_of_residue(const pdb_residue &);
// 	geom::coord get_carbon_alpha_coord_of_residue(const pdb_residue &);
// 	geom::coord get_carbon_coord_of_residue(const pdb_residue &);
// 	geom::coord get_carbon_beta_coord_of_residue(const pdb_residue &);
// 	geom::coord get_oxygen_coord_of_residue(const pdb_residue &);

// 	BOOST_CHECK( true );
// 	const coord      c{  22.093, 133.611,  36.006 };
// 	const coord      o{  21.020, 133.627,  36.607 };
// 	const coord      n{  18.386, 132.630,  37.199 };

// 	const coord prev_c{  18.001, 131.369,  37.293 };
// 	const coord prev_o{  16.908, 131.048,  37.733 };
// 	const coord prev_c_to_o = prev_o - prev_c;
// 	const coord      h = n - ( prev_c_to_o / length( prev_c_to_o ) );
// 	// const coord      h { 19.281, 132.893,  36.8387 };
	
// 	// ATOM    149  N   SER A1019      19.532 130.693  35.496  1.00  8.91
// 	// ATOM    151  C   SER A1019      18.001 131.369  37.293  1.00  8.78
// 	// ATOM    152  O   SER A1019      16.908 131.048  37.733  1.00  9.57
// 	// ATOM    153  CB  SER A1019      20.136 130.128  37.776  1.00  9.19
// 	// ATOM    154  OG  SER A1019      21.016 129.101  37.430  1.00 10.29

// 	const double dist_on = distance_between_points( o, n );
// 	const double dist_ch = distance_between_points( c, h );
// 	const double dist_oh = distance_between_points( o, h );
// 	const double dist_cn = distance_between_points( c, n );
// 	const double inv_dist_on = 1.0 / dist_on;
// 	const double inv_dist_ch = 1.0 / dist_ch;
// 	const double inv_dist_oh = 1.0 / dist_oh;
// 	const double inv_dist_cn = 1.0 / dist_cn;

// 	const double energy = 0.42 * 0.2 * 332 * ( inv_dist_on + inv_dist_ch - inv_dist_oh - inv_dist_cn );

// 	// Created residue ALA with hydrogen 19.281, 132.893, 36.8387

// 	std::cerr << "dist_c      is : " << c           << "\n";
// 	std::cerr << "dist_o      is : " << o           << "\n";
// 	std::cerr << "dist_n      is : " << n           << "\n";
// 	std::cerr << "dist_h      is : " << h           << "\n";
// 	std::cerr << "dist_on     is : " << dist_on     << "\n";
// 	std::cerr << "dist_ch     is : " << dist_ch     << "\n";
// 	std::cerr << "dist_oh     is : " << dist_oh     << "\n";
// 	std::cerr << "dist_cn     is : " << dist_cn     << "\n";
// 	std::cerr << "inv_dist_on is : " << inv_dist_on << "\n";
// 	std::cerr << "inv_dist_ch is : " << inv_dist_ch << "\n";
// 	std::cerr << "inv_dist_oh is : " << inv_dist_oh << "\n";
// 	std::cerr << "inv_dist_cn is : " << inv_dist_cn << "\n";

// 	std::cerr << "energy is        : " << energy          << "\n";
// 	std::cerr << "should be        : " << -0.968158347676 << "\n";

// 	// ATOM    133  C   GLY A1016      22.093 133.611  36.006  1.00  8.83
// 	// ATOM    134  O   GLY A1016      21.020 133.627  36.607  1.00  9.33

// 	// ATOM    155  N   ALA A1020      18.386 132.630  37.199  1.00  9.12
// 	// ATOM    157  C   ALA A1020      16.225 133.709  36.827  1.00  9.61
// 	// ATOM    158  O   ALA A1020      15.148 133.956  37.406  1.00 10.41
// }