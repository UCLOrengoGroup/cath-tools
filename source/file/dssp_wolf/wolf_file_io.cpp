/// \file
/// \brief The wolf_file_io class definitions

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

#include "wolf_file_io.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "exception/runtime_error_exception.hpp"
#include "file/dssp_wolf/wolf_file.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/rotation.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc_type.hpp"

#include <cstdio>
#include <iostream>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::file;
using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::numeric_cast;

/// \brief TODOCUMENT
///
/// \relates protein
wolf_file cath::file::read_wolf(const path &arg_wolf_filename ///< TODOCUMENT
                                ) {
	const double WOLF_ROTATION_TOLERANCE(0.00021);

	char   c;
	char   residue_string[6];
	char   pdb_name_insert;

	// Open residue data file
	FILE *wolf_infile;
	if ( ( wolf_infile = fopen(arg_wolf_filename.string().c_str(),"r") ) == nullptr )  {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to open wolf file \"" + arg_wolf_filename.string() + "\" for reading"));
	}

	size_t length = 0;
	bool foundlength = false;
	while ( !foundlength ) {
		c = numeric_cast<char>( getc( wolf_infile ) );
		if ( c == EOF ) {
			cerr << "ssap: WOLF file ended unexpectedly!" << endl;
			fclose( wolf_infile );
			exit(1);
		}

		if ( c != ' ' ) {
			while( c = numeric_cast<char>( getc( wolf_infile ) ), c != '\n' ) {
				if (c==EOF) {
					cerr << "ssap: WOLF file ended unexpectedly!" << endl;
					fclose( wolf_infile );
					exit(1);
				}
			}
		}
		else {
			if (!fscanf( wolf_infile, "%zu", &length )) {
				cerr << "ssap: WOLF file ended unexpectedly!" << endl;
				fclose( wolf_infile );
				exit(1);
			}
			foundlength = true;
		}
	}

	// read total length of protein data in file and allocate memory
	while ( c = numeric_cast<char>( getc( wolf_infile ) ), c!= '#'  );
	while ( c = numeric_cast<char>( getc( wolf_infile ) ), c!= '\n' );

	residue_vec residues;
	residues.reserve(length);
	// Read through file saving data for selected chain
	for (size_t residue_ctr = 0; residue_ctr < length; ++residue_ctr) {

		int       number; // The number of the residue in the WOLF file perhaps?
		int       pdb_name_number;
		char      buffer[1000];

		// Read Kabsch number and residue number in structure databank
		buffer[0] = '\n';
		buffer[1] = 0;

		if (!fgets(buffer, 999, wolf_infile)) {
			cerr << "WOLF file ended unexpectedly!" << endl;
			fclose( wolf_infile );
			exit(1);
		}

		//cerr << "About to sscanf the following string : \"" << buffer << "\"" << endl;
		sscanf(buffer,"%d %s", &number, residue_string);
		pdb_name_insert = ' ';
		sscanf(residue_string,"%d%c", &pdb_name_number, &pdb_name_insert);

		char chain_char;
		sscanf( buffer+11, "%c", &chain_char );

		// Read amino acid
		char amino_acid_char;
		sscanf(buffer+13, "%c", &amino_acid_char);
		if (amino_acid_char == '!') {
			continue;
		}

		// Save disulphide bonds
		if ( islower(amino_acid_char) ) {
			cerr << "IMPORTANT NOTICE : ******** THIS BIT OF CODE IS NOT COMPLETELY USELESS - HURRAH - PLEASE TELL SOMEONE *******" << endl;
			amino_acid_char ='C';
		}

		int    access;
		int    dummy_hydrogen_bond_int;
		double dummy_hydrogen_bond_float;
		double phi_in_degrees;
		double psi_in_degrees;
		double dummy_alpha;
		double dummy_kappa;
		double dummy_tco;
		double ca_x, ca_y, ca_z;
		double cb_x, cb_y, cb_z;
		double frame_x1, frame_y1, frame_z1;
		double frame_x2, frame_y2, frame_z2;
		double frame_x3, frame_y3, frame_z3;

		sscanf(
			buffer+34,
			"%d\
			%d,%lg%d,%lg %d,%lg%d,%lg\
			%lg%lg%lg%lg%lg %lg%lg%lg %lg%lg%lg\
			%d%lg%d%lg %d%lg%d%lg %d%lg%d%lg\
			%lg%lg%lg %lg%lg%lg %lg%lg%lg",

			// Solvent accessibility
			&access,

			// Hydrogen bond information
			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,
			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,

			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,
			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,

			// Torsional and backbone angles
			&dummy_tco, &dummy_kappa, &dummy_alpha, &phi_in_degrees, &psi_in_degrees,

			// Coordinates of Ca and Cb atoms of residues
			&ca_x, &ca_y, &ca_z,
			&cb_x, &cb_y, &cb_z,

			// Hydrogen bond information (3 extra pairs over Kabsch & Sanders data)
			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,
			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,

			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,
			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,

			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,
			&dummy_hydrogen_bond_int, &dummy_hydrogen_bond_float,

			//  Frames of reference for residue
			&frame_x1, &frame_y1, &frame_z1,
			&frame_x2, &frame_y2, &frame_z2,
			&frame_x3, &frame_y3, &frame_z3
		);

//		const rotation frame_rotation = tidy_rotation(
//			frame_x1, frame_x2, frame_x3,
//			frame_y1, frame_y2, frame_y3,
//			frame_z1, frame_z2, frame_z3,
//			WOLF_ROTATION_TOLERANCE
//		);
//		const rotation frame_rotation = transpose_copy(tidy_rotation(
//			frame_x1, frame_y1, frame_z1,
//			frame_x2, frame_y2, frame_z2,
//			frame_x3, frame_y3, frame_z3,
//			WOLF_ROTATION_TOLERANCE
//		));
		const rotation frame_rotation(
			frame_x1, frame_y1, frame_z1,
			frame_x2, frame_y2, frame_z2,
			frame_x3, frame_y3, frame_z3,
			WOLF_ROTATION_TOLERANCE
		);

		const coord carbon_alpha_coord ( ca_x, ca_y, ca_z );
		const coord carbon_beta_coord  ( cb_x, cb_y, cb_z );

		// Set backbone and torsional angles
		const auto shifted_phi = shift_copy( make_angle_from_degrees<double>( numeric_cast<int>( phi_in_degrees ) ), zero_angle<double>(), angle_endpoint_loc::USE_EITHER );
		const auto shifted_psi = shift_copy( make_angle_from_degrees<double>( numeric_cast<int>( psi_in_degrees ) ), zero_angle<double>(), angle_endpoint_loc::USE_EITHER );
//		dummy_alpha = (dummy_alpha <   0) ? (360 + dummy_alpha) : dummy_alpha;
//		dummy_kappa = (dummy_kappa <   0) ? (360 + dummy_kappa) : dummy_kappa;
//		dummy_kappa = (dummy_kappa > 180) ? (360 - dummy_kappa) : dummy_kappa;

		const residue_id res_id{
			chain_label{ chain_char },
			make_residue_name_with_non_insert_char( pdb_name_number, pdb_name_insert, ' ' )
		};

		const int sec_struc_number(0);
		const sec_struc_type the_sec_struc_type(sec_struc_type::COIL);
		residues.push_back(residue(
			res_id,
			amino_acid( amino_acid_char ),
			carbon_alpha_coord,
			carbon_beta_coord,

			sec_struc_number,
			the_sec_struc_type,

			frame_rotation,
			shifted_phi,
			shifted_psi,
			numeric_cast<size_t>( access )
		));
	}

	fclose( wolf_infile );

	return wolf_file(residues);
}
