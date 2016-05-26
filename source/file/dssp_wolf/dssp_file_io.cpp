/// \file
/// \brief The dssp_file_io definitions

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

#include "dssp_file_io.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>

#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "exception/runtime_error_exception.h"
#include "file/dssp_wolf/dssp_file.h"
#include "structure/chain_label.h"
#include "structure/protein/residue.h"

#include <cmath>
#include <fstream>

using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::algorithm::is_space;
using boost::algorithm::starts_with;
using boost::algorithm::trim_copy;
using boost::lexical_cast;

/// \brief Parse a dssp_file object from an istream
///
/// \relates dssp_file
dssp_file cath::file::read_dssp_file(const path &arg_dssp_file ///< The DSSP file from which to parse a dssp_file object
                                     ) {
	ifstream my_dssp_istream;
	open_ifstream( my_dssp_istream, arg_dssp_file );
	const dssp_file the_dssp_file = read_dssp( my_dssp_istream );
	my_dssp_istream.close();
	return the_dssp_file;
}

/// \brief Parse a dssp_file object from an istream
///
/// \relates dssp_file
dssp_file cath::file::read_dssp(istream &arg_istream ///< The istream from which to parse a dssp_file object
                                ) {
	// Read the first line and check it looks vaguely like a DSSP header line
	string dssp_first_line;
	getline(arg_istream, dssp_first_line);
	if ( !starts_with(dssp_first_line, "==== ") || !icontains(dssp_first_line, "DSSP") || dssp_first_line.length() > 300 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("First line is not a DSSP header line"));
	}

	// Scan through the header
	bool found_column_headings( false );
	while ( ! found_column_headings ) {
		string dssp_header_line;
		getline(arg_istream, dssp_header_line);
		found_column_headings = starts_with(dssp_header_line, "  #  RESIDUE AA ");
	}

	// Read residues
	size_t residue_ctr = 1;
	residue_vec new_residues;
	while ( ! arg_istream.eof()) {
		string dssp_residue_line;
		getline(arg_istream, dssp_residue_line);

		// If this line contains nothing but whitespace, then skip it
		if (all(dssp_residue_line, is_space())) {
			continue;
		}

		const bool dssp_entry_is_null = ( dssp_residue_line.at(13) == '!' );
		// If this line indicates a problem residue, then increment the residue counter and skip it
//		if (dssp_entry_is_null) {
////			cerr << "Chain ? : " << residue::null_residue << endl;
//			new_residues.push_back(residue::NULL_RESIDUE);
//			++residue_ctr;
//			continue;
//		}

//		cerr << "Parse residue " << residue_ctr << " from \"" << dssp_residue_line << "\"" << endl;
		const size_chain_residue_tuple  residue_details      = parse_dssp_residue_line( dssp_residue_line );
		const size_t                   &parsed_residue_index = get<0>( residue_details );
		const chain_label              &the_chain_label      = get<1>( residue_details );
		const residue                  &parsed_residue       = get<2>( residue_details );

		if ( ! dssp_entry_is_null && parsed_residue_index != residue_ctr ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Error in DSSP sequential residue numbers"));
		}
//		cerr << "Chain " << the_chain_label << " : " << parsed_residue << endl;

		// Some PDBs (eg 4tsw) may have erroneous consecutive duplicate residues.
		// Though that's a bit rubbish, it shouldn't break the whole comparison
		// so if that's detected, just warn and move on (without appending to new_residues).
		if ( ! dssp_entry_is_null && ! new_residues.empty() && new_residues.back().get_pdb_residue_name() == parsed_residue.get_pdb_residue_name() ) {
			BOOST_LOG_TRIVIAL( warning ) << "Whilst parsing DSSP file, found conflicting consecutive entries for residue \""
				<< parsed_residue.get_pdb_residue_name()
				<< "\" on chain '"
				<< the_chain_label
				<< "' (with amino acids \""
				<< new_residues.back().get_amino_acid().get_code()
				<< "\" and then \""
				<< parsed_residue.get_amino_acid().get_code()
				<< "\") - ignoring latter entry";
		}
		else {
			new_residues.push_back( parsed_residue );
		}
		++residue_ctr;
	}

	return dssp_file( new_residues );
}

/// \brief TODOCUMENT
size_chain_residue_tuple cath::file::parse_dssp_residue_line(const string &arg_dssp_residue_line ///< The DSSP residue line to parse
                                                             ) {
	try {
		const bool dssp_entry_is_null = ( arg_dssp_residue_line.at( 13 ) == '!' );
		if ( dssp_entry_is_null ) {
			return make_tuple( 0_z, chain_label( ' ' ), residue::NULL_RESIDUE );
		}
		                                                                                         // Comments with DSSP format documentation
		                                                                                         // (http://swift.cmbi.ru.nl/gv/dssp/)
		const size_t        seq_res_num    =     stoul( arg_dssp_residue_line.substr(  0, 5)) ;  //   1 -   5    sequential resnumber, including chain breaks as extra residues
		const string        res_name_str   = trim_copy( arg_dssp_residue_line.substr(  5, 6)) ;  //   6 -  11    original PDB resname, not nec. sequential, may contain letters
		const char         &chain_char     =            arg_dssp_residue_line.at    ( 11   )  ;  //  12
		const char         &amino_acid_c   =            arg_dssp_residue_line.at    ( 13   )  ;  //  14          amino acid sequence in one letter code
		const char         &sstruc_code    =            arg_dssp_residue_line.at    ( 16   )  ;  //  17          secondary structure summary based on columns 19-38
//		const char         &turn3_helix    =            arg_dssp_residue_line.at    ( 18   )  ;  //  19          3-turns/helix
//		const char         &turn4_helix    =            arg_dssp_residue_line.at    ( 19   )  ;  //  20          4-turns/helix
//		const char         &turn5_helix    =            arg_dssp_residue_line.at    ( 20   )  ;  //  21          5-turns/helix
//		const char         &geom_bend      =            arg_dssp_residue_line.at    ( 21   )  ;  //  22          geometrical bend
//		const char         &chirality      =            arg_dssp_residue_line.at    ( 22   )  ;  //  23          chirality
//		const char         &beta_brid_1    =            arg_dssp_residue_line.at    ( 23   )  ;  //  24          beta bridge label
//		const char         &beta_brid_2    =            arg_dssp_residue_line.at    ( 24   )  ;  //  25          beta bridge label
//		const size_t        beta_part_1    =     stoul( arg_dssp_residue_line.substr( 25, 4)) ;  //  26 -  29    beta bridge partner resnum
//		const size_t        beta_part_2    =     stoul( arg_dssp_residue_line.substr( 29, 4)) ;  //  30 -  33    beta bridge partner resnum
//		const char         &sheet_label    =            arg_dssp_residue_line.at    ( 33   )  ;  //  34          beta sheet label
		const size_t        solv_access    =     stoul( arg_dssp_residue_line.substr( 34, 4)) ;  //  35 -  38    solvent accessibility
		// There are other four other columns here (N-H-->O, O-->H-N, N-H-->O and O-->H-N)
//		const double        tco            =      stod( arg_dssp_residue_line.substr( 85, 6)) ;  //  86 -  91    cosine of angle between C=O of residue I and C=O of residue I-1. For α-helices, TCO is near +1, for β-sheets TCO is near -1. Not used for structure definition.
//		const double        kappa          =      stod( arg_dssp_residue_line.substr( 91, 6)) ;  //  92 -  97    virtual bond angle (bend angle) defined by the three Cα atoms of residues I-2,I,I+2. Used to define bend (structure code 'S').
//		const double        alpha          =      stod( arg_dssp_residue_line.substr( 97, 6)) ;  //  98 - 103    virtual torsion angle (dihedral angle) defined by the four Cα atoms of residues I-1,I,I+1,I+2.Used to define chirality (structure code '+' or '-').
		const double        phi_in_degrees =      stod( arg_dssp_residue_line.substr(103, 6)) ;  // 104 - 109    IUPAC peptide backbone torsion angles
		const double        psi_in_degrees =      stod( arg_dssp_residue_line.substr(109, 6)) ;  // 110 - 115    IUPAC peptide backbone torsion angles
		const double        carbon_a_x     =      stod( arg_dssp_residue_line.substr(115, 7)) ;  // 116 - 122    echo of Ca atom coordinates
		const double        carbon_a_y     =      stod( arg_dssp_residue_line.substr(122, 7)) ;  // 123 - 129    echo of Ca atom coordinates
		const double        carbon_a_z     =      stod( arg_dssp_residue_line.substr(129, 7)) ;  // 130 - 136    echo of Ca atom coordinates

		// Parse the residue name string into a residue_name
		const residue_name res_name = make_residue_name( res_name_str );

		// Convert the DSSP secondary structure summary code into a sec_struc_type
		const sec_struc_type sec_struc = sstruc_code == 'H' ? sec_struc_type::ALPHA_HELIX
		                               : sstruc_code == 'E' ? sec_struc_type::BETA_STRAND
		                                                    : sec_struc_type::COIL;

//		cerr << "seq_res_num    : " << seq_res_num    << endl;
//		cerr << "res_name       : " << res_name       << endl;
//		cerr << "chain_char     : " << chain_char     << endl;
//		cerr << "amino_acid     : " << amino_acid     << endl;
//		cerr << "sec_struc      : " << sec_struc      << endl;
//		cerr << "turn3_helix    : " << turn3_helix    << endl;
//		cerr << "turn4_helix    : " << turn4_helix    << endl;
//		cerr << "turn5_helix    : " << turn5_helix    << endl;
//		cerr << "geom_bend      : " << geom_bend      << endl;
//		cerr << "chirality      : " << chirality      << endl;
//		cerr << "beta_brid_1    : " << beta_brid_1    << endl;
//		cerr << "beta_brid_2    : " << beta_brid_2    << endl;
//		cerr << "beta_part_1    : " << beta_part_1    << endl;
//		cerr << "beta_part_2    : " << beta_part_2    << endl;
//		cerr << "sheet_label    : " << sheet_label    << endl;
//		cerr << "solv_access    : " << solv_access    << endl;
//		cerr << "tco            : " << tco            << endl;
//		cerr << "kappa          : " << kappa          << endl;
//		cerr << "alpha          : " << alpha          << endl;
//		cerr << "phi_in_degrees : " << phi_in_degrees << endl;
//		cerr << "psi_in_degrees : " << psi_in_degrees << endl;
//		cerr << "carbon_a_x     : " << carbon_a_x     << endl;
//		cerr << "carbon_a_y     : " << carbon_a_y     << endl;
//		cerr << "carbon_a_z     : " << carbon_a_z     << endl;

		const auto shifted_phi = make_angle_from_degrees<double>( round( phi_in_degrees > 0 ? phi_in_degrees : 360.0 + phi_in_degrees) );
		const auto shifted_psi = make_angle_from_degrees<double>( round( psi_in_degrees > 0 ? psi_in_degrees : 360.0 + psi_in_degrees) );

		return make_tuple(
			seq_res_num,
			chain_label( chain_char ),
			residue(
				res_name,
				amino_acid(amino_acid_c),
				coord( carbon_a_x, carbon_a_y, carbon_a_z ),
				coord::ORIGIN_COORD,
				0,
				sec_struc,
				rotation::IDENTITY_ROTATION(),
				shifted_phi,
				shifted_psi,
				solv_access
			)
		);
	}
	catch ( const boost::bad_lexical_cast & ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to cast a column whilst parsing a DSSP residue record.\nRecord was \"" + arg_dssp_residue_line + "\""));
	}
}
