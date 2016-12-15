ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   1) /// \file
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   2) /// \brief The dssp_file_io definitions
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   3) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   4) /// \copyright
7790957b (Tony Lewis 2015-07-27 12:27:10 +0100   5) /// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   6) /// Copyright (C) 2011, Orengo Group, University College London
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   7) ///
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   8) /// This program is free software: you can redistribute it and/or modify
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100   9) /// it under the terms of the GNU General Public License as published by
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  10) /// the Free Software Foundation, either version 3 of the License, or
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  11) /// (at your option) any later version.
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  12) ///
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  13) /// This program is distributed in the hope that it will be useful,
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  14) /// but WITHOUT ANY WARRANTY; without even the implied warranty of
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  15) /// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  16) /// GNU General Public License for more details.
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  17) ///
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  18) /// You should have received a copy of the GNU General Public License
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  19) /// along with this program.  If not, see <http://www.gnu.org/licenses/>.
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  20) 
13b40c2b (Tony Lewis 2016-12-02 17:31:51 +0000  21) #include "dssp_file_io.hpp"
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  22) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  23) #include <boost/algorithm/string/classification.hpp>
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  24) #include <boost/algorithm/string/predicate.hpp>
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  25) #include <boost/algorithm/string/trim.hpp>
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  26) #include <boost/lexical_cast.hpp>
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100  27) #include <boost/log/trivial.hpp>
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  28) 
13b40c2b (Tony Lewis 2016-12-02 17:31:51 +0000  29) #include "common/file/open_fstream.hpp"
13b40c2b (Tony Lewis 2016-12-02 17:31:51 +0000  30) #include "common/size_t_literal.hpp"
13b40c2b (Tony Lewis 2016-12-02 17:31:51 +0000  31) #include "exception/runtime_error_exception.hpp"
13b40c2b (Tony Lewis 2016-12-02 17:31:51 +0000  32) #include "file/dssp_wolf/dssp_file.hpp"
13b40c2b (Tony Lewis 2016-12-02 17:31:51 +0000  33) #include "structure/chain_label.hpp"
13b40c2b (Tony Lewis 2016-12-02 17:31:51 +0000  34) #include "structure/protein/residue.hpp"
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  35) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  36) #include <cmath>
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  37) #include <fstream>
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  38) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  39) using namespace boost::algorithm;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  40) using namespace boost::filesystem;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  41) using namespace cath;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  42) using namespace cath::common;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  43) using namespace cath::file;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  44) using namespace cath::geom;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  45) using namespace std;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  46) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  47) using boost::algorithm::is_space;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  48) using boost::algorithm::starts_with;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  49) using boost::algorithm::trim_copy;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  50) using boost::lexical_cast;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  51) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  52) /// \brief Parse a dssp_file object from an istream
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  53) ///
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  54) /// \relates dssp_file
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  55) dssp_file cath::file::read_dssp_file(const path &arg_dssp_file ///< The DSSP file from which to parse a dssp_file object
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  56)                                      ) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  57) 	ifstream my_dssp_istream;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  58) 	open_ifstream( my_dssp_istream, arg_dssp_file );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  59) 	const dssp_file the_dssp_file = read_dssp( my_dssp_istream );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  60) 	my_dssp_istream.close();
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  61) 	return the_dssp_file;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  62) }
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  63) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  64) /// \brief Parse a dssp_file object from an istream
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  65) ///
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  66) /// \relates dssp_file
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  67) dssp_file cath::file::read_dssp(istream &arg_istream ///< The istream from which to parse a dssp_file object
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  68)                                 ) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  69) 	// Read the first line and check it looks vaguely like a DSSP header line
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  70) 	string dssp_first_line;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  71) 	getline(arg_istream, dssp_first_line);
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  72) 	if ( !starts_with(dssp_first_line, "==== ") || !icontains(dssp_first_line, "DSSP") || dssp_first_line.length() > 300 ) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  73) 		BOOST_THROW_EXCEPTION(runtime_error_exception("First line is not a DSSP header line"));
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  74) 	}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  75) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  76) 	// Scan through the header
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  77) 	bool found_column_headings( false );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  78) 	while ( ! found_column_headings ) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  79) 		string dssp_header_line;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  80) 		getline(arg_istream, dssp_header_line);
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  81) 		found_column_headings = starts_with(dssp_header_line, "  #  RESIDUE AA ");
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  82) 	}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  83) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  84) 	// Read residues
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  85) 	size_t residue_ctr = 1;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  86) 	residue_vec new_residues;
8c49d7aa (Tony Lewis 2016-05-26 17:59:20 +0100  87) 	while ( ! arg_istream.eof()) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  88) 		string dssp_residue_line;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  89) 		getline(arg_istream, dssp_residue_line);
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  90) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  91) 		// If this line contains nothing but whitespace, then skip it
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  92) 		if (all(dssp_residue_line, is_space())) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  93) 			continue;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  94) 		}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  95) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  96) 		const bool dssp_entry_is_null = ( dssp_residue_line.at(13) == '!' );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  97) 		// If this line indicates a problem residue, then increment the residue counter and skip it
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  98) //		if (dssp_entry_is_null) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100  99) ////			cerr << "Chain ? : " << residue::null_residue << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 100) //			new_residues.push_back(residue::NULL_RESIDUE);
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 101) //			++residue_ctr;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 102) //			continue;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 103) //		}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 104) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 105) //		cerr << "Parse residue " << residue_ctr << " from \"" << dssp_residue_line << "\"" << endl;
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 106) 		const size_residue_pair  residue_details      = parse_dssp_residue_line( dssp_residue_line );
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 107) 		const size_t            &parsed_residue_index = residue_details.first;
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 108) 		const residue           &parsed_residue       = residue_details.second;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 109) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 110) 		if ( ! dssp_entry_is_null && parsed_residue_index != residue_ctr ) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 111) 			BOOST_THROW_EXCEPTION(runtime_error_exception("Error in DSSP sequential residue numbers"));
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 112) 		}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 113) //		cerr << "Chain " << the_chain_label << " : " << parsed_residue << endl;
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 114) 
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 115) 		// Some PDBs (eg 4tsw) may have erroneous consecutive duplicate residues.
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 116) 		// Though that's a bit rubbish, it shouldn't break the whole comparison
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 117) 		// so if that's detected, just warn and move on (without appending to new_residues).
614fe2c6 (Tony Lewis 2016-06-27 16:43:42 +0100 118) 		const bool this_is_not_null = ! dssp_entry_is_null;
614fe2c6 (Tony Lewis 2016-06-27 16:43:42 +0100 119) 		const bool has_prev         = ! new_residues.empty();
614fe2c6 (Tony Lewis 2016-06-27 16:43:42 +0100 120) 		const bool prev_is_not_null = ( has_prev && ! is_null_residue( new_residues.back() ) );
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 121) 		if ( this_is_not_null && prev_is_not_null && new_residues.back().get_pdb_residue_id() == parsed_residue.get_pdb_residue_id() ) {
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 122) 			BOOST_LOG_TRIVIAL( warning ) << "Whilst parsing DSSP file, found conflicting consecutive entries for residue \""
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 123) 				<< parsed_residue.get_pdb_residue_id()
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 124) 				<< "\" (with amino acids \""
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 125) 				<< new_residues.back().get_amino_acid().get_code()
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 126) 				<< "\" and then \""
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 127) 				<< parsed_residue.get_amino_acid().get_code()
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 128) 				<< "\") - ignoring latter entry";
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 129) 		}
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 130) 		else {
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 131) 			new_residues.push_back( parsed_residue );
86abc80d (Tony Lewis 2016-05-26 18:43:27 +0100 132) 		}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 133) 		++residue_ctr;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 134) 	}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 135) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 136) 	return dssp_file( new_residues );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 137) }
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 138) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 139) /// \brief TODOCUMENT
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 140) size_residue_pair cath::file::parse_dssp_residue_line(const string &arg_dssp_residue_line ///< The DSSP residue line to parse
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 141)                                                       ) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 142) 	try {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 143) 		const bool dssp_entry_is_null = ( arg_dssp_residue_line.at( 13 ) == '!' );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 144) 		if ( dssp_entry_is_null ) {
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 145) 			return { 0_z, residue::NULL_RESIDUE };
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 146) 		}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 147) 		                                                                                         // Comments with DSSP format documentation
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 148) 		                                                                                         // (http://swift.cmbi.ru.nl/gv/dssp/)
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 149) 		const size_t        seq_res_num    =     stoul( arg_dssp_residue_line.substr(  0, 5)) ;  //   1 -   5    sequential resnumber, including chain breaks as extra residues
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 150) 		const string        res_name_str   = trim_copy( arg_dssp_residue_line.substr(  5, 6)) ;  //   6 -  11    original PDB resname, not nec. sequential, may contain letters
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 151) 		const char         &chain_char     =            arg_dssp_residue_line.at    ( 11   )  ;  //  12
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 152) 		const char         &amino_acid_c   =            arg_dssp_residue_line.at    ( 13   )  ;  //  14          amino acid sequence in one letter code
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 153) 		const char         &sstruc_code    =            arg_dssp_residue_line.at    ( 16   )  ;  //  17          secondary structure summary based on columns 19-38
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 154) //		const char         &turn3_helix    =            arg_dssp_residue_line.at    ( 18   )  ;  //  19          3-turns/helix
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 155) //		const char         &turn4_helix    =            arg_dssp_residue_line.at    ( 19   )  ;  //  20          4-turns/helix
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 156) //		const char         &turn5_helix    =            arg_dssp_residue_line.at    ( 20   )  ;  //  21          5-turns/helix
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 157) //		const char         &geom_bend      =            arg_dssp_residue_line.at    ( 21   )  ;  //  22          geometrical bend
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 158) //		const char         &chirality      =            arg_dssp_residue_line.at    ( 22   )  ;  //  23          chirality
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 159) //		const char         &beta_brid_1    =            arg_dssp_residue_line.at    ( 23   )  ;  //  24          beta bridge label
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 160) //		const char         &beta_brid_2    =            arg_dssp_residue_line.at    ( 24   )  ;  //  25          beta bridge label
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 161) //		const size_t        beta_part_1    =     stoul( arg_dssp_residue_line.substr( 25, 4)) ;  //  26 -  29    beta bridge partner resnum
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 162) //		const size_t        beta_part_2    =     stoul( arg_dssp_residue_line.substr( 29, 4)) ;  //  30 -  33    beta bridge partner resnum
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 163) //		const char         &sheet_label    =            arg_dssp_residue_line.at    ( 33   )  ;  //  34          beta sheet label
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 164) 		const size_t        solv_access    =     stoul( arg_dssp_residue_line.substr( 34, 4)) ;  //  35 -  38    solvent accessibility
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 165) 		// There are other four other columns here (N-H-->O, O-->H-N, N-H-->O and O-->H-N)
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 166) //		const double        tco            =      stod( arg_dssp_residue_line.substr( 85, 6)) ;  //  86 -  91    cosine of angle between C=O of residue I and C=O of residue I-1. For α-helices, TCO is near +1, for β-sheets TCO is near -1. Not used for structure definition.
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 167) //		const double        kappa          =      stod( arg_dssp_residue_line.substr( 91, 6)) ;  //  92 -  97    virtual bond angle (bend angle) defined by the three Cα atoms of residues I-2,I,I+2. Used to define bend (structure code 'S').
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 168) //		const double        alpha          =      stod( arg_dssp_residue_line.substr( 97, 6)) ;  //  98 - 103    virtual torsion angle (dihedral angle) defined by the four Cα atoms of residues I-1,I,I+1,I+2.Used to define chirality (structure code '+' or '-').
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 169) 		const double        phi_in_degrees =      stod( arg_dssp_residue_line.substr(103, 6)) ;  // 104 - 109    IUPAC peptide backbone torsion angles
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 170) 		const double        psi_in_degrees =      stod( arg_dssp_residue_line.substr(109, 6)) ;  // 110 - 115    IUPAC peptide backbone torsion angles
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 171) 		const double        carbon_a_x     =      stod( arg_dssp_residue_line.substr(115, 7)) ;  // 116 - 122    echo of Ca atom coordinates
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 172) 		const double        carbon_a_y     =      stod( arg_dssp_residue_line.substr(122, 7)) ;  // 123 - 129    echo of Ca atom coordinates
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 173) 		const double        carbon_a_z     =      stod( arg_dssp_residue_line.substr(129, 7)) ;  // 130 - 136    echo of Ca atom coordinates
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 174) 
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 175) 		// Parse the residue name string into a residue_id
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 176) 		const residue_id res_id{ chain_label( chain_char ), make_residue_name( res_name_str ) };
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 177) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 178) 		// Convert the DSSP secondary structure summary code into a sec_struc_type
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 179) 		const sec_struc_type sec_struc = sstruc_code == 'H' ? sec_struc_type::ALPHA_HELIX
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 180) 		                               : sstruc_code == 'E' ? sec_struc_type::BETA_STRAND
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 181) 		                                                    : sec_struc_type::COIL;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 182) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 183) //		cerr << "seq_res_num    : " << seq_res_num    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 184) //		cerr << "res_name       : " << res_name       << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 185) //		cerr << "chain_char     : " << chain_char     << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 186) //		cerr << "amino_acid     : " << amino_acid     << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 187) //		cerr << "sec_struc      : " << sec_struc      << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 188) //		cerr << "turn3_helix    : " << turn3_helix    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 189) //		cerr << "turn4_helix    : " << turn4_helix    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 190) //		cerr << "turn5_helix    : " << turn5_helix    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 191) //		cerr << "geom_bend      : " << geom_bend      << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 192) //		cerr << "chirality      : " << chirality      << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 193) //		cerr << "beta_brid_1    : " << beta_brid_1    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 194) //		cerr << "beta_brid_2    : " << beta_brid_2    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 195) //		cerr << "beta_part_1    : " << beta_part_1    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 196) //		cerr << "beta_part_2    : " << beta_part_2    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 197) //		cerr << "sheet_label    : " << sheet_label    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 198) //		cerr << "solv_access    : " << solv_access    << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 199) //		cerr << "tco            : " << tco            << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 200) //		cerr << "kappa          : " << kappa          << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 201) //		cerr << "alpha          : " << alpha          << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 202) //		cerr << "phi_in_degrees : " << phi_in_degrees << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 203) //		cerr << "psi_in_degrees : " << psi_in_degrees << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 204) //		cerr << "carbon_a_x     : " << carbon_a_x     << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 205) //		cerr << "carbon_a_y     : " << carbon_a_y     << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 206) //		cerr << "carbon_a_z     : " << carbon_a_z     << endl;
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 207) 
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 208) 		const auto shifted_phi = make_angle_from_degrees<double>( round( phi_in_degrees > 0 ? phi_in_degrees : 360.0 + phi_in_degrees) );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 209) 		const auto shifted_psi = make_angle_from_degrees<double>( round( psi_in_degrees > 0 ? psi_in_degrees : 360.0 + psi_in_degrees) );
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 210) 
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 211) 		return {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 212) 			seq_res_num,
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 213) 			
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 214) 			residue(
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 215) 				res_id,
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 216) 				amino_acid(amino_acid_c),
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 217) 				coord( carbon_a_x, carbon_a_y, carbon_a_z ),
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 218) 				coord::ORIGIN_COORD,
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 219) 				0,
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 220) 				sec_struc,
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 221) 				rotation::IDENTITY_ROTATION(),
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 222) 				shifted_phi,
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 223) 				shifted_psi,
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 224) 				solv_access
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 225) 			)
dd8438d9 (Tony Lewis 2016-12-07 13:37:35 +0000 226) 		};
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 227) 	}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 228) 	catch ( const boost::bad_lexical_cast & ) {
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 229) 		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to cast a column whilst parsing a DSSP residue record.\nRecord was \"" + arg_dssp_residue_line + "\""));
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 230) 	}
ffec4263 (Tony Lewis 2015-07-21 11:25:24 +0100 231) }
