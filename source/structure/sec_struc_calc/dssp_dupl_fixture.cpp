/// \file
/// \brief The dssp_dupl_fixture class definitions

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

#include "dssp_dupl_fixture.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp> // ** USING THIS? **
#include <boost/log/trivial.hpp>

#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"

#include <fstream>
#include <iostream> // ****** TEMPORARY ******

using namespace cath;
using namespace cath::common;
using namespace cath::sec;

using boost::algorithm::is_space;
using boost::algorithm::starts_with;
using boost::algorithm::trim_copy;
using boost::filesystem::path;
using boost::icontains;
using boost::is_any_of;
using boost::none;
using boost::token_compress_on;
using std::get;
using std::getline;
using std::ifstream;
using std::istream;
using std::pair;
using std::string;

/// \brief TODOCUMENT
dssp_dupl_res_vec dssp_dupl_fixture::parse_dssp_for_calc_testing(const path &arg_dssp_file ///< TODOCUMENT
                                                                 ) {
	ifstream my_dssp_istream;
	open_ifstream( my_dssp_istream, arg_dssp_file );
	const auto parsed_results = parse_dssp_for_calc_testing( my_dssp_istream );
	my_dssp_istream.close();
	return parsed_results;
}

dssp_dupl_res_vec dssp_dupl_fixture::parse_dssp_for_calc_testing(istream &arg_dssp_stream ///< TODOCUMENT
                                                                 ) {
	// Read the first line and check it looks vaguely like a DSSP header line
	string dssp_first_line;
	getline(arg_dssp_stream, dssp_first_line);
	if ( !starts_with(dssp_first_line, "==== ") || !icontains(dssp_first_line, "DSSP") || dssp_first_line.length() > 300 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("First line is not a DSSP header line"));
	}


	// Scan through the header
	bool found_column_headings( false );
	while ( ! found_column_headings ) {
		string dssp_header_line;
		getline( arg_dssp_stream, dssp_header_line );
		found_column_headings = starts_with(dssp_header_line, "  #  RESIDUE AA ");
	}

	// Read residues
	size_t residue_ctr = 1;
	dssp_dupl_res_vec dssp_dupl_residues;
	while ( ! arg_dssp_stream.eof()) {
		string dssp_residue_line;
		getline( arg_dssp_stream, dssp_residue_line );

		// If this line contains nothing but whitespace, then skip it
		if ( all( dssp_residue_line, is_space() ) ) {
			continue;
		}

		const bool dssp_entry_is_null = ( dssp_residue_line.at(13) == '!' );

		const auto           residue_details      = parse_dssp_residue_line( dssp_residue_line );
		const size_t        &parsed_residue_index = get<0>( residue_details );
		const dssp_dupl_res &dssp_dupl_residue    = get<1>( residue_details );

		if ( ! dssp_entry_is_null && parsed_residue_index != residue_ctr ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				"Error in DSSP sequential residue numbers "
				+ std::to_string( parsed_residue_index )
				+ " "
				+ std::to_string( residue_ctr          )
			));
		}

		// Some PDBs (eg 4tsw) may have erroneous consecutive duplicate residues.
		// Though that's a bit rubbish, it shouldn't break the whole comparison
		// so if that's detected, just warn and move on (without appending to new_residues).
		const bool this_is_not_null = ! dssp_entry_is_null;
		const bool has_prev         = ! dssp_dupl_residues.empty();
		const bool prev_is_not_null = ( has_prev && dssp_dupl_residues.back().residue_index > 0 );
		if ( this_is_not_null && prev_is_not_null && dssp_dupl_residues.back().pdb_residue_name== dssp_dupl_residue.pdb_residue_name ) {
			BOOST_LOG_TRIVIAL( warning ) << "Whilst parsing DSSP file, found conflicting consecutive entries for residue \""
				<< dssp_dupl_residue.pdb_residue_name
				<< "\" - ignoring latter entry";
		}
		else {
			dssp_dupl_residues.push_back( dssp_dupl_residue );
		}
		++residue_ctr;
	}

	return dssp_dupl_residues;
}

/// \brief TODOCUMENT
pair<size_t, dssp_dupl_res> dssp_dupl_fixture::parse_dssp_residue_line(const string &arg_dssp_residue_line ///< The DSSP residue line to parse
                                                                       ) {
	try {
		const bool dssp_entry_is_null = ( arg_dssp_residue_line.at( 13 ) == '!' );
		if ( dssp_entry_is_null ) {
			return { 0_z, dssp_dupl_res{} };
		}
		const size_t seq_res_num  =     stoul( arg_dssp_residue_line.substr(  0, 5)) ;  //   1 -   5    sequential resnumber, including chain breaks as extra residues
		const string res_name_str = trim_copy( arg_dssp_residue_line.substr(  5, 6)) ;  //   6 -  11    original PDB resname, not nec. sequential, may contain letters

		const string h_bond_a_str = trim_copy( arg_dssp_residue_line.substr( 39, 11 ) );
		const string h_bond_b_str = trim_copy( arg_dssp_residue_line.substr( 50, 11 ) );
		const string h_bond_c_str = trim_copy( arg_dssp_residue_line.substr( 61, 11 ) );
		const string h_bond_d_str = trim_copy( arg_dssp_residue_line.substr( 72, 11 ) );

		const auto   h_bond_a     = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 39, 11 ) ) );
		const auto   h_bond_b     = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 50, 11 ) ) );
		const auto   h_bond_c     = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 61, 11 ) ) );
		const auto   h_bond_d     = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 72, 11 ) ) );

		const residue_name res_name = make_residue_name( res_name_str );

		std::cerr << "res_name is     : \"" << res_name     << "\"\n";
		std::cerr << "h_bond_a_str is : " << ( h_bond_a ? h_bond_a_str : "no-bond" ) << "\n";
		std::cerr << "h_bond_b_str is : " << ( h_bond_b ? h_bond_b_str : "no-bond" ) << "\n";
		std::cerr << "h_bond_c_str is : " << ( h_bond_c ? h_bond_c_str : "no-bond" ) << "\n";
		std::cerr << "h_bond_d_str is : " << ( h_bond_d ? h_bond_d_str : "no-bond" ) << "\n";

		return { seq_res_num, dssp_dupl_res{} };
	}
	catch ( const boost::bad_lexical_cast & ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to cast a column whilst parsing a DSSP residue record.\nRecord was \"" + arg_dssp_residue_line + "\""));
	}
}

/// \brief TODOCUMENT
dsspfile_h_bond_opt dssp_dupl_fixture::parse_dsspfile_bond(const string &arg_h_bond_string ///< TODOCUMENT
                                                           ) {
	const str_vec parts = split_build<str_vec>( arg_h_bond_string, is_any_of( "," ), token_compress_on );
	if ( parts.size() != 2 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Did not find two parts in DSSP file h-bond"));
	}
	const int    offset = stoi( trim_copy( parts.front() ) );
	const double energy = stod( trim_copy( parts.front() ) );
	if ( offset == 0 ) {
		if ( energy != 0.0 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Whilst try to parse H-bond from DSSP file data, non-zero energy for zero offset"));
		}
		return none;
	}
	else {
		return dsspfile_h_bond{ offset, energy };
	}
}