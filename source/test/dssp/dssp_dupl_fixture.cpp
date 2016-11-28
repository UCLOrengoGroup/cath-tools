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
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/irange.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "common/algorithm/copy_build.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/debug_numeric_cast.h"
#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "structure/sec_struc_calc/bifur_hbond_list.h"
#include "test/global_test_constants.h"

#include <fstream>

using namespace cath;
using namespace cath::common;
using namespace cath::sec;
using namespace std::literals::string_literals;

using boost::adaptors::filtered;
using boost::algorithm::is_space;
using boost::algorithm::starts_with;
using boost::algorithm::trim_copy;
using boost::filesystem::path;
using boost::format;
using boost::icontains;
using boost::irange;
using boost::is_any_of;
using boost::make_optional;
using boost::none;
using boost::test_tools::check_is_close;
using boost::test_tools::percent_tolerance;
using boost::token_compress_on;
using std::get;
using std::getline;
using std::ifstream;
using std::istream;
using std::min;
using std::pair;
using std::string;

/// \brief Map the specified dsspfile_hbond_opt parsed from a DSSP file with specified supplementary data
///        to a hbond_half_opt
///
/// \relates dsspfile_hbond
hbond_half_opt cath::sec::mapped_dsspfile_hbond(const dsspfile_hbond_opt &arg_dsspfile_hbond_opt,        ///< The dsspfile_hbond_opt as parsed from a DSSP file
                                                const size_t             &arg_residue_index,             ///< The index of the residue in the DSSP file
                                                const size_vec           &arg_normal_index_of_dssp_index ///< A mapping from normal index to the DSSP index
                                                ) {
	// No hbond in -> no hbond out
	if ( ! arg_dsspfile_hbond_opt ) {
		return none;
	}

	// Calculate the DSSP index of arg_dsspfile_hbond_opt and check it's non-negative
	const int dssp_index = static_cast<int>( arg_residue_index ) + arg_dsspfile_hbond_opt->first - 1;
	if ( dssp_index < 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("DSSP bond destination residue index is less than 0 (" + ::std::to_string( dssp_index ) + ")"));
	}

	// Check the DSSP index is within range
	const size_t dssp_index_size = static_cast<size_t>( dssp_index );
	if ( dssp_index_size >= arg_normal_index_of_dssp_index.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"DSSP hbond destination index is "
			+ ::std::to_string( dssp_index_size )
			+ ", which is greater than or equal to the number of entries ("
			+ ::std::to_string( arg_normal_index_of_dssp_index.size() )
			+ ")"
		));
	}

	// Generate the hbond_half
	return { hbond_half{
		debug_numeric_cast< hbond_partner_t >( arg_normal_index_of_dssp_index[ dssp_index_size ] ),
		arg_dsspfile_hbond_opt->second
	} };
}

/// \brief Generate a string describing the difference between the specified DSSP parsed residue and calculated residue
///        (in the context described by the specified string) or return none if there isn't any difference
///
/// \relates hbond_half
str_opt cath::sec::difference_string(const string         &arg_context_str,        ///< A string describing the context of this comparison
                                     const hbond_half_opt &arg_dsspfile_hbond_opt, ///< The hbond_half_opt parsed from the DSSP file (and converted via mapped_dsspfile_hbond())
                                     const hbond_half_opt &arg_hbond_half_opt      ///< The hbond_half_opt calculated independently
                                     ) {
	if ( ( ! arg_dsspfile_hbond_opt ) != ( ! arg_hbond_half_opt ) ) {
		if ( arg_dsspfile_hbond_opt ) {
			return "DSSP has a " + arg_context_str + " bond where one hasn't been calculated"s;
		}
		else {
			return "Calculated a " + arg_context_str + " bond (" + to_string( *arg_hbond_half_opt )  + ") where DSSP doesn't have one"s;
		}
	}
	if ( arg_dsspfile_hbond_opt && arg_hbond_half_opt ) {
		if ( arg_dsspfile_hbond_opt->index != arg_hbond_half_opt->index ) {
			return "DSSP bond destination residue index of "
				+ ::std::to_string( arg_dsspfile_hbond_opt->index )
				+ " doesn't match calculated index of "
				+ ::std::to_string( arg_hbond_half_opt->index );
		}

		const auto rounded_calc_energy = stod( ( format("%3.1f") % arg_hbond_half_opt->energy ).str() );
		if ( ! check_is_close( arg_dsspfile_hbond_opt->energy, rounded_calc_energy, boost::test_tools::percent_tolerance( 0.0001 ) ) ) {
			return "DSSP calculated bond energy of "
				+ ::std::to_string( arg_dsspfile_hbond_opt->energy )
				+ " doesn't match calculated energy of "
				+ ::std::to_string( rounded_calc_energy )
				+ " (rounded from "
				+ ::std::to_string( arg_hbond_half_opt->energy )
				+ ")";
		}
	}
	return none;
}

/// \brief Generate a string describing the difference between the specified dssp_dupl_res parsed from a DSSP file
///        and the specified bifur_hbond or return none if there isn't any difference
///
/// \relates dssp_dupl_res
str_opt cath::sec::difference_string(const dssp_dupl_res &arg_dssp_dupl_res,             ///< A dssp_dupl_res parsed from a DSSP file
                                     const bifur_hbond   &arg_bifur_hbond,               ///< The dssp_dupl_res to be compared to the dssp_dupl_res
                                     const size_vec      &arg_normal_index_of_dssp_index ///< A mapping from normal index to the DSSP index 
                                     ) {
	const auto &res_idx = arg_dssp_dupl_res.residue_index;
	const str_opt nh_1st = difference_string( "nh-1st", mapped_dsspfile_hbond( arg_dssp_dupl_res.hbonds_this_nh_1st_2nd.first,  res_idx, arg_normal_index_of_dssp_index ), arg_bifur_hbond.get_bound_pair_for_this_nh().first  );
	const str_opt nh_2nd = difference_string( "nh-2nd", mapped_dsspfile_hbond( arg_dssp_dupl_res.hbonds_this_nh_1st_2nd.second, res_idx, arg_normal_index_of_dssp_index ), arg_bifur_hbond.get_bound_pair_for_this_nh().second );
	const str_opt co_1st = difference_string( "co-1st", mapped_dsspfile_hbond( arg_dssp_dupl_res.hbonds_this_co_1st_2nd.first,  res_idx, arg_normal_index_of_dssp_index ), arg_bifur_hbond.get_bound_pair_for_this_co().first  );
	const str_opt co_2nd = difference_string( "co-2nd", mapped_dsspfile_hbond( arg_dssp_dupl_res.hbonds_this_co_1st_2nd.second, res_idx, arg_normal_index_of_dssp_index ), arg_bifur_hbond.get_bound_pair_for_this_co().second );

	const string suffix_string = " at DSSP residue "
		+ to_string( arg_dssp_dupl_res.pdb_residue_name )
		+ " (index "
		+ ::std::to_string( res_idx )
		+ ")";
	return
		nh_1st ? make_optional( *nh_1st + suffix_string ):
		nh_2nd ? make_optional( *nh_2nd + suffix_string ):
		co_1st ? make_optional( *co_1st + suffix_string ):
		co_2nd ? make_optional( *co_2nd + suffix_string ):
		         none;
}

/// \brief Generate a string describing the difference between the specified dssp_dupl_ress parsed from a DSSP file
///        and the specified bifur_hbond_list or return none if there isn't any difference
///
/// \relates dssp_dupl_res
str_opt cath::sec::difference_string(const dssp_dupl_res_vec &arg_dssp_dupl_res_vec, ///< A vector of dssp_dupl_res parsed from a DSSP file
                                     const bifur_hbond_list  &arg_bifur_hbond_list   ///< The dssp_dupl_res_list to be compared to the dssp_dupl_res vector
                                     ) {
	const auto num_dssp_dupl_res = arg_dssp_dupl_res_vec.size();
	const auto arg_bifur_hbonds  = arg_bifur_hbond_list.size();

	const auto dssp_non_null_indices = copy_build<size_vec>(
		irange( 0_z, num_dssp_dupl_res )
			| filtered(
				[&] (const size_t &x) {
					return ! arg_dssp_dupl_res_vec[ x ].pdb_residue_name.get_is_null_residue_name();
				}
			)
	);

	const size_t num_non_null_residues = dssp_non_null_indices.size();

	const str_opt num_prob = make_optional(
		( num_non_null_residues != arg_bifur_hbonds ),
		"Number of (non-null) DSSP hbond residues ("
			+ ::std::to_string( num_non_null_residues )
			+ "), doesn't match the number of calculated hbond residues ("
			+ ::std::to_string( arg_bifur_hbonds )
			+ "). "
	);

	const auto normal_index_of_dssp_index = [&] () {
		size_vec result( num_dssp_dupl_res, 0 );
		for (const size_t &x : irange( 0_z, dssp_non_null_indices.size() ) ) {
			result[ dssp_non_null_indices[ x ] ] = x;
		}
		return result;
	}();

	for (const size_t &index : irange( 0_z, min( num_dssp_dupl_res, arg_bifur_hbonds ) ) ) {
		const auto diff_str = difference_string( arg_dssp_dupl_res_vec[ dssp_non_null_indices[ index ] ], arg_bifur_hbond_list[ index ], normal_index_of_dssp_index );
		if ( diff_str ) {
			return ( num_prob ? *num_prob : string{} ) + *diff_str;
		}
	}

	return num_prob
		? make_optional( *num_prob + " No differences spotted until the end of the shortest" )
		: none;
}

/// \brief Parse a DSSP file from the specified file into a dssp_dupl_res_vec for testing
dssp_dupl_res_vec dssp_dupl_fixture::parse_dssp_for_calc_testing(const path &arg_dssp_file ///< The file from which to parse the data
                                                                 ) {
	ifstream my_dssp_istream;
	open_ifstream( my_dssp_istream, arg_dssp_file );
	const auto parsed_results = parse_dssp_for_calc_testing( my_dssp_istream );
	my_dssp_istream.close();
	return parsed_results;
}

/// \brief Parse a DSSP file from the specified istream into a dssp_dupl_res_vec for testing
dssp_dupl_res_vec dssp_dupl_fixture::parse_dssp_for_calc_testing(istream &arg_dssp_stream ///< The istream from which to parse the data
                                                                 ) {
	// Read the first line and check it looks vaguely like a DSSP header line
	string dssp_first_line;
	getline(arg_dssp_stream, dssp_first_line);
	if ( !starts_with(dssp_first_line, "==== ") || !icontains(dssp_first_line, "DSSP") || dssp_first_line.length() > 300 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("First line is not a DSSP header line"));
	}


	// Scan through the header
	bool found_column_headings( false );
	while ( ! found_column_headings && ! arg_dssp_stream.eof() ) {
		string dssp_header_line;
		getline( arg_dssp_stream, dssp_header_line );
		found_column_headings = starts_with(dssp_header_line, "  #  RESIDUE AA ");
	}

	// Read residues
	size_t residue_ctr = 1;
	dssp_dupl_res_vec dssp_dupl_residues;
	while ( ! arg_dssp_stream.eof() ) {
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

/// \brief Parse data from the specified DSSP line
pair<size_t, dssp_dupl_res> dssp_dupl_fixture::parse_dssp_residue_line(const string &arg_dssp_residue_line ///< The DSSP residue line to parse
                                                                       ) {
	// try {
		const bool dssp_entry_is_null = ( arg_dssp_residue_line.at( 13 ) == '!' );
		if ( dssp_entry_is_null ) {
			return { 0_z, dssp_dupl_res{} };
		}
		const size_t seq_res_num  =     stoul( arg_dssp_residue_line.substr(  0, 5)) ;  //   1 -   5    sequential resnumber, including chain breaks as extra residues
		const string res_name_str = trim_copy( arg_dssp_residue_line.substr(  5, 6)) ;  //   6 -  11    original PDB resname, not nec. sequential, may contain letters

		const string hbond_a_str  = trim_copy( arg_dssp_residue_line.substr( 39, 11 ) );
		const string hbond_b_str  = trim_copy( arg_dssp_residue_line.substr( 50, 11 ) );
		const string hbond_c_str  = trim_copy( arg_dssp_residue_line.substr( 61, 11 ) );
		const string hbond_d_str  = trim_copy( arg_dssp_residue_line.substr( 72, 11 ) );

		const auto   hbond_a      = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 39, 11 ) ) );
		const auto   hbond_b      = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 50, 11 ) ) );
		const auto   hbond_c      = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 61, 11 ) ) );
		const auto   hbond_d      = parse_dsspfile_bond( trim_copy( arg_dssp_residue_line.substr( 72, 11 ) ) );

		const residue_name res_name = make_residue_name( res_name_str );

		return {
			seq_res_num,
			dssp_dupl_res{
				make_pair( hbond_a, hbond_c ),
				make_pair( hbond_b, hbond_d ),
				seq_res_num,
				res_name
			},
		};
	// }
	// catch ( const boost::bad_lexical_cast & ) {
	// 	BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to cast a column whilst parsing a DSSP residue record.\nRecord was \"" + arg_dssp_residue_line + "\""));
	// }
}

/// \brief Parse an hbond from the specified string from within a DSSP line
dsspfile_hbond_opt dssp_dupl_fixture::parse_dsspfile_bond(const string &arg_hbond_string ///< The string containing the DSSP h-bond data (eg "-2,-2.6")
                                                          ) {
	const str_vec parts = split_build<str_vec>( arg_hbond_string, is_any_of( "," ), token_compress_on );
	if ( parts.size() != 2 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Did not find two parts in DSSP file h-bond"));
	}
	const int    offset = stoi( trim_copy( parts.front() ) );
	const double energy = stod( trim_copy( parts.back () ) );
	if ( offset == 0 ) {
		if ( energy != 0.0 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Whilst try to parse H-bond from DSSP file data, non-zero energy for zero offset"));
		}
		return none;
	}
	else {
		return dsspfile_hbond{ offset, energy };
	}
}

/// \brief Test constant for the DSSP hbond test data subdirectory
const path & dssp_dupl_fixture::DSSP_HBOND_TEST_DATA_DIR() {
	static const path dssp_hbond_test_data_dir( global_test_constants::TEST_SOURCE_DATA_DIR() / "dssp" / "hbond" );
	return dssp_hbond_test_data_dir;
}
