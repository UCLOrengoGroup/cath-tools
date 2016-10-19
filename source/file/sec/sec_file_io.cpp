/// \file
/// \brief The sec_file_io definitions

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

#include "sec_file_io.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/file/open_fstream.h"
#include "common/lexical_cast_line.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "exception/runtime_error_exception.h"
#include "file/sec/sec_file.h"
#include "file/sec/sec_file_record.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

#include <fstream>
#include <iostream>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::file::detail;
using namespace cath::geom;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::token_compress_on;
using boost::algorithm::trim_copy;
using boost::lexical_cast;

/// \brief Read a sec file into a sec_file object
///
/// \relates sec_file
sec_file cath::file::read_sec(const path &arg_sec_filename ///< The file from which to parse the sec data
                              ) {
	ifstream my_sec_istream;
	open_ifstream(my_sec_istream, arg_sec_filename);
	const sec_file the_sec_file = read_sec( my_sec_istream );
	my_sec_istream.close();
	return the_sec_file;
}

/// \brief Read a sec file from an istream into a sec_file object
///
/// \relates sec_file
///
/// An example sec file (1cukA02.sec, circa CATH v4.0.0) :
///
///     9
///          1    S *     4     7  -11.40   1.59  -8.55    -0.19 -0.56 -0.80
///          71.69 -179.41 -22.77     121.00 -156.16 -78.05     -74.07 27.49 114.77     138.09 -141.65 -103.06     -32.10 -13.79 -29.50     174.84 175.75 176.85     128.14 130.89 48.75     -44.59 -72.60 -114.25
///          2    S *     7    10  -10.92   9.91  -5.50     0.07 -0.96  0.29
///          49.31 23.25 -55.29     -145.76 -153.11 137.54     66.39 37.75 -80.29     -103.79 165.61 -6.73     103.15 -4.84 -160.38     56.45 -49.70 71.52     -116.28 106.81 -91.48
///          3    S *    10    13  -12.57  16.34 -10.48     0.57 -0.33  0.75
///          164.93 -176.36 -167.18     17.08 14.50 -25.01     -153.10 142.36 48.56     53.84 -28.09 -105.10      7.14 -72.95 126.80     -165.59 83.56 -36.19
///          4    S *    16    21   -9.69  10.38 -10.75    -0.56  0.52 -0.65
///          -147.85 -169.14 142.17     41.97 -41.28 -144.27     -111.09 148.27 62.08     -157.79 103.41 -66.02     29.48 -100.08 130.98
///          5    S *    24    29   -7.75   8.94 -14.45     0.78 -0.07  0.62
///          -170.18 127.86 73.56     36.75 -42.59 -80.09     -9.94 -87.46 151.81     177.33 69.06 -11.18
///          6    H *    32    35  -23.59  10.63 -15.68     0.01 -0.05 -1.00
///          -153.06 -170.45 -153.65     160.24 144.68 78.25     -12.49 -58.80 -84.75
///          7    S *    44    48  -15.86   3.62  -8.05     0.14  0.49  0.86
///          -46.70 -44.86 -128.10     140.57 111.65 68.90
///          8    S *    48    54  -11.34   1.01 -22.49    -0.57 -0.24  0.79
///          -172.73 156.51 -162.99
///          9    S *    57    63  -11.17   5.11 -20.75     0.86  0.09 -0.51
sec_file cath::file::read_sec(istream &arg_istream ///< The istream from which to parse the sec file
                              ) {
	// Grab the number of secondary structures from the first line
	size_t num_sec_strucs;
	try {
		num_sec_strucs = lexical_cast_trimmed_line<size_t>(arg_istream);
	}
	catch (const boost::bad_lexical_cast &) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Invalid sec file - first line does not consist of a valid number of secondary structures"));
	}

	// Prepare vectors of sec_file_record and sec_struc_planar_angles_vec to populate
	vector<sec_file_record> new_sec_records;
	new_sec_records.reserve(num_sec_strucs);

	sec_struc_planar_angles_vec_vec new_planar_angles;
	new_planar_angles.reserve(num_sec_strucs);

	// Read data for secondary structures
	for (size_t sec_record_ctr = 0; sec_record_ctr < num_sec_strucs; ++sec_record_ctr) {
		// Parse the next line as a sec main line
		string sec_main_line;
		getline(arg_istream, sec_main_line);
		pair<size_t, sec_file_record> index_and_sec_record = parse_sec_main_line( sec_main_line );
		const size_t    &new_sec_record_index = index_and_sec_record.first;
		sec_file_record &new_sec_record       = index_and_sec_record.second;

		// Check that the secondary structure's number is as expected
		if (new_sec_record_index != sec_record_ctr + 1) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				"Invalid sec file : entry number " + lexical_cast<string>(sec_record_ctr + 1) +
				" is labelled as secondary structure number " + lexical_cast<string>(new_sec_record_index)
			));
		}

		// Add new_sec_record to the back of new_sec_records
		new_sec_records.push_back(new_sec_record);

		// Parse the next line as a sec angles line, check the number of entries
		// and then add them to the back new_planar_angles
		string sec_angles_line;
		getline(arg_istream, sec_angles_line);
		const sec_struc_planar_angles_vec parsed_planar_angles = parse_sec_angles_line(sec_angles_line);
		if (num_sec_strucs - sec_record_ctr - 1 != parsed_planar_angles.size()) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("The number of planar angles parsed from sec file was not as expected"));
		}
		if (!parsed_planar_angles.empty()) {
			new_planar_angles.push_back(parsed_planar_angles);
		}
	}

	// If the input stream is not at the end, try gobbling up all whitespace and then throw an error if it's still not at the end
	if (!arg_istream.eof()) {
		ws(arg_istream);
		if (!arg_istream.eof()) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Further data remains in the input stream after successfully parsing a sec file from it"));
		}
	}

	return sec_file(new_sec_records, new_planar_angles);
}

/// \brief Parse a main line from a sec file out of the specified string
///
/// \relates sec_struc
///
/// \todo Consider moving part of this code into a sec_file_record extraction operator
pair<size_t, sec_file_record> cath::file::detail::parse_sec_main_line(const string &arg_sec_line_string ///< String containing the sec file main line
                                                                      ) {
	// Prepare variables to hold the parsed data
	size_t         sec_struc_number;
	sec_struc_type helix_or_strand;
	char           unused_strand_type;
	size_t         start_residue_num, stop_residue_num;
	double         midpoint_x,  midpoint_y,  midpoint_z;
	double         unit_dirn_x, unit_dirn_y, unit_dirn_z;

	// Put the string into an istringstream and then parse out the parts
	// (and do this in a try block so that cast conversion failures can be caught)
	istringstream sec_line_ss(arg_sec_line_string);
	try {
		sec_line_ss >> sec_struc_number;
		sec_line_ss >> helix_or_strand;
		sec_line_ss >> unused_strand_type;
		sec_line_ss >> start_residue_num >> stop_residue_num;
		sec_line_ss >> midpoint_x  >> midpoint_y  >> midpoint_z;
		sec_line_ss >> unit_dirn_x >> unit_dirn_y >> unit_dirn_z;
	}
	catch (const boost::bad_lexical_cast &) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to cast a column whilst parsing a main lie record from a sec file, which probably means the line's malformed.\nRecord was \"" + arg_sec_line_string + "\""));
	}

	// Return a pair containing the index and a newly populated sec_file_record
	return make_pair(
		sec_struc_number,
		sec_file_record(
			start_residue_num,
			stop_residue_num,
			helix_or_strand,
			coord( midpoint_x,  midpoint_y,  midpoint_z ),
			coord( unit_dirn_x, unit_dirn_y, unit_dirn_z )
		)
	);
}

/// \brief Parse sec_struc_planar angles out a sec file planar angles line
///
/// \relates sec_struc
sec_struc_planar_angles_vec cath::file::parse_sec_angles_line(const string &arg_sec_line_string ///< String containing the sec file planar angles line
                                                              ) {
	// Trim any whitespace of the ends of the string and just return an empty vector if that's empty
	const string trimmed_sec_angles_string(trim_copy(arg_sec_line_string));
	if (trimmed_sec_angles_string.empty()) {
		return sec_struc_planar_angles_vec();
	}

	// Split the line on whitespace and throw if the number of parts doesn't divide by 3
	str_deq line_parts = split_build<str_deq>( trimmed_sec_angles_string, is_any_of( " " ), token_compress_on );
	if ( line_parts.size() % 3 != 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse planar angles from line in sec file because number of parts in line is not a multiple of three"));
	}

	// Work though triplets of parts
	sec_struc_planar_angles_vec planar_angles;
	while (!line_parts.empty()) {
		assert(line_parts.size() % 3 == 0);

		// Parse out three angles
		// Do it in a try block to allow any cast conversions to be caught
		try {
			const double planar_angle_x       = stod( line_parts.front() );
			line_parts.pop_front();
			const double planar_angle_minus_y = stod( line_parts.front() );
			line_parts.pop_front();
			const double planar_angle_z       = stod( line_parts.front() );
			line_parts.pop_front();
			if ( ! isfinite( planar_angle_x ) || ! isfinite( planar_angle_minus_y ) || ! isfinite( planar_angle_z ) ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to get a finite number from a column whilst parsing a planar angles value from a sec file, which probably means the line has an infinite or 'nan' (not-a-number) value.\nRecord was \"" + arg_sec_line_string + "\""));
			}
			planar_angles.push_back(sec_struc_planar_angles(
				planar_angle_x,
				planar_angle_minus_y,
				planar_angle_z
			));
		}
		catch (const boost::bad_lexical_cast &) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to cast a column whilst parsing a planar angles value from a sec file, which probably means the line's malformed.\nRecord was \"" + arg_sec_line_string + "\""));
		}
	}

	// Return the result
	return planar_angles;
}

