/// \file
/// \brief The sec_calc definitions

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

#include "sec_calc.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/adaptor/equal_grouped.hpp"
#include "common/boost_addenda/range/back.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/size_t_literal.hpp"
#include "file/sec/sec_file.hpp"
#include "file/sec/sec_file_record.hpp"
#include "structure/geometry/line.hpp"
#include "structure/geometry/pca.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sec;
using namespace std::literals::string_literals;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::format;
using boost::integer_range;
using boost::irange;
using boost::sub_range;
using std::make_pair;
using std::string;

/// \brief Round a coord in the way they are in sec files to allow for tests to compare
///
/// \relates coord
coord cath::sec::round_like_sec_file_copy(const coord &prm_sec_coord ///< The coord to round
                                          ) {
	return {
		stod( ( format("%6.2f") % prm_sec_coord.get_x() ).str() ),
		stod( ( format("%6.2f") % prm_sec_coord.get_y() ).str() ),
		stod( ( format("%6.2f") % prm_sec_coord.get_z() ).str() )
	};
}

/// \brief Round a sec_file_record in the way they are in sec files to allow for tests to compare
///
/// \relates sec_file_record
sec_file_record cath::sec::round_like_sec_file_copy(const sec_file_record &prm_sec_file_record ///< The sec_file_record to round
                                                    ) {
	return {
		prm_sec_file_record.get_start_residue_num(),
		prm_sec_file_record.get_stop_residue_num(),
		prm_sec_file_record.get_type(),
		round_like_sec_file_copy( prm_sec_file_record.get_midpoint()  ),
		round_like_sec_file_copy( prm_sec_file_record.get_unit_dirn() )
	};
}

/// \brief Calculate the sec information for the secondary structure between the specified start and stop
///        indices on the specified protein
sec_file_record cath::sec::calculate_sec_record(const protein &prm_protein,         ///< The protein to query
                                                const size_t  &prm_start_res_index, ///< The index of the first  residue of the secondary structure
                                                const size_t  &prm_stop_res_index   ///< The index of the second residue of the secondary structure
                                                ) {
	// for (const size_t &res_ctr : irange( prm_start_res_index, prm_stop_res_index + 1_z ) ) {
	// 	std::cerr << prm_protein.get_residue_ref_of_index( res_ctr ) << "\n";
	// }

	const coord_list core_coords{ transform_build<coord_vec>(
		irange( prm_start_res_index + 1_z, prm_stop_res_index ),
		[&] (const size_t &prm_triplet_mid_res_index) {
			return prosec_axis_point_of_residue_triple( prm_protein, prm_triplet_mid_res_index );
		}
	) };

	// for (const coord &core_coord : core_coords) {
	// 	std::cerr << "Core coord : " << core_coord << "\n";
	// }
	// std::cerr << "\n";

	const auto lobf = line_of_best_fit( core_coords );
	const auto line_approach_first = closest_point_on_line_to_point( lobf, front( core_coords ) );
	const auto line_approach_last  = closest_point_on_line_to_point( lobf, back ( core_coords ) );
	const auto midpoint            = ( line_approach_first + line_approach_last ) / 2.0;
	// std::cerr << "line_approach_first : " << line_approach_first << "\n";
	// std::cerr << "line_approach_last  : " << line_approach_last  << "\n";
	// std::cerr << "midpoint            : " << midpoint            << "\n";
	// std::cerr << "Unit direction      : " << ( - normalise_copy( lobf.get_dirn() ) ) << "\n";
	return {
		prm_start_res_index + 1u,
		prm_stop_res_index  + 1u,
		prm_protein.get_residue_ref_of_index( prm_start_res_index + 1 ).get_sec_struc_type(),
		midpoint,
		( - normalise_copy( lobf.get_dirn() ) )
	};
}

/// \brief Return the prosec_sec_type corresponding to the specified sec_struc_type
///
/// \relates prosec_sec_type
prosec_sec_type cath::sec::make_prosec_sec_type(const sec_struc_type &prm_sec_struc ///< The sec_struc_type to convert
                                                ) {
	switch( prm_sec_struc ) {
		case ( sec_struc_type::ALPHA_HELIX ) : { return prosec_sec_type::ALPHA_HELIX; }
		case ( sec_struc_type::BETA_STRAND ) : { return prosec_sec_type::BETA_STRAND; }
		case ( sec_struc_type::COIL        ) : {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot make_prosec_sec_type() of sec_struc_type::COIL"));
			return prosec_sec_type::ALPHA_HELIX; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
	BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of sec_struc_type not recognised whilst making prosec_sec_type"));
}

/// \brief Get the prosec axis point of the middle of the three specified consecutive residues in the same secondary structure
///
/// The result is somewhere between the central location and the midpoint of the two straddling locations.
/// The fraction along hat line is determined by the ends_midpoint_weight_of_sec_type() for the secondary
/// structure type
coord cath::sec::prosec_axis_point_of_residue_triple(const residue &prm_residue_before, ///< The preceding residue
                                                     const residue &prm_residue,        ///< The residue
                                                     const residue &prm_residue_after   ///< The following residue
                                                     ) {
	return prosec_axis_point_of_residue_triple(
		prm_residue_before . get_carbon_alpha_coord(),
		prm_residue        . get_carbon_alpha_coord(),
		prm_residue_after  . get_carbon_alpha_coord(),
		make_prosec_sec_type( prm_residue.get_sec_struc_type() )
	);
}

/// \brief Get the prosec axis point of the middle of the three specified locations of consecutive carbon-alphas
///        that are part of the specified type of secondary structure
///
/// The result is somewhere between the central location and the midpoint of the two straddling locations.
/// The fraction along hat line is determined by the ends_midpoint_weight_of_sec_type() for the secondary
/// structure type
coord cath::sec::prosec_axis_point_of_residue_triple(const protein &prm_protein,          ///< The protein
                                                     const size_t  &prm_mid_residue_index ///< The index of the middle residue
                                                     ) {
	if ( prm_mid_residue_index == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}
	return prosec_axis_point_of_residue_triple(
		prm_protein.get_residue_ref_of_index( prm_mid_residue_index - 1_z ),
		prm_protein.get_residue_ref_of_index( prm_mid_residue_index       ),
		prm_protein.get_residue_ref_of_index( prm_mid_residue_index + 1_z )
	);
}

/// \brief Get the start/stop residue indices of the secondary structures of length >= MIN_SEC_STRUC_LENGTH
size_size_pair_vec cath::sec::get_sec_starts_and_stops(const protein &prm_protein ///< The protein to query
                                                       ) {
	constexpr size_t MIN_SEC_STRUC_LENGTH = 4_z;

	return copy_build<size_size_pair_vec>(
		indices( prm_protein.get_length() )
			// Group according to secondary structure corresponding to index in protein
			| equal_grouped( [&] (const size_t &x, const size_t &y) {
				return (
					prm_protein.get_residue_ref_of_index( x ).get_sec_struc_type()
					!=
					prm_protein.get_residue_ref_of_index( y ).get_sec_struc_type()
				);
			} )
			// Remove entries that are coil or less than MIN_SEC_STRUC_LENGTH residues long
			| filtered( [&] (const sub_range<integer_range<size_t>> &x) {
				return (
					back( x ) + 1_z - front( x ) >= MIN_SEC_STRUC_LENGTH
					&&
					prm_protein.get_residue_ref_of_index( front( x ) ).get_sec_struc_type() != sec_struc_type::COIL
				);
			} )
			// Convert to a pair of start/stop indices
			| transformed( [] (const sub_range<integer_range<size_t>> &x) {
				return make_pair( front( x ), back( x ) );
			} )
	);
}

/// \brief Generate a vector of sec_file_record objects for the specified protein
sec_file_record_vec cath::sec::get_sec_records(const protein &prm_protein ///< The protein to analyse
                                               ) {
	return transform_build<sec_file_record_vec>(
		get_sec_starts_and_stops( prm_protein ),
		[&] (const size_size_pair &x) {
			return calculate_sec_record( prm_protein, x.first, x.second );
		}
	);
}

/// \brief Generate a vector of sec_file_record objects for the specified protein
sec_file cath::sec::get_sec_file(const protein &prm_protein ///< The protein to analyse
                                 ) {
	return make_sec_file_with_calced_planar_angles( get_sec_records( prm_protein ) );
}

/// \brief Get some simple PyMOL script text to draw the specified sec_file_record on the specified protein
///
/// An index is required to differentiate pseudo-atoms
string cath::sec::get_pymol_script_text(const protein         &prm_protein,         ///< The protein on which the secondary structure belongs
                                        const sec_file_record &prm_sec_file_record, ///< The sec_file_record to represent in PyMOL script
                                        const size_t          &prm_index            ///< The index of the sec_file_record (arbitrary but must be unique within a run of PyMOL)
                                        ) {
	const auto &midpoint         = prm_sec_file_record.get_midpoint();
	const line  the_line         = { midpoint, prm_sec_file_record.get_unit_dirn() };
	const auto &start_coord      = prm_protein.get_residue_ref_of_index( prm_sec_file_record.get_start_residue_num() ).get_carbon_alpha_coord();
	const auto &stop_coord       = prm_protein.get_residue_ref_of_index( prm_sec_file_record.get_stop_residue_num () ).get_carbon_alpha_coord();
	const auto  closest_to_start = closest_point_on_line_to_point( the_line, start_coord );
	const auto  closest_to_stop  = closest_point_on_line_to_point( the_line, stop_coord  );
	const auto &the_ss           = prm_sec_file_record.get_type();
	const auto  ss_str           = to_string( the_ss );
	const auto  colour_str       = ( the_ss == sec_struc_type::ALPHA_HELIX ) ? "orange"s    :
	                               ( the_ss == sec_struc_type::BETA_STRAND ) ? "pink"s      :
	                                                                           "chocolate"s ;

	const auto id_str            = std::to_string( midpoint.get_x() * midpoint.get_y() * midpoint.get_z() );

	return
	      "pseudoatom " + ss_str + "start, resi=" + std::to_string( prm_index ) + ", color=black, pos=["
		+ std::to_string( closest_to_start.get_x() ) + ", " + std::to_string( closest_to_start.get_y() ) + ", " + std::to_string( closest_to_start.get_z() ) + "]\n"
		+ "pseudoatom " + ss_str + "end,   resi=" + std::to_string( prm_index ) + ", color=" + colour_str + ", pos=["
		+ std::to_string( closest_to_stop .get_x() ) + ", " + std::to_string( closest_to_stop .get_y() ) + ", " + std::to_string( closest_to_stop .get_z() ) + "]\n"
		+ "distance " + ss_str + "line, " + ss_str + "start///" + std::to_string( prm_index ) + "/, " + ss_str + "end///" + std::to_string( prm_index ) + "/";
}

/// \brief Get some simple PyMOL script text to draw the specified sec_file_records on the specified protein
///
/// An index is required to differentiate pseudo-atoms
string cath::sec::get_pymol_script_text(const protein             &prm_protein,         ///< The protein on which the secondary structures belong
                                        const sec_file_record_vec &prm_sec_file_records ///< The sec_file_records to represent in PyMOL script
                                        ) {
	return join(
		indices( prm_sec_file_records.size() )
			| transformed( [&] (const size_t &x) {
				return get_pymol_script_text( prm_protein, prm_sec_file_records[ x ], x );
			} ),
		"\n" ) + R"(
hide
show_as lines, ( name n,ca,c )
show_as spheres,  Hstart
show_as spheres,  Hend
show_as spheres,  Sstart
show_as spheres,  Send
show_as dashes,   Hline
show_as dashes,   Sline
set sphere_scale, 0.5
set line_width,   2.5
set dash_width,   0.7
set dash_radius,  0.1
set dash_gap,     0.0
set dash_color,   black
set label_size,   25
label n. CA, resi
)";
}

// Notes:
// sec_struc_type (structure/protein/sec_struc_type.hpp) currently has three types:
//  * ALPHA_HELIX
//  * BETA_STRAND
//  * COIL
//
// The stream operators for that use H and S (and " " for output) respectively
//
// Sec files use:
//  * "H"
//  * "S"
//
// DSSP files use:
//  * H = α-helix
//  * B = residue in isolated β-bridge
//  * E = extended strand, participates in β ladder
//  * G = 3-helix (310 helix)
//  * I = 5 helix (π-helix)
//  * T = hydrogen bonded turn
//  * S = bend
//
// (of which prosec appears to use E (as S), G, H and I (as H))
//
// 1a8hA has examples of all of these types
//
// cath-tools' parse_dssp_residue_line() parses DSSP data in with:
//   sstruc_code == 'H' ? sec_struc_type::ALPHA_HELIX
// : sstruc_code == 'E' ? sec_struc_type::BETA_STRAND
//                      : sec_struc_type::COIL;
//


// pseudoatom secpointstart, color=tv_blue, pos=[24.495, 12.602, 15.071]
// pseudoatom secpointend,   color=tv_blue, pos=[25.860,  6.308, 13.730]
// distance secline, secpointstart, secpointend
// pseudoatom corepoint, color=red, pos=[24.464816, 12.605688, 15.025432]
// pseudoatom corepoint, color=red, pos=[25.248030,  9.404682, 14.483892]
// pseudoatom corepoint, color=red, pos=[25.828824,  6.310958, 13.682542]
// hide
// show_as lines, (name n,ca,c)
// show_as spheres, secpointstart
// show_as spheres, secpointend
// show_as spheres, corepoint
// show_as dashes, secline
// set sphere_scale, 0.01
// set dash_width,  0.7
// set dash_radius, 0.0025
// set dash_gap, 0.0
// set dash_color, black
// label n. CA, resi


// Example of 1hdoA00 (for which residue names match residue indices):
// 
// Res DSSP Guess Sec Flag
// =======================
//   1                    
//   2                    
//   3                    
//   4                    
//   5   E    E    E      
//   6   E    E    E      
//   7   E    E    E      
//   8   E    E    E      
//   9   E    E    E      
//  10   S                
//  11   T                
//  12   T                
//  13   S                
//  14   H    H    H      
//  15   H    H    H      
//  16   H    H    H      
//  17   H    H    H      
//  18   H    H    H      
//  19   H    H    H      
//  20   H    H    H      
//  21   H    H    H      
//  22   H    H    H      
//  23   H    H    H      
//  24   H    H    H      
//  25   H    H    H      
//  26   T                
//  27   T                
//  28                    
//  29   E    E    E      
//  30   E    E    E      
//  31   E    E    E      
//  32   E    E    E      
//  33   E    E    E      
//  34   E    E    E      
//  35   S                
//  36                    
//  37   G                
//  38   G                
//  39   G                
//  40   S                
//  41                    
//  42   S                
//  43   S                
//  44   S                
//  45                    
//  46                    
//  47                    
//  48   S                
//  49   E    E    E      
//  50   E    E    E      
//  51   E    E    E      
//  52   E    E    E      
//  53   S                
//  54                    
//  55   T                
//  56   T                
//  57   S                
//  58   H    H    H      
//  59   H    H    H      
//  60   H    H    H      
//  61   H    H    H      
//  62   H    H    H      
//  63   H    H    H      
//  64   H    H    H      
//  65   H    H    H      
//  66   T                
//  67   T                
//  68                    
//  69   S                
//  70   E    E    E      
//  71   E    E    E      
//  72   E    E    E      
//  73   E    E    E      
//  74                    
//  75                    
//  76                    
//  77                    
//  78   T                
//  79   T                
//  80                    
//  81                    
//  82   S                
//  83                    
//  84                    
//  85                    
//  86   H    H    H      
//  87   H    H    H      
//  88   H    H    H      
//  89   H    H    H      
//  90   H    H    H      
//  91   H    H    H      
//  92   H    H    H      
//  93   H    H    H      
//  94   H    H    H      
//  95   H    H    H      
//  96   H    H    H      
//  97   H    H    H      
//  98   H    H    H      
//  99   H    H    H      
// 100   H    H    H      
// 101   H    H    H      
// 102   T                
// 103                    
// 104                    
// 105   E    E    E      
// 106   E    E    E      
// 107   E    E    E      
// 108   E    E    E      
// 109   E    E    E      
// 110                    
// 111                    
// 112   G                
// 113   G                
// 114   G                
// 115   T                
// 116   S                
// 117                    
// 118   T                
// 119   T                
// 120                    
// 121   S                
// 122                    
// 123   G                
// 124   G                
// 125   G                
// 126   H    H    H      
// 127   H    H    H      
// 128   H    H    H      
// 129   H    H    H      
// 130   H    H    H      
// 131   H    H    H      
// 132   H    H    H      
// 133   H    H    H      
// 134   H    H    H      
// 135   H    H    H      
// 136   H    H    H      
// 137   H    H    H      
// 138   H    H    H      
// 139   H    H    H      
// 140   H    H    H      
// 141   H    H    H      
// 142   T                
// 143                    
// 144   S                
// 145   E    E    E      
// 146   E    E    E      
// 147   E    E    E      
// 148   E    E    E      
// 149   E    E    E      
// 150                    
// 151                    
// 152   S         E  ****
// 153   E    E    E      
// 154   E    E    E      
// 155   E    E    E      
// 156             E  ****
// 157                    
// 158                    
// 159                    
// 160                    
// 161   S                
// 162                    
// 163             E  ****
// 164   E    E    E      
// 165   E    E    E      
// 166   E    E    E      
// 167   S         E  ****
// 168   S                
// 169   S                
// 170                    
// 171   S                
// 172   S                
// 173                    
// 174   S          E ****
// 175   E    E     E     
// 176   E    E     E     
// 177   E    E     E     
// 178   H    H    HE ****
// 179   H    H    H      
// 180   H    H    H      
// 181   H    H    H      
// 182   H    H    H      
// 183   H    H    H      
// 184   H    H    H      
// 185   H    H    H      
// 186   H    H    H      
// 187   H    H    H      
// 188   T                
// 189   T                
// 190   S                
// 191                    
// 192   S                
// 193   T                
// 194   T                
// 195   T                
// 196   T                
// 197                    
// 198   E    E    E      
// 199   E    E    E      
// 200   E    E    E      
// 201   E    E    E      
// 202   E    E    E      
// 203                    
// 204                    
// 205                    
// =======================
// Res DSSP Guess Sec Flag