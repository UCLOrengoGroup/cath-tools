/// \file
/// \brief The residue class definitions

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

#include "residue.hpp"

#include <iomanip>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/throw_exception.hpp>

#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include "cath/common/difference.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/structure/protein/sec_struc_type.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::geom;

using ::boost::lexical_cast;
using ::boost::numeric_cast;
using ::std::get;
using ::std::ostream;
using ::std::ostringstream;
using ::std::right;
using ::std::setw;
using ::std::string;
using ::std::string_view;
using ::std::tuple;

#include "cath/structure/geometry/quat_rot.hpp"

// {
// 	// Useful information breaking down the 232/240 bytes currently required for a residue
// 	template <size_t T> class SD;

// 	struct backward_residue_id final {
// 		residue_name a;
// 		chain_label b;
// 	};

// 	struct broken_residue_id final {
// 		int res_num = 0;
// 		::std::optional<char> insert;
// 		bool is_null_residue_name;
// 		chain_label b;
// 	};

// 	SD< sizeof( residue               ) > size_of_residue;
// 	SD< sizeof( residue_id            ) > size_of_residue_id;
// 	SD< sizeof( amino_acid            ) > size_of_amino_acid;
// 	SD< sizeof( geom::coord           ) > size_of_coord;
// 	SD< sizeof( size_t                ) > size_of_size_t;
// 	SD< sizeof( sec_struc_type        ) > size_of_sec_struc_type;
// 	SD< sizeof( rotation              ) > size_of_rotation;
// 	SD< sizeof( doub_angle            ) > size_of_doub_angle;

// 	SD< sizeof( quat_rot_impl<double> ) > size_of_doub_quat_rot;

// 	SD< sizeof( residue_name          ) > size_of_residue_name;
// 	SD< sizeof( chain_label           ) > size_of_chain_label;
// 	SD< sizeof( int                   ) > size_of_int;
// 	SD< sizeof( optional<char>        ) > size_of_char_opt;
// 	SD< sizeof( bool                  ) > size_of_bool;
// 	SD< sizeof( backward_residue_id   ) > size_of_backward_residue_id;
// 	SD< sizeof( broken_residue_id     ) > size_of_broken_residue_id;
// }

/// Total solvent accessible surface areas
constexpr int ACCESSI( const char &x ) {
	// clang-format off
	switch ( x ) {
		case( 'A' ) : { return  85; }
		case( 'B' ) : { return 179; }
		case( 'C' ) : { return 139; }
		case( 'D' ) : { return 153; }
		case( 'E' ) : { return 208; }
		case( 'F' ) : { return 219; }
		case( 'G' ) : { return  83; }
		case( 'H' ) : { return 198; }
		case( 'I' ) : { return 170; }
		case( 'J' ) : { return 999; }
		case( 'K' ) : { return 215; }
		case( 'L' ) : { return 182; }
		case( 'M' ) : { return 199; }
		case( 'N' ) : { return 179; }
		case( 'O' ) : { return 999; }
		case( 'P' ) : { return 110; }
		case( 'Q' ) : { return 203; }
		case( 'R' ) : { return 259; }
		case( 'S' ) : { return 131; }
		case( 'T' ) : { return 172; }
		case( 'U' ) : { return 999; }
		case( 'V' ) : { return 151; }
		case( 'W' ) : { return 249; }
		case( 'X' ) : { return 999; }
		case( 'Y' ) : { return 246; }
		case( 'Z' ) : { return 208; }
		default : {
			BOOST_THROW_EXCEPTION( out_of_range_exception( "Unhandled char" ) );
		}
	}
	// clang-format on
}

///// \brief Throw if amino acid is invalid
//void residue::check_amino_acid() const {
//	const char amino_acid = get_amino_acid();
//	if (amino_acid != '0' && (!isalpha(amino_acid) || !isupper(amino_acid))) {
//		BOOST_THROW_EXCEPTION(invalid_argument_exception(
//			"Residue's amino acid character '"
//			+ string{ amino_acid }
//			+ "' is not an alpha character"
//		));
//	}
//}

/// \brief TODOCUMENT
///
/// \relates residue
string cath::ssap_legacy_alignment_left_side_string(const residue &prm_residue ///< TODOCUMENT
                                                    ) {
	ostringstream aln_ss;
	aln_ss        <<    setw(4) << pdb_number_or_value_if_null( prm_residue, 0 );
	aln_ss << " " << ( prm_residue.get_sec_struc_type() != sec_struc_type::COIL ? lexical_cast<string>( prm_residue.get_sec_struc_type() ) : "0" );
	aln_ss << " " <<   pdb_insert_or_value_if_null_or_absent( prm_residue, '0' );
	aln_ss << " " <<   get_amino_acid_letter_tolerantly( prm_residue );
	return aln_ss.str();
}

/// \brief Return the SSAP legacy alignment left-side string for a gap
///
/// \relates residue
string cath::ssap_legacy_alignment_left_side_gap_string() {
	return "   0 0 0 0";
}

/// \brief TODOCUMENT
///
/// \relates residue
string cath::ssap_legacy_alignment_right_side_string(const residue &prm_residue ///< TODOCUMENT
                                                     ) {
	ostringstream aln_ss;
	aln_ss        <<   get_amino_acid_letter_tolerantly( prm_residue );
	aln_ss << " " <<   pdb_insert_or_value_if_null_or_absent( prm_residue, '0' );
	aln_ss << " " << ( prm_residue.get_sec_struc_type() != sec_struc_type::COIL ? lexical_cast<string>( prm_residue.get_sec_struc_type() ) : "0" );
	aln_ss << " " <<   setw(4) << pdb_number_or_value_if_null( prm_residue, 0 );
	return aln_ss.str();
}

/// \brief Return the SSAP legacy alignment right-side string for a gap
///
/// \relates residue
string cath::ssap_legacy_alignment_right_side_gap_string() {
	return "0 0 0    0";
}

/// \brief TODOCUMENT
///
/// \relates residue
string cath::get_pdb_residue_id_string(const residue &prm_residue ///< TODOCUMENT
                                       ) {
	return to_string( prm_residue.get_pdb_residue_id() );
}

/// \brief TODOCUMENT
///
/// \relates residue
int cath::get_accessi_of_residue(const residue &prm_residue ///< TODOCUMENT
                                 ) {
	const char one_letter_amino_acid = get_amino_acid_letter_tolerantly( prm_residue );
	return ACCESSI( one_letter_amino_acid ) - numeric_cast<int>( prm_residue.get_access() );
}

/// \brief TODOCUMENT
///
/// \relates residue
ostream & cath::operator<<(ostream       &prm_os,     ///< TODOCUMENT
                           const residue &prm_residue ///< TODOCUMENT
                           ) {
	prm_os << "residue[" << right << setw( 4 ) << pdb_number( prm_residue );
	prm_os                                     << pdb_insert_or_value_if_null_or_absent( prm_residue, ' ' );
	prm_os << ", SS:"                          << prm_residue.get_sec_struc_type();
	prm_os << " ("       << right << setw( 2 ) << prm_residue.get_sec_struc_number();
	prm_os << "), AA:"                         << get_amino_acid_letter_tolerantly( prm_residue );
	prm_os << ", ACC:"   << right << setw( 2 ) << prm_residue.get_access();
	prm_os << ", PHI:"   << right << setw( 4 ) << prm_residue.get_phi_angle();
	prm_os << ", PSI:"   << right << setw( 4 ) << prm_residue.get_psi_angle();
	prm_os << ", CA "                          << prm_residue.get_carbon_alpha_coord();
	prm_os << ", CB "                          << prm_residue.get_carbon_beta_coord();
	prm_os << ", Frame "                       << prm_residue.get_frame();
	prm_os << "]";
	return prm_os;
}

/// \brief TODOCUMENT
rotation cath::construct_residue_frame(const coord &prm_nitrogen_position,     ///< TODOCUMENT
                                       const coord &prm_carbon_alpha_position, ///< TODOCUMENT
                                       const coord &prm_carbon_position        ///< TODOCUMENT
                                       ) {
	const coord unit_n_to_ca = normalise_copy( prm_carbon_alpha_position - prm_nitrogen_position );
	const coord unit_n_to_c  = normalise_copy( prm_carbon_position       - prm_nitrogen_position );
	return rotation_to_x_axis_and_x_y_plane(unit_n_to_c, cross_product( unit_n_to_ca, unit_n_to_c ) );

//	const coord unit_perp_to_n_ca_c_plane           = normalise_copy( cross_product( unit_n_to_ca, unit_n_to_c               ) );
//	const coord unit_in_n_ca_c_plane_perp_to_n_to_c = normalise_copy( cross_product( unit_n_to_c,  unit_perp_to_n_ca_c_plane ) );
//	cerr << "Unit N->C                                     : " << unit_n_to_c                         << endl;
//	cerr << "Unit perpendicular to N-C-CA plane            : " << unit_perp_to_n_ca_c_plane           << endl;
//	cerr << "Unit in N-C-CA plane, perpendicular to N to C : " << unit_in_n_ca_c_plane_perp_to_n_to_c << endl;
//	cerr << "Rotation matrix                               : " << bob                                 << endl;
//	cerr << endl;
}

/// \brief Combine two residue objects representing the same residue as parsed from a DSSP file and PDB file
///
/// This populates the accessibility and secondary structure from the DSSP and everything else from the PDB
residue cath::combine_residues_from_dssp_and_pdb(const residue                  &prm_dssp_residue,       ///< The residue that has been parsed from a DSSP file,
                                                 const residue                  &prm_pdb_residue,        ///< An equivalent residue that has been parsed from a pdb file (and converted to a protein object via)
                                                 const dssp_skip_angle_skipping &prm_pdb_skipped_angles, ///< Whether the PDB tried to break PHI/PSI angles like DSSP (so it's worth printing info where that's failed)
                                                 const residue_makeup           &prm_residue_makeup      ///< The makeup of the PDB residue (so that when this is residue_makeup::SOME_NON_PROPER_AMINO_ACIDS the function knows to silently tolerate a DSSP AA of 'X' where the residue from the PDB has a sensible AA)
                                                 ) {
	// Check that these residues match (or that the DSSP residue is an error residue)
	const residue_id &dssp_residue_id = prm_dssp_residue.get_pdb_residue_id();
	const residue_id &pdb_residue_id  = prm_pdb_residue.get_pdb_residue_id();
	if ( dssp_residue_id != pdb_residue_id && prm_dssp_residue != NULL_RESIDUE ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"When combining a DSSP residue with a PDB residue, DSSP residue \""
			+ to_string( dssp_residue_id )
			+ "\" does not match PDB residue \""
			+ to_string( pdb_residue_id  )
			+ "\""
		));
	}

	// Check that the two residues' amino acids match
	const amino_acid dssp_amino_acid = prm_dssp_residue.get_amino_acid();
	const amino_acid pdb_amino_acid  = prm_pdb_residue.get_amino_acid();
	if ( dssp_amino_acid != pdb_amino_acid ) {
		// If DSSP is UNK/X and PDB is ASX/B, GLX/Z, PYL/O or SEC/U or a non-proper amino acid from
		// HETATM records, then just accept it (and go with the PDB decision)
		//
		// \todo Should this also handle XLE/J in the same way?
		const bool dssp_is_unk   = ( dssp_amino_acid == amino_acid{ 'X' } );
		const bool pdb_is_proper = is_proper_amino_acid( pdb_amino_acid );
		const bool pdb_is_asx    = ( pdb_amino_acid  == amino_acid{ 'B' } );
		const bool pdb_is_glx    = ( pdb_amino_acid  == amino_acid{ 'Z' } );
		const bool pdb_is_pyl    = ( pdb_amino_acid  == amino_acid{ 'O' } );
		const bool pdb_is_sec    = ( pdb_amino_acid  == amino_acid{ 'U' } );

		if ( dssp_is_unk && ! pdb_is_proper ) {
			::spdlog::debug( R"(The amino acid "{}" parsed from a PDB for residue "{}" does not match amino)"
			                 R"( acid "{}" parsed from a DSSP but this is fine because it's a HETATM amino acid)",
			                 get_code_string( pdb_amino_acid ),
			                 to_string( pdb_residue_id ),
			                 get_code_string( dssp_amino_acid ) );
		}
		else if ( dssp_is_unk && ( pdb_is_pyl || pdb_is_sec || pdb_is_asx || pdb_is_glx ) ) {
			::spdlog::info( R"(The amino acid "{}" parsed from a PDB for residue "{}" does not match amino acid "{}" )"
			                "parsed from a DSSP - continuing because it looks like nothing worse than DSSP not "
			                "recognising a non-standard amino-acid code",
			                get_code_string( pdb_amino_acid ),
			                to_string( pdb_residue_id ),
			                get_code_string( dssp_amino_acid ) );
		}
		else if ( ! dssp_is_unk || prm_residue_makeup == residue_makeup::ALL_PROPER_AMINO_ACIDS ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Amino acid \""
				+ get_code_string( dssp_amino_acid )
				+ "\" parsed from DSSP for residue \""
				+ to_string( pdb_residue_id )
				+ "\" does not match amino acid \""
				+ get_code_string( pdb_amino_acid )
				+ "\" parsed from PDB"
			));
		}
	}

	// Check that the PDB and DSSP give PHI/PSI angles within 1 degree of each other
	//
	// \todo Simplify this code by creating an enum class that represents PHI/PSI
	//       and that can be used as a get_phi_or_psi_angle() function on a residue
	//       and then just loop over those two enum class values.
	using str_ang_ang_tpl = tuple<string, doub_angle, doub_angle>;
	const auto phi_psi_data = {
		str_ang_ang_tpl{ "phi", prm_pdb_residue.get_phi_angle(), prm_dssp_residue.get_phi_angle() },
		str_ang_ang_tpl{ "psi", prm_pdb_residue.get_psi_angle(), prm_dssp_residue.get_psi_angle() },
	};
	for (const auto &phi_psi_entry : phi_psi_data) {
		const string     &angle_name = get<0>( phi_psi_entry );
		const doub_angle &pdb_angle  = get<1>( phi_psi_entry );
		const doub_angle &dssp_angle = get<2>( phi_psi_entry );
		if ( angle_in_degrees( wrapped_difference( pdb_angle, dssp_angle ) ) > 1.0 ) {
			if ( angle_in_degrees( dssp_angle ) == 360.0 ) {
				if ( prm_pdb_skipped_angles == dssp_skip_angle_skipping::BREAK_ANGLES ) {
					::spdlog::info( "The {} angle calculated for residue {} ({}) from the PDB conflicts with the DSSP "
					                "angle of 360.0, which is typically used by DSSP where it detects a break in the "
					                "chain (perhaps because it has rejected a neighbouring residue)",
					                angle_name,
					                pdb_residue_id,
					                pdb_angle );
				}
			}
			else {
				::spdlog::warn(
				  "Whilst combining PDB and DSSP files, at residue {} detected conflicting {} angles: {} and {}",
				  pdb_residue_id,
				  angle_name,
				  pdb_angle,
				  dssp_angle );
			}
		}
	}

	// Check that the PDB and DSSP give x/y/z coordinate values within 0.5001 of each other
	using char_doub_doub_tpl = tuple<char, double, double>;
	const auto &pdb_ca_coord  = prm_pdb_residue.get_carbon_alpha_coord();
	const auto &dssp_ca_coord = prm_dssp_residue.get_carbon_alpha_coord();
	const auto coord_data = {
		char_doub_doub_tpl{ 'x', pdb_ca_coord.get_x(), dssp_ca_coord.get_x() },
		char_doub_doub_tpl{ 'y', pdb_ca_coord.get_y(), dssp_ca_coord.get_y() },
		char_doub_doub_tpl{ 'z', pdb_ca_coord.get_z(), dssp_ca_coord.get_z() },
	};
	for (const auto &x_y_z_entry : coord_data) {
		const char   &dim        = get<0>( x_y_z_entry );
		const double &pdb_value  = get<1>( x_y_z_entry );
		const double &dssp_value = get<2>( x_y_z_entry );
		if ( difference( pdb_value, dssp_value ) > 0.05001 ) {
			::spdlog::warn(
			  "Whilst combining PDB and DSSP files, at residue {} detected conflicting {} coordinate: {} and {}",
			  pdb_residue_id,
			  dim,
			  pdb_value,
			  dssp_value );
		}
	}

	// Return a new residue that takes the accessibility and secondary struct from the DSSP residue
	// and everything else from the PDB residue
	const residue new_residue(
		pdb_residue_id,
		prm_pdb_residue.get_amino_acid(), //< Use the PDB's AA (see notes about PYL/SEC above)
		prm_pdb_residue.get_carbon_alpha_coord(),
		prm_pdb_residue.get_carbon_beta_coord(),
		prm_dssp_residue.get_sec_struc_number(),
		prm_dssp_residue.get_sec_struc_type(),
		prm_pdb_residue.get_frame(),
		prm_pdb_residue.get_phi_angle(),
		prm_pdb_residue.get_psi_angle(),
		prm_dssp_residue.get_access()
	);

	// Temporary debug statements
//	cerr << "DSSP residue : " << prm_dssp_residue << endl;
//	cerr << " PDB residue : " << prm_pdb_residue  << endl;
//	cerr << " new residue : " << new_residue      << endl;

	return new_residue;
}
