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

#include "residue.h"

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/throw_exception.hpp>

#include "exception/invalid_argument_exception.h"
#include "structure/protein/sec_struc_type.h"

#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <utility>


using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::lexical_cast;
using boost::math::isnormal;
using boost::numeric_cast;

/// \brief TODOCUMENT
///
/// Total solvent accessible surface areas
const map<char, int> ACCESSI = {
	{'A',  85}, {'B', 179}, {'C', 139}, {'D', 153}, {'E', 208},
	{'F', 219}, {'G',  83}, {'H', 198}, {'I', 170}, {'J', 999},
	{'K', 215}, {'L', 182}, {'M', 199}, {'N', 179}, {'O', 999},
	{'P', 110}, {'Q', 203}, {'R', 259}, {'S', 131}, {'T', 172},
	{'U', 999}, {'V', 151}, {'W', 249}, {'X', 999}, {'Y', 246},
	{'Z', 208}
};

/// \brief TODOCUMENT
const residue residue::NULL_RESIDUE(
	residue_name(  0  ),
	amino_acid  ( 'X' ),
	coord::ORIGIN_COORD,
	coord::ORIGIN_COORD,
	0, // Secondary structure number
	sec_struc_type::COIL,
	rotation::IDENTITY_ROTATION(),
	DEFAULT_PHI_PSI(),
	DEFAULT_PHI_PSI(),
	0 // Accessibility
);

/// \brief Throw if angle is not in range [0, 360]
void residue::check_phi_psi_angle(const doub_angle &arg_angle
                                  ) {
//	using boost::math::isfinite;
	if ( arg_angle < zero_angle<double>() || arg_angle > one_revolution<double>() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Phi/psi angle ("
			+ lexical_cast<string>( arg_angle )
			+ ") is not in valid range [0, 360]"
		));
	}
}

///// \brief Throw if amino acid is invalid
//void residue::check_amino_acid() const {
//	const char amino_acid = get_amino_acid();
//	if (amino_acid != '0' && (!isalpha(amino_acid) || !isupper(amino_acid))) {
//		BOOST_THROW_EXCEPTION(invalid_argument_exception(
//			"Residue's amino acid character '"
//			+ string(1, amino_acid)
//			+ "' is not an alpha character"
//		));
//	}
//}

/// \brief Ctor for residue.
residue::residue(const residue_name   &arg_residue_name,       ///< TODOCUMENT
                 const amino_acid     &arg_amino_acid,         ///< TODOCUMENT
                 const coord          &arg_carbon_alpha_coord, ///< Coordinates of the carbon alpha atom
                 const coord          &arg_carbon_beta_coord,  ///< Coordinates of the carbon beta atom
                 const size_t         &arg_sec_struc_number,   ///< TODOCUMENT
                 const sec_struc_type &arg_sec_struc_type,     ///< TODOCUMENT
                 const rotation       &arg_frame,              ///< TODOCUMENT
                 const doub_angle     &arg_phi,                ///< TODOCUMENT
                 const doub_angle     &arg_psi,                ///< TODOCUMENT
                 const size_t         &arg_access              ///< TODOCUMENT
                 ) : the_residue_name  ( arg_residue_name       ),
                     the_amino_acid    ( arg_amino_acid         ),
                     carbon_alpha_coord( arg_carbon_alpha_coord ),
                     carbon_beta_coord ( arg_carbon_beta_coord  ),
                     sec_struc_number  ( arg_sec_struc_number   ),
                     the_sec_struc_type( arg_sec_struc_type     ),
                     frame             ( arg_frame              ),
                     phi_angle         ( arg_phi                ),
                     psi_angle         ( arg_psi                ),
                     access            ( arg_access             ) {
	check_phi_psi_angle( get_phi_angle() );
	check_phi_psi_angle( get_psi_angle() );
}

/// \brief TODOCUMENT
void residue::set_amino_acid(const amino_acid &arg_amino_acid ///< TODOCUMENT
                             ) {
	the_amino_acid = arg_amino_acid;
}

/// \brief TODOCUMENT
void residue::set_residue_sec_struc_number(const size_t &arg_sec_struc_number ///< TODOCUMENT
                                           ) {
	sec_struc_number = arg_sec_struc_number;
}

/// \brief TODOCUMENT
void residue::set_sec_struc_type(const sec_struc_type &arg_sec_struc_type ///< TODOCUMENT
                                 ) {
	the_sec_struc_type = arg_sec_struc_type;
}

/// \brief TODOCUMENT
const residue_name & residue::get_pdb_residue_name() const {
	return the_residue_name;
}

/// \brief TODOCUMENT
amino_acid residue::get_amino_acid() const {
	return the_amino_acid;
}

/// \brief TODOCUMENT
size_t residue::get_sec_struc_number() const {
	return sec_struc_number;
}

/// \brief TODOCUMENT
sec_struc_type residue::get_sec_struc_type() const {
	return the_sec_struc_type;
}

/// \brief TODOCUMENT
const doub_angle & residue::get_phi_angle() const {
	return phi_angle;
}

/// \brief TODOCUMENT
const doub_angle & residue::get_psi_angle() const {
	return psi_angle;
}

/// \brief TODOCUMENT
size_t residue::get_access() const {
	return access;
}

/// \brief The default angles in degrees that's used for undetermined phi/psi angles
///
/// \todo Shouldn't the undetermined angles be handled more explicitly?
///       The DSSP/WOLF files have the angles in the more natural range of [-180, 180]
///       and the currently always get shifted into the (0, 360] range
///       if they were left alone, that would leave 360.0 as a special "undetermined" value
///       but this would require changes in residues_have_similar_area_angle_props()
doub_angle residue::DEFAULT_PHI_PSI() {
	return one_revolution<double>();
}

/// \brief Equality operator for residue
///
/// operator!= is provided by boost::equality_comparable
///
/// \relates residue
bool cath::operator==(const residue &arg_residue_1, ///< The first  residue to compare
                      const residue &arg_residue_2  ///< The second residue to compare
                      ) {
	// Search for possible differences
	if ( arg_residue_1.get_pdb_residue_name()     != arg_residue_2.get_pdb_residue_name()     ) {
		return false;
	}
	if ( arg_residue_1.get_amino_acid()           != arg_residue_2.get_amino_acid()           ) {
		return false;
	}
	if ( arg_residue_1.get_carbon_alpha_coord()   != arg_residue_2.get_carbon_alpha_coord()   ) {
		return false;
	}
	if ( arg_residue_1.get_carbon_beta_coord()    != arg_residue_2.get_carbon_beta_coord()    ) {
		return false;
	}
	if ( arg_residue_1.get_sec_struc_number()     != arg_residue_2.get_sec_struc_number()     ) {
		return false;
	}
	if ( arg_residue_1.get_sec_struc_type()       != arg_residue_2.get_sec_struc_type()       ) {
		return false;
	}
	if ( arg_residue_1.get_frame()                != arg_residue_2.get_frame()                ) {
		return false;
	}
	if ( arg_residue_1.get_phi_angle() != arg_residue_2.get_phi_angle() ) {
		return false;
	}
	if ( arg_residue_1.get_psi_angle() != arg_residue_2.get_psi_angle() ) {
		return false;
	}
	if ( arg_residue_1.get_access()               != arg_residue_2.get_access()               ) {
		return false;
	}

	// If no differences have been found, then return false
	return true;
}

/// \brief TODOCUMENT
///
/// \relates residue
int cath::get_pdb_name_number(const residue &arg_residue ///< TODOCUMENT
                              ) {
	return arg_residue.get_pdb_residue_name().get_residue_number();
}

/// \brief TODOCUMENT
///
/// \relates residue
bool cath::has_pdb_name_insert(const residue &arg_residue ///< TODOCUMENT
                               ) {
	return has_insert_code( arg_residue.get_pdb_residue_name() );
}

/// \brief TODOCUMENT
///
/// \relates residue
char cath::get_pdb_name_insert(const residue &arg_residue ///< TODOCUMENT
                               ) {
	if ( ! has_pdb_name_insert( arg_residue ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get PDB residue insert code of residue that has no insert code"));
	}
	return get_insert_code( arg_residue.get_pdb_residue_name() );
}

/// \brief TODOCUMENT
///
/// \relates residue
bool cath::residue_matches_residue_name(const residue      &arg_residue,     ///< TODOCUMENT
                                        const residue_name &arg_residue_name ///< TODOCUMENT
                                        ) {
	return ( arg_residue.get_pdb_residue_name() == arg_residue_name );
}

/// \brief TODOCUMENT
///
/// \relates residue
char cath::get_amino_acid_letter(const residue &arg_residue ///< TODOCUMENT
                                 ) {
	return arg_residue.get_amino_acid().get_letter();
}

/// \brief TODOCUMENT
///
/// \relates residue
string cath::get_amino_acid_code(const residue &arg_residue ///< TODOCUMENT
                                 ) {
	return arg_residue.get_amino_acid().get_code();
}

/// \brief TODOCUMENT
///
/// \relates residue
string cath::get_amino_acid_name(const residue &arg_residue ///< TODOCUMENT
                                 ) {
	return arg_residue.get_amino_acid().get_name();
}

/// \brief TODOCUMENT
///
/// \relates residue
string cath::ssap_legacy_alignment_left_side_string(const residue &arg_residue ///< TODOCUMENT
                                                    ) {
	ostringstream aln_ss;
	aln_ss        <<    setw(4) << get_pdb_name_number( arg_residue );
	aln_ss << " " << ( arg_residue.get_sec_struc_type() != sec_struc_type::COIL ? lexical_cast<string>( arg_residue.get_sec_struc_type() ) : "0" );
	aln_ss << " " << ( has_pdb_name_insert  ( arg_residue )     ?                       get_pdb_name_insert( arg_residue ) : '0' );
	aln_ss << " " <<   get_amino_acid_letter( arg_residue );
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
string cath::ssap_legacy_alignment_right_side_string(const residue &arg_residue ///< TODOCUMENT
                                                     ) {
	ostringstream aln_ss;
	aln_ss        <<   get_amino_acid_letter( arg_residue );
	aln_ss << " " << ( has_pdb_name_insert  ( arg_residue )     ?                       get_pdb_name_insert( arg_residue ) : '0' );
	aln_ss << " " << ( arg_residue.get_sec_struc_type() != sec_struc_type::COIL ? lexical_cast<string>( arg_residue.get_sec_struc_type() ) : "0" );
	aln_ss << " " <<   setw(4) << get_pdb_name_number( arg_residue );
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
string cath::get_pdb_residue_name_string(const residue &arg_residue ///< TODOCUMENT
                                         ) {
	return lexical_cast<string>( arg_residue.get_pdb_residue_name() );
}

/// \brief TODOCUMENT
///
/// \relates residue
int cath::get_accessi_of_residue(const residue &arg_residue ///< TODOCUMENT
                                 ) {
	const char one_letter_amino_acid = get_amino_acid_letter( arg_residue );
	return ACCESSI.find( one_letter_amino_acid )->second - numeric_cast<int>( arg_residue.get_access() );
}

/// \brief TODOCUMENT
///
/// \relates residue
ostream & cath::operator<<(ostream       &arg_os,     ///< TODOCUMENT
                           const residue &arg_residue ///< TODOCUMENT
                           ) {
	arg_os << "residue[" << right << setw( 4 ) << get_pdb_name_number( arg_residue );
	arg_os                                     << ( has_pdb_name_insert( arg_residue ) ? get_pdb_name_insert( arg_residue ) : ' ' );
	arg_os << ", SS:"                          << arg_residue.get_sec_struc_type();
	arg_os << " ("       << right << setw( 2 ) << arg_residue.get_sec_struc_number();
	arg_os << "), AA:"                         << get_amino_acid_letter( arg_residue );
	arg_os << ", ACC:"   << right << setw( 2 ) << arg_residue.get_access();
	arg_os << ", PHI:"   << right << setw( 4 ) << arg_residue.get_phi_angle();
	arg_os << ", PSI:"   << right << setw( 4 ) << arg_residue.get_psi_angle();
	arg_os << ", CA "                          << arg_residue.get_carbon_alpha_coord();
	arg_os << ", CB "                          << arg_residue.get_carbon_beta_coord();
	arg_os << ", Frame "                       << arg_residue.get_frame();
	arg_os << "]";
	return arg_os;
}

/// \brief TODOCUMENT
rotation cath::construct_residue_frame(const coord &arg_nitrogen_position,     ///< TODOCUMENT
                                       const coord &arg_carbon_alpha_position, ///< TODOCUMENT
                                       const coord &arg_carbon_position        ///< TODOCUMENT
                                       ) {
	const coord unit_n_to_ca = normalise_copy( arg_carbon_alpha_position - arg_nitrogen_position );
	const coord unit_n_to_c  = normalise_copy( arg_carbon_position       - arg_nitrogen_position );
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
residue cath::combine_residues_from_dssp_and_pdb(const residue &arg_dssp_residue, ///< The residue that has been parsed from a DSSP file,
                                                 const residue &arg_pdb_residue   ///< An equivalent residue that has been parsed from a pdb file (and converted to a protein object via)
                                                 ) {
	// Check that these residues match (or that the DSSP residue is an error residue)
	const residue_name &dssp_residue_name = arg_dssp_residue.get_pdb_residue_name();
	const residue_name &pdb_residue_name  = arg_pdb_residue.get_pdb_residue_name();
	if ( dssp_residue_name != pdb_residue_name && arg_dssp_residue != residue::NULL_RESIDUE ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"When combining a DSSP residue with a PDB residue, DSSP residue \""
			+ lexical_cast<string>( dssp_residue_name )
			+ "\" does not match PDB residue \""
			+ lexical_cast<string>( pdb_residue_name  )
			+ "\""
		));
	}

	// Check that the two residues' amino acids match
	const amino_acid dssp_amino_acid = arg_dssp_residue.get_amino_acid();
	const amino_acid pdb_amino_acid  = arg_pdb_residue.get_amino_acid();
	if ( dssp_amino_acid != pdb_amino_acid ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Amino acid \""
			+ dssp_amino_acid.get_code()
			+ "\" parsed from DSSP for residue \""
			+ lexical_cast<string>( pdb_residue_name )
			+ "\" does not match amino acid \""
			+ pdb_amino_acid.get_code()
			+ "\" parsed from PDB"
		));
	}

	// Return a new residue that takes the accessibility and secondary struct from the DSSP residue
	// and everything else from the PDB residue
	const residue new_residue(
		pdb_residue_name,
		arg_pdb_residue.get_amino_acid(),
		arg_pdb_residue.get_carbon_alpha_coord(),
		arg_pdb_residue.get_carbon_beta_coord(),
		arg_dssp_residue.get_sec_struc_number(),
		arg_dssp_residue.get_sec_struc_type(),
		arg_pdb_residue.get_frame(),
		arg_pdb_residue.get_phi_angle(),
		arg_pdb_residue.get_psi_angle(),
		arg_dssp_residue.get_access()
	);

	// Temporary debug statements
//	cerr << "DSSP residue : " << arg_dssp_residue << endl;
//	cerr << " PDB residue : " << arg_pdb_residue  << endl;
//	cerr << " new residue : " << new_residue      << endl;

	return new_residue;
}

/// \brief Check whether the specified residue is a null residue
///
/// This currently does a full equality comparison with residue::NULL_RESIDUE but
/// could probably be made more efficient if there is call for that.
///
/// \relates residue
bool cath::is_null_residue(const residue &arg_residue ///< TODOCUMENT
                           ) {
	return (arg_residue == residue::NULL_RESIDUE);
}

/// \brief Wipe the residue's secondary structure number and set its label to sec_struc_type::COIL
///
/// \relates residue
void cath::wipe_secondary_structure(residue &arg_residue ///< The residue to be modified
                                    ) {
	arg_residue.set_sec_struc_type          ( sec_struc_type::COIL );
	arg_residue.set_residue_sec_struc_number( 0                    );
}
