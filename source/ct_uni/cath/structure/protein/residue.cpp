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

#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/throw_exception.hpp>

#include "cath/common/difference.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/structure/protein/sec_struc_type.hpp"

#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;

using boost::lexical_cast;
using boost::numeric_cast;
using std::get;
using std::map;
using std::ostream;
using std::ostringstream;
using std::right;
using std::setw;
using std::string;
using std::tuple;

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
// 		boost::optional<char> insert;
// 		bool is_null_residue_name;
// 		chain_label b;
// 	};

// 	using boost::optional;

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
	residue_id{},
	amino_acid  ( 'X' ),
	ORIGIN_COORD,
	ORIGIN_COORD,
	0, // Secondary structure number
	sec_struc_type::COIL,
	rotation::IDENTITY_ROTATION(),
	DEFAULT_PHI_PSI(),
	DEFAULT_PHI_PSI(),
	0 // Accessibility
);

/// \brief Throw if angle is not in range [0, 360]
void residue::check_phi_psi_angle(const doub_angle &prm_angle
                                  ) {
//	using boost::math::isfinite;
	if ( prm_angle < zero_angle<double>() || prm_angle > one_revolution<double>() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Phi/psi angle ("
			+ lexical_cast<string>( prm_angle )
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
//			+ string{ amino_acid }
//			+ "' is not an alpha character"
//		));
//	}
//}

/// \brief Ctor for residue.
residue::residue(const residue_id     &prm_residue_id,         ///< TODOCUMENT
                 amino_acid            prm_amino_acid,         ///< TODOCUMENT
                 coord                 prm_carbon_alpha_coord, ///< Coordinates of the carbon alpha atom
                 coord                 prm_carbon_beta_coord,  ///< Coordinates of the carbon beta atom
                 const size_t         &prm_sec_struc_number,   ///< TODOCUMENT
                 const sec_struc_type &prm_sec_struc_type,     ///< TODOCUMENT
                 rotation              prm_frame,              ///< TODOCUMENT
                 doub_angle            prm_phi,                ///< TODOCUMENT
                 doub_angle            prm_psi,                ///< TODOCUMENT
                 const size_t         &prm_access              ///< TODOCUMENT
                 ) : the_residue_id    ( prm_residue_id                       ),
                     the_amino_acid    ( std::move ( prm_amino_acid         ) ),
                     carbon_alpha_coord( std::move ( prm_carbon_alpha_coord ) ),
                     carbon_beta_coord ( std::move ( prm_carbon_beta_coord  ) ),
                     sec_struc_number  ( prm_sec_struc_number                 ),
                     the_sec_struc_type( prm_sec_struc_type                   ),
                     frame             ( std::move ( prm_frame              ) ),
                     phi_angle         ( std::move ( prm_phi                ) ),
                     psi_angle         ( std::move ( prm_psi                ) ),
                     access            ( prm_access                           ) {
	check_phi_psi_angle( get_phi_angle() );
	check_phi_psi_angle( get_psi_angle() );
}

/// \brief Setter for the amino acid
residue & residue::set_amino_acid(const amino_acid &prm_amino_acid ///< The amino acid to set
                                  ) {
	the_amino_acid = prm_amino_acid;
	return *this;
}

/// \brief Setter for the residue sec_struc number
residue & residue::set_residue_sec_struc_number(const size_t &prm_sec_struc_number ///< The residue sec_struc number to set
                                                ) {
	sec_struc_number = prm_sec_struc_number;
	return *this;
}

/// \brief Setter for the sec_struc_type
residue & residue::set_sec_struc_type(const sec_struc_type &prm_sec_struc_type ///< The sec_struc_type to set
                                      ) {
	the_sec_struc_type = prm_sec_struc_type;
	return *this;
}

/// \brief Setter for the accessibility (calculated in a DSSP/wolf manner)
residue & residue::set_access(const size_t &prm_accessibility ///< The accessibility (calculated in a DSSP/wolf manner) to set
                              ) {
	access = prm_accessibility;
	return *this;
}

/// \brief TODOCUMENT
const residue_id & residue::get_pdb_residue_id() const {
	return the_residue_id;
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
bool cath::operator==(const residue &prm_residue_1, ///< The first  residue to compare
                      const residue &prm_residue_2  ///< The second residue to compare
                      ) {
	return (
		( prm_residue_1.get_pdb_residue_id()     == prm_residue_2.get_pdb_residue_id()     )
		&&
		( prm_residue_1.get_amino_acid()         == prm_residue_2.get_amino_acid()         )
		&&
		( prm_residue_1.get_carbon_alpha_coord() == prm_residue_2.get_carbon_alpha_coord() )
		&&
		( prm_residue_1.get_carbon_beta_coord()  == prm_residue_2.get_carbon_beta_coord()  )
		&&
		( prm_residue_1.get_sec_struc_number()   == prm_residue_2.get_sec_struc_number()   )
		&&
		( prm_residue_1.get_sec_struc_type()     == prm_residue_2.get_sec_struc_type()     )
		&&
		( prm_residue_1.get_frame()              == prm_residue_2.get_frame()              )
		&&
		( prm_residue_1.get_phi_angle()          == prm_residue_2.get_phi_angle()          )
		&&
		( prm_residue_1.get_psi_angle()          == prm_residue_2.get_psi_angle()          )
		&&
		( prm_residue_1.get_access()             == prm_residue_2.get_access()             )
	);
}

/// \brief TODOCUMENT
///
/// \relates residue
const chain_label & cath::get_chain_label(const residue &prm_residue ///< TODOCUMENT
                                          ) {
	return prm_residue.get_pdb_residue_id().get_chain_label();
}

/// \brief TODOCUMENT
///
/// \relates residue
const residue_name & cath::get_pdb_residue_name(const residue &prm_residue ///< TODOCUMENT
                                                ) {
	return prm_residue.get_pdb_residue_id().get_residue_name();
}


/// \brief TODOCUMENT
///
/// \relates residue
const int & cath::pdb_number(const residue &prm_residue ///< TODOCUMENT
                             ) {
	return get_pdb_residue_name( prm_residue ).residue_number();
}

/// \brief TODOCUMENT
///
/// \relates residue
int cath::pdb_number_or_value_if_null(const residue &prm_residue, ///< TODOCUMENT
                                      const int     &prm_value
                                      ) {
	return residue_number_or_value_if_null( get_pdb_residue_name( prm_residue ), prm_value );
}

/// \brief TODOCUMENT
///
/// \relates residue
bool cath::has_pdb_insert(const residue &prm_residue ///< TODOCUMENT
                          ) {
	return has_insert( get_pdb_residue_name( prm_residue ) );
}

/// \brief TODOCUMENT
///
/// \relates residue
bool cath::has_pdb_insert_or_value_if_null(const residue &prm_residue, ///< TODOCUMENT
                                           const bool    &prm_value    ///< TODOCUMENT
                                           ) {
	return has_insert_or_value_if_null( get_pdb_residue_name( prm_residue ), prm_value );
}

/// \brief TODOCUMENT
///
/// \relates residue
const char & cath::pdb_insert(const residue &prm_residue ///< TODOCUMENT
                              ) {
	if ( ! has_pdb_insert( prm_residue ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get PDB residue insert code of residue that has no insert code"));
	}
	return insert( get_pdb_residue_name( prm_residue ) );
}

/// \brief TODOCUMENT
///
/// \relates residue
char cath::pdb_insert_or_value_if_null(const residue &prm_residue, ///< TODOCUMENT
                                       const char    &prm_value    ///< TODOCUMENT
                                       ) {
	return insert_or_value_if_null( get_pdb_residue_name( prm_residue ), prm_value );
}

/// \brief TODOCUMENT
///
/// \relates residue
char cath::pdb_insert_or_value_if_null_or_absent(const residue &prm_residue, ///< TODOCUMENT
                                                 const char    &prm_value    ///< TODOCUMENT
                                                 ) {
	return insert_or_value_if_null_or_absent( get_pdb_residue_name( prm_residue ), prm_value );
}


/// \brief TODOCUMENT
///
/// \relates residue
bool cath::residue_matches_residue_id(const residue    &prm_residue,   ///< TODOCUMENT
                                      const residue_id &prm_residue_id ///< TODOCUMENT
                                      ) {
	return ( prm_residue.get_pdb_residue_id() == prm_residue_id );
}

/// \brief TODOCUMENT
///
/// \relates residue
char cath::get_amino_acid_letter_tolerantly(const residue &prm_residue ///< TODOCUMENT
                                            ) {
	return prm_residue.get_amino_acid().get_letter_tolerantly();
}

/// \brief TODOCUMENT
///
/// \relates residue
char_3_arr cath::get_amino_acid_code(const residue &prm_residue ///< TODOCUMENT
                                     ) {
	return prm_residue.get_amino_acid().get_code();
}

/// \brief TODOCUMENT
///
/// \relates residue
string cath::get_amino_acid_name(const residue &prm_residue ///< TODOCUMENT
                                 ) {
	return prm_residue.get_amino_acid().get_name();
}

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
	return ACCESSI.find( one_letter_amino_acid )->second - numeric_cast<int>( prm_residue.get_access() );
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
	if ( dssp_residue_id != pdb_residue_id && prm_dssp_residue != residue::NULL_RESIDUE ) {
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
			BOOST_LOG_TRIVIAL( debug ) << "The amino acid \""
				<< get_code_string( pdb_amino_acid )
				<< "\" parsed from a PDB for residue \""
				<< to_string( pdb_residue_id )
				<< "\" does not match amino acid \""
				<< get_code_string( dssp_amino_acid )
				<< "\" parsed from a DSSP but this is fine because it's a HETATM amino acid";
		}
		else if ( dssp_is_unk && ( pdb_is_pyl || pdb_is_sec || pdb_is_asx || pdb_is_glx ) ) {
			BOOST_LOG_TRIVIAL( info ) << "The amino acid \""
				<< get_code_string( pdb_amino_acid )
				<< "\" parsed from a PDB for residue \""
				<< to_string( pdb_residue_id )
				<< "\" does not match amino acid \""
				<< get_code_string( dssp_amino_acid )
				<< "\" parsed from a DSSP - continuing because it looks "
				<< "like nothing worse than DSSP not recognising a non-standard amino-acid code";
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
					BOOST_LOG_TRIVIAL( info ) << "The "
						<< angle_name
						<< " angle calculated for residue "
						<< pdb_residue_id
						<< " ("
						<< pdb_angle
						<< ") from the PDB conflicts with the DSSP angle of 360.0"
						<< ", which is typically used by DSSP where it detects a"
						<< " break in the chain (perhaps because it has rejected a neighbouring residue)";
				}
			}
			else {
				BOOST_LOG_TRIVIAL( warning ) << "Whilst combining PDB and DSSP files, at residue "
					<< pdb_residue_id
					<< " detected conflicting "
					<< angle_name
					<< " angles: "
					<< pdb_angle
					<< " and "
					<< dssp_angle;
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
			BOOST_LOG_TRIVIAL( warning ) << "Whilst combining PDB and DSSP files, at residue "
				<< pdb_residue_id
				<< " detected conflicting "
				<< dim
				<< " coordinate: "
				<< pdb_value
				<< " and "
				<< dssp_value;
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

/// \brief Check whether the specified residue is a null residue
///
/// This currently does a full equality comparison with residue::NULL_RESIDUE but
/// could probably be made more efficient if there is call for that.
///
/// \relates residue
bool cath::is_null_residue(const residue &prm_residue ///< TODOCUMENT
                           ) {
	return (prm_residue == residue::NULL_RESIDUE);
}

/// \brief Wipe the residue's secondary structure number and set its label to sec_struc_type::COIL
///
/// \relates residue
void cath::wipe_secondary_structure(residue &prm_residue ///< The residue to be modified
                                    ) {
	prm_residue.set_sec_struc_type          ( sec_struc_type::COIL );
	prm_residue.set_residue_sec_struc_number( 0                    );
}
