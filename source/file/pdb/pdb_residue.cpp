/// \file
/// \brief The pdb_residue class definitions

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

#include "pdb_residue.h"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>

#include "common/cpp14/cbegin_cend.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb_atom.h"
#include "structure/protein/residue.h"

#include <cmath>

using namespace boost::math::constants;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::algorithm::any_of;
using boost::algorithm::join;
using boost::lexical_cast;
using boost::none;

/// \brief Ctor for pdb_residue
pdb_residue::pdb_residue(const chain_label  &arg_chain_label,  ///< TODOCUMENT
                         const residue_name &arg_residue_name, ///< TODOCUMENT
                         const pdb_atom_vec &arg_atoms         ///< TODOCUMENT
                         ) : the_chain_label ( arg_chain_label  ),
                             the_residue_name( arg_residue_name ),
                             atoms           ( arg_atoms        ) {
}

/// \brief Ctor for pdb_residue
pdb_residue::pdb_residue(const chain_label   &arg_chain_label,  ///< TODOCUMENT
                         const residue_name  &arg_residue_name, ///< TODOCUMENT
                         pdb_atom_vec       &&arg_atoms         ///< TODOCUMENT
                         ) : the_chain_label ( arg_chain_label        ),
                             the_residue_name( arg_residue_name       ),
                             atoms           ( std::move( arg_atoms ) ) {
}

/// \brief TODOCUMENT
chain_label pdb_residue::get_chain_label() const {
	return the_chain_label;
}

/// \brief TODOCUMENT
residue_name pdb_residue::get_residue_name() const {
	return the_residue_name;
}

/// \brief TODOCUMENT
size_t pdb_residue::get_num_atoms() const {
	return atoms.size();
}

/// \brief TODOCUMENT
const pdb_atom & pdb_residue::get_atom_cref_of_index(const size_t &arg_index ///< TODOCUMENT
                                                     ) const {
	return atoms[arg_index];
}

/// \brief TODOCUMENT
void pdb_residue::set_chain_label(const chain_label &arg_chain_label ///< TODOCUMENT
                                  ) {
	the_chain_label = arg_chain_label;
}

/// \brief TODOCUMENT
void pdb_residue::rotate(const rotation &arg_rotation ///< TODOCUMENT
                         ) {
	for (pdb_atom &my_pdb_atom : atoms) {
		my_pdb_atom.rotate(arg_rotation);
	}
}

/// \brief TODOCUMENT
void pdb_residue::operator+=(const coord &arg_coord ///< TODOCUMENT
                             ) {
	for (pdb_atom &my_pdb_atom : atoms) {
		my_pdb_atom += arg_coord;
	}
}

/// \brief TODOCUMENT
void pdb_residue::operator-=(const coord &arg_coord ///< TODOCUMENT
                             ) {
	for (pdb_atom &my_pdb_atom : atoms) {
		my_pdb_atom -= arg_coord;
	}
}

/// \brief TODOCUMENT
pdb_residue::const_iterator pdb_residue::begin() const {
	return common::cbegin( atoms );
}

/// \brief TODOCUMENT
pdb_residue::const_iterator pdb_residue::end() const {
	return common::cend( atoms );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
const amino_acid & cath::file::get_amino_acid(const pdb_residue &arg_residue ///< The pdb_residue to query
                                              ) {
	// If any proper amino acids can be found, then return the first such one
	for (const pdb_atom &the_atom : arg_residue) {
		if ( the_atom.get_amino_acid().is_proper_amino_acid() ) {
			return the_atom.get_amino_acid();
		}
	}
	// Otherwise, return the amino acid of the first
	return arg_residue.get_atom_cref_of_index( 0 ).get_amino_acid();
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
opt_char cath::file::get_amino_acid_letter(const pdb_residue &arg_residue ///< The pdb_residue to query
                                           ) {
	const amino_acid &the_amino_acid = get_amino_acid( arg_residue );
	return the_amino_acid.is_proper_amino_acid() ? opt_char( the_amino_acid.get_letter() )
	                                             : none;
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
string cath::file::get_amino_acid_code(const pdb_residue &arg_residue ///< The pdb_residue to query
                                       ) {
	return get_amino_acid(arg_residue).get_code();
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
string cath::file::get_amino_acid_name(const pdb_residue &arg_residue ///< The pdb_residue to query
                                       ) {
	return get_amino_acid(arg_residue).get_name();
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::has_atom_of_id_of_residue(const pdb_residue &arg_pdb_residue, ///< The residue from which to extract the atom coordinates
                                           const string      &arg_pdb_atom_id  ///< TODOCUMENT
                                           ) {
	return any_of(
		arg_pdb_residue,
		[&] (const pdb_atom &x) { return x.get_element_type() == arg_pdb_atom_id; }
	);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::has_nitrogen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                               ) {
	return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_NITROGEN );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::has_carbon_alpha_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                                   ) {
	return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON_ALPHA );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::has_carbon_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                             ) {
	return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::has_carbon_beta_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                                  ) {
	return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON_BETA );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::has_oxygen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                             ) {
	return has_atom_of_id_of_residue(arg_pdb_residue, pdb_atom::PDB_ID_OXYGEN);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::is_backbone_complete(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                      ) {
	// If the residue has N, CA and C then add it to new_pdb_residues
	const bool has_n  = has_nitrogen_coord_of_residue     ( arg_pdb_residue );
	const bool has_ca = has_carbon_alpha_coord_of_residue ( arg_pdb_residue );
	const bool has_c  = has_carbon_coord_of_residue       ( arg_pdb_residue );
	return ( has_n && has_ca && has_c );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
coord cath::file::get_atom_of_id_of_residue(const pdb_residue &arg_pdb_residue, ///< The residue from which to extract the atom coordinates
                                            const string      &arg_pdb_atom_id  ///< TODOCUMENT
                                            ) {
	const size_t num_atoms = arg_pdb_residue.get_num_atoms();
	for (size_t atom_ctr = 0; atom_ctr < num_atoms; ++atom_ctr) {
		const pdb_atom the_atom = arg_pdb_residue.get_atom_cref_of_index(atom_ctr);
		if ( the_atom.get_element_type() == arg_pdb_atom_id ) {
			return the_atom.get_coord();
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception(
		"Cannot find atom of type "
		+ arg_pdb_atom_id
		+ " within the "
		+ lexical_cast<string>( num_atoms )
		+ " atom(s) of residue "
		+ lexical_cast<string>( arg_pdb_residue.get_residue_name() )
	));
	return coord::ORIGIN_COORD; // To appease Eclipse's syntax parser
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
coord cath::file::get_nitrogen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                                ) {
	return get_atom_of_id_of_residue(arg_pdb_residue, pdb_atom::PDB_ID_NITROGEN);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
coord cath::file::get_carbon_alpha_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                                    ) {
	return get_atom_of_id_of_residue(arg_pdb_residue, pdb_atom::PDB_ID_CARBON_ALPHA);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
coord cath::file::get_carbon_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                              ) {
	return get_atom_of_id_of_residue(arg_pdb_residue, pdb_atom::PDB_ID_CARBON);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
coord cath::file::get_carbon_beta_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                                   ) {
	return get_atom_of_id_of_residue(arg_pdb_residue, pdb_atom::PDB_ID_CARBON_BETA);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
coord cath::file::get_oxygen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                              ) {
	return get_atom_of_id_of_residue(arg_pdb_residue, pdb_atom::PDB_ID_OXYGEN);
}

/// \brief Fake the location of a CB atom, typically used for glycine residues
///
/// \relates pdb_residue
coord cath::file::fake_carbon_beta_coord_of_residue(const pdb_residue &arg_residue ///< TODOCUMENT
                                                    ) {
	// Prepare a typical CA to CB vector once a residue has been rotated under its ssap_frame
	const double TYPICAL_CA_TO_CB_DISTANCE = 1.527;
	const coord  TYPICAL_CA_TO_CB_UNDER_SSAP_FRAME(
		0,
		- sqrt(2.0/3.0) * TYPICAL_CA_TO_CB_DISTANCE,
		  sqrt(1.0/3.0) * TYPICAL_CA_TO_CB_DISTANCE
	);

	// Grab the CA position and SSAP frame from this residue
	const coord    ca               = get_carbon_alpha_coord_of_residue(arg_residue);
	const rotation ssap_frame       = get_ssap_frame_of_residue(arg_residue);

	//Rotate the TYPICAL_CA_TO_CB_UNDER_SSAP_FRAME by the inverse of this residue's SSAP frame,
	// add this residue's CA position and then return the result
	const coord rotated_ca_to_cb = rotate_copy(transpose_copy(ssap_frame), TYPICAL_CA_TO_CB_UNDER_SSAP_FRAME);
	const coord final_cb         = ca + rotated_ca_to_cb;
	return final_cb;
}

/// \brief Get the carbon beta coord of the residue if it has one or fake it otherwise (ie for glycine)
///
/// \relates pdb_residue
coord cath::file::get_or_predict_carbon_beta_coord_of_residue(const pdb_residue &arg_residue ///< TODOCUMENT
                                                              ) {
	return has_carbon_beta_coord_of_residue(arg_residue) ? get_carbon_beta_coord_of_residue(  arg_residue )
	                                                     : fake_carbon_beta_coord_of_residue( arg_residue );
}

/// \brief Get the SSAP-frame of the residue
///
/// TODOCUMENT NOW
///
/// \relates pdb_residue
rotation cath::file::get_ssap_frame_of_residue(const pdb_residue &arg_residue ///< The residue to query
                                               ) {
	// Calculate the frame of the residue
	const coord n_coord  = get_nitrogen_coord_of_residue(     arg_residue );
	const coord ca_coord = get_carbon_alpha_coord_of_residue( arg_residue );
	const coord c_coord  = get_carbon_coord_of_residue(       arg_residue );
	return construct_residue_frame(n_coord, ca_coord, c_coord);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
ostream & cath::file::write_pdb_file_entry(ostream           &arg_os,         ///< TODOCUMENT
                                           const pdb_residue &arg_residue ///< TODOCUMENT
                                           ) {
	const size_t num_atoms = arg_residue.get_num_atoms();
	for (size_t atom_ctr = 0; atom_ctr < num_atoms; ++atom_ctr) {
		const pdb_atom &atom = arg_residue.get_atom_cref_of_index(atom_ctr);
		write_pdb_file_entry(
			arg_os,
			arg_residue.get_chain_label(),
			arg_residue.get_residue_name(),
			atom
		);
		arg_os << "\n";
	}
	return arg_os;
}


/// \brief Get the psi angle of this residue and the phi angle of the next residue, both in positive radians within (0, 2*pi]
///
/// \relates pdb_residue
doub_angle_doub_angle_pair cath::file::get_psi_of_this_and_phi_of_next(const pdb_residue &arg_this_residue, ///< This residue
                                                                       const pdb_residue &arg_next_residue  ///< The next residue
                                                                       ) {
	// Grab the positions of this residue
	const coord this_n ( get_nitrogen_coord_of_residue    ( arg_this_residue ) );
	const coord this_ca( get_carbon_alpha_coord_of_residue( arg_this_residue ) );
	const coord this_c ( get_carbon_coord_of_residue      ( arg_this_residue ) );

	// ...and the next one
	const coord next_n ( get_nitrogen_coord_of_residue    ( arg_next_residue ) );
	const coord next_ca( get_carbon_alpha_coord_of_residue( arg_next_residue ) );
	const coord next_c ( get_carbon_coord_of_residue      ( arg_next_residue ) );

	// Compute the angles
	const auto psi_of_this = dihedral_angle_between_four_points( this_n, this_ca, this_c,  next_n );
	const auto phi_of_next = dihedral_angle_between_four_points( this_c, next_n,  next_ca, next_c );

	// Ensure angles are positive (ie between 0 and 2*pi, using 2*pi for endpoint values)
	const auto shifted_psi_of_this = shift_copy( psi_of_this, zero_angle<double>(), angle_endpoint_loc::USE_UPPER );
	const auto shifted_phi_of_next = shift_copy( phi_of_next, zero_angle<double>(), angle_endpoint_loc::USE_UPPER );

	// Return the computed angles
	return make_pair( shifted_psi_of_this, shifted_phi_of_next );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
residue cath::file::build_residue_of_pdb_residue(const pdb_residue &arg_residue, ///< TODOCUMENT
                                                 const doub_angle  &arg_phi,     ///< TODOCUMENT
                                                 const doub_angle  &arg_psi,     ///< TODOCUMENT
                                                 const size_t      &arg_access   ///< TODOCUMENT
                                                 ) {
	// Calculate the CA, CB and ssap_frame of the residue
	const coord    ca_coord   = get_carbon_alpha_coord_of_residue          ( arg_residue );
	const coord    cb_coord   = get_or_predict_carbon_beta_coord_of_residue( arg_residue );
	const rotation ssap_frame = get_ssap_frame_of_residue                  ( arg_residue );

	// Build a residue object and return it
	return residue(
		arg_residue.get_residue_name(),
		get_amino_acid( arg_residue ),
		ca_coord,
		cb_coord,
		0,
		sec_struc_type::COIL,
		ssap_frame,
		arg_phi,
		arg_psi,
		arg_access
	);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
ostream & cath::file::operator<<(ostream               &arg_os,          ///< TODOCUMENT
                                 const pdb_residue_vec &arg_pdb_residues ///< TODOCUMENT
                                 ) {
	str_vec pdb_residue_strings;
	pdb_residue_strings.reserve( arg_pdb_residues.size() );
	for (const pdb_residue &the_residue : arg_pdb_residues) {
		pdb_residue_strings.push_back( lexical_cast<string>( the_residue ) );
	}
	arg_os << join( pdb_residue_strings, ", ");
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
ostream & cath::file::operator<<(ostream           &arg_os, ///< TODOCUMENT
                                 const pdb_residue &arg_pdb_residue ///< TODOCUMENT
                                 ) {
	arg_os << "Residue[";
	arg_os << arg_pdb_residue.get_chain_label();
	arg_os << ", ";
	arg_os << arg_pdb_residue.get_residue_name();
	arg_os << ", ";
	arg_os << get_carbon_alpha_coord_of_residue(arg_pdb_residue);
// 	arg_os << ", ";
// 	arg_os << arg_pdb_residue.get_amino_acid();
	arg_os << "]";
	return arg_os;
}

