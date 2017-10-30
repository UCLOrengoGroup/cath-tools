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

#include "pdb_residue.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/boost_addenda/range/max_proj_element.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "structure/protein/residue.hpp"

#include <cmath>

using namespace boost::math::constants;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::adaptors::reversed;
using boost::adaptors::filtered;
using boost::algorithm::any_of;
using boost::algorithm::join;
using boost::lexical_cast;
using boost::none;

/// \brief TODOCUMENT
///
/// \relates pdb_residue
char_opt cath::file::get_letter_if_amino_acid(const pdb_residue &arg_residue ///< The pdb_residue to query
                                              ) {
	return arg_residue.get_amino_acid().get_letter_if_amino_acid();
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
char cath::file::get_amino_acid_letter_tolerantly(const pdb_residue &arg_residue ///< The pdb_residue to query
                                                  ) {
	return arg_residue.get_amino_acid().get_letter_tolerantly();
}

/// \brief Get the three-letter-code char_3_arr for the amino acid in the specified pdb_residue
///
/// \relates pdb_residue
char_3_arr cath::file::get_amino_acid_code(const pdb_residue &arg_residue ///< The pdb_residue to query
                                           ) {
	return arg_residue.get_amino_acid().get_code();
}

/// \brief Get the three-letter-code string for the amino acid in the specified pdb_residue
///
/// \relates pdb_residue
string cath::file::get_amino_acid_code_string(const pdb_residue &arg_residue ///< The pdb_residue to query
                                              ) {
	return get_code_string( arg_residue.get_amino_acid() );
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
string cath::file::get_amino_acid_name(const pdb_residue &arg_residue ///< The pdb_residue to query
                                       ) {
	return arg_residue.get_amino_acid().get_name();
}

/// \brief Whether the specified pdb_residue contains any pdb_atoms that aren't proper amino-acids
///        (eg from HETATM records)
///
/// \relates pdb_residue
residue_makeup cath::file::contains_non_proper_amino_acids(const pdb_residue &arg_residue ///< The pdb_residue to query
                                                           ) {
	return
		any_of(
			arg_residue,
			[] (const pdb_atom &x) {
				return ! is_proper_amino_acid( x.get_amino_acid() );
			}
		)
		? residue_makeup::SOME_NON_PROPER_AMINO_ACIDS
		: residue_makeup::ALL_PROPER_AMINO_ACIDS;
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
bool cath::file::is_backbone_complete(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
                                      ) {
	// If the residue has N, CA and C then add it to new_pdb_residues
	const bool has_n  = arg_pdb_residue.has_nitrogen();
	const bool has_ca = arg_pdb_residue.has_carbon_alpha();
	const bool has_c  = arg_pdb_residue.has_carbon();
	return ( has_n && has_ca && has_c );
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
	const coord    ca               = get_carbon_alpha_coord   ( arg_residue );
	const rotation ssap_frame       = get_ssap_frame_of_residue( arg_residue );

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
	return arg_residue.has_carbon_beta() ? get_carbon_beta_coord            ( arg_residue )
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
	const coord n_coord  = get_nitrogen_coord    ( arg_residue );
	const coord ca_coord = get_carbon_alpha_coord( arg_residue );
	const coord c_coord  = get_carbon_coord      ( arg_residue );
	return construct_residue_frame(n_coord, ca_coord, c_coord);
}

/// \brief Return the maximum distance between the specified pdb_residue's atoms and the specified coord
///        or 0.0 if the pdb_residue is empty
///
/// \relates pdb_residue
double cath::file::max_dist_from_coord(const pdb_residue &arg_pdb_residue, ///< The pdb_residue to query
                                       const coord       &arg_coord        ///< The coord to which the distances should be measured
                                       ) {
	if ( arg_pdb_residue.empty() ) {
		return 0.0;
	}
	return max_proj(
		arg_pdb_residue,
		less<>{},
		[&] (const pdb_atom &atom) {
			return distance_between_points( arg_coord, atom.get_coord() );
		}
	);
}

/// \brief TODOCUMENT
///
/// \relates pdb_residue
ostream & cath::file::write_pdb_file_entry(ostream           &arg_os,     ///< TODOCUMENT
                                           const pdb_residue &arg_residue ///< TODOCUMENT
                                           ) {
	const size_t num_atoms = arg_residue.get_num_atoms();
	for (const size_t &atom_ctr : indices( num_atoms ) ) {
		const pdb_atom &atom = arg_residue.get_atom_cref_of_index(atom_ctr);
		write_pdb_file_entry(
			arg_os,
			arg_residue.get_residue_id(),
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
	const coord this_n ( get_nitrogen_coord    ( arg_this_residue ) );
	const coord this_ca( get_carbon_alpha_coord( arg_this_residue ) );
	const coord this_c ( get_carbon_coord      ( arg_this_residue ) );

	// ...and the next one
	const coord next_n ( get_nitrogen_coord    ( arg_next_residue ) );
	const coord next_ca( get_carbon_alpha_coord( arg_next_residue ) );
	const coord next_c ( get_carbon_coord      ( arg_next_residue ) );

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
	const coord    ca_coord   = get_carbon_alpha_coord( arg_residue );
	const coord    cb_coord   = get_or_predict_carbon_beta_coord_of_residue( arg_residue );
	const rotation ssap_frame = get_ssap_frame_of_residue                  ( arg_residue );

	// Build a residue object and return it
	return residue(
		arg_residue.get_residue_id(),
		arg_residue.get_amino_acid(),
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

/// \brief Whether DSSP will skip this residue
///
/// This attempts to match the DSSP criteria for skipping a residue:
///  * has at least one N, CA, C and O atom that meets
///    * alt_locn_is_dssp_accepted() if the first atom meets alt_locn_is_dssp_accepted()
///    * or `( x.get_alt_locn() == ' ' )` otherwise
///
/// I've found an issue in DSSP that, for example, causes it to ignore a valid residue with
/// mixed altlocn of 'A'/'B' ('A' first) to be skipped iff it appears after a residue
/// with all rejected altlocn (eg 'B'). I'm in contact with the DSSP people about this and
/// may be sending them a pull-request with a fix soon. I'll see what happens with that
/// before spending any time trying to replicate that behaviour here.
///
/// \relates pdb_residue
bool cath::file::dssp_will_skip_residue(const pdb_residue &arg_pdb_residue ///< The pdb_residue to test
                                        ) {
	// Determine whether there is a first atom with non-standard alt-locn
	const bool has_nonstd_1st_altloc = ! arg_pdb_residue.empty() && ! alt_locn_is_dssp_accepted( front( arg_pdb_residue ) );

	// Create a closure for checking whether a PDB atom's altlocn is acceptable given has_nonstd_1st_altloc
	const auto is_accepted_locn = [&] (const pdb_atom &x) {
		return has_nonstd_1st_altloc
			? ( x.get_alt_locn() == ' ' )
			: alt_locn_is_dssp_accepted( x );
	};

	// Loop over the atoms, detecting whether any of them is a valid N/CA/C/O
	//
	// Could probably tidy this up with a bit of metaprogramming
	bool has_n  = false;
	bool has_ca = false;
	bool has_c  = false;
	bool has_o  = false;
	for (const pdb_atom &the_atom : arg_pdb_residue | filtered( is_accepted_locn ) ) {
		const coarse_element_type coarse_element = get_coarse_element_type( the_atom );
		if ( ! has_n  && coarse_element == coarse_element_type::NITROGEN     ) {
			has_n  = true;
			continue;
		}
		if ( ! has_ca && coarse_element == coarse_element_type::CARBON_ALPHA ) {
			has_ca = true;
			continue;
		}
		if ( ! has_c  && coarse_element == coarse_element_type::CARBON       ) {
			has_c  = true;
			continue;
		}
		if ( ! has_o  && coarse_element == coarse_element_type::OXYGEN       ) {
			has_o  = true;
			continue;
		}
	}

	// Return whether no acceptable ATOM was found for any of N/CA/C/O
	return ! ( has_n && has_ca && has_c && has_o );
}

/// \brief Get all coordinates of all the pdb_atoms in the specified pdb_residues
///
/// \relates pdb_residue
coord_vec cath::file::get_all_coords(const pdb_residue_vec &arg_pdb_residues ///< The pdb_residues to query
                                     ) {
	coord_vec results;
	for (const pdb_residue &res : arg_pdb_residues) {
		for (const pdb_atom &atom : res) {
			results.push_back( atom.get_coord() );
		}
	}
	return results;
}

/// \brief Get all coordinates of all the pdb_atoms in the specified pdb_residues
///
/// \relates pdb_residue
coord_coord_linkage_pair_vec cath::file::get_all_coords_with_linkage(const pdb_residue_vec &arg_pdb_residues ///< The pdb_residues to query
                                                                     ) {
	coord_coord_linkage_pair_vec results;
	for (const pdb_residue &res : arg_pdb_residues) {
		for (const pdb_atom &atom : res) {
			results.emplace_back(
				atom.get_coord(),
				is_water( atom )
					? coord_linkage::ADD_ONLY
					: coord_linkage::ADD_AND_LINK
			);
		}
	}
	return results;
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
	arg_os << arg_pdb_residue.get_residue_id();
	if ( arg_pdb_residue.has_carbon_alpha() ) {
		arg_os << ", ";
		arg_os << get_carbon_alpha_coord( arg_pdb_residue );
	}
// 	arg_os << ", ";
// 	arg_os << arg_pdb_residue.get_amino_acid();
	arg_os << "]";
	return arg_os;
}

