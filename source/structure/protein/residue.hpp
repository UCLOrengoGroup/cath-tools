/// \file
/// \brief The residue class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_RESIDUE_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_RESIDUE_H

#include <boost/operators.hpp>

#include "file/pdb/dssp_skip_policy.hpp"
#include "structure/geometry/angle.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/rotation.hpp"
#include "structure/protein/amino_acid.hpp"
#include "structure/protein/sec_struc_type.hpp"
#include "structure/residue_id.hpp"

#include <cstddef>
#include <iosfwd>

namespace cath {

	/// \brief TODOCUMENT
	///
	/// \todo This currently hard-codes the member+getter for two atoms coordinates: CA and CB
	///       These two should remain and remain mandatory but they should be complemented by
	///       N, C and O but these should be optional (with has_atom_coord method).
	///       This looks like a strong candidate for a template-based solution...
	///       Perhaps something like a Boost Fusion container that stores the coords based on enum values
	///       and template method accessors like `the_residue.has_atom_coord<ATOM_O>()` and
	///       `the_residue.get_atom_coord<ATOM_O>()`, and some static method of detecting which values
	///       are mandatory and which optional.
	///       The code that does this should be encapsulated in an implementation class because this class
	///       may be independently useful, eg for potentially applying common_atom_selection_policy to
	///       altering the distances used in SSAP scores.
	class residue final : private boost::equality_comparable<residue> {
	private:
		/// \brief TODOCUMENT
		residue_id     the_residue_id;

		/// \brief Single character representing the amino acid ('0' or an upper case letter)
		amino_acid     the_amino_acid;

		/// \brief Coordinates of the carbon alpha atom
		geom::coord    carbon_alpha_coord;

		/// \brief Coordinates of the carbon beta atom
		geom::coord    carbon_beta_coord;

		/// \brief TODOCUMENT
		size_t         sec_struc_number;

		/// \brief TODOCUMENT
		sec_struc_type the_sec_struc_type; ///< \todo Check: is it a problem that the sec files have values of H/S but the source compares this to H/E.
		                                   ///<       At least sec_struc.get_type() is used to check that one sec_struc's type matches another.

		/// \brief TODOCUMENT
		geom::rotation   frame;

		/// \brief The phi angle of this residue, in degrees within (0, 360.0] (possibly [0, 360.0] ? )
		geom::doub_angle phi_angle;

		/// \brief The psi angle of this residue, in degrees within (0, 360.0] (possibly [0, 360.0] ? )
		geom::doub_angle psi_angle;

		/// \brief TODOCUMENT
		size_t           access;

	private:
		static void check_phi_psi_angle(const geom::doub_angle &);

	public:
		residue(const residue_id &,
		        const amino_acid &,
		        const geom::coord &,
		        const geom::coord &,
		        const size_t &,
		        const sec_struc_type &,
		        const geom::rotation &,
		        const geom::doub_angle &,
		        const geom::doub_angle &,
		        const size_t &);

		void initialise();

		void set_amino_acid(const amino_acid &);
		void set_residue_sec_struc_number(const size_t &);
		void set_sec_struc_type(const sec_struc_type &);

		const residue_id & get_pdb_residue_id() const;
		amino_acid get_amino_acid() const;
		inline const geom::coord & get_carbon_alpha_coord() const;
		inline const geom::coord & get_carbon_beta_coord() const;

		size_t get_sec_struc_number() const;
		sec_struc_type get_sec_struc_type() const;

		inline const geom::rotation & get_frame() const;
		const geom::doub_angle & get_phi_angle() const;
		const geom::doub_angle & get_psi_angle() const;
		size_t get_access() const;

		static geom::doub_angle DEFAULT_PHI_PSI();

		/// \todo Sort out indexing and then get rid of this because a blank residue doesn't really make any sense
		static const residue NULL_RESIDUE;
	};

	/// \brief TODOCUMENT
	const geom::coord & residue::get_carbon_alpha_coord() const {
		return carbon_alpha_coord;
	}

	/// \brief TODOCUMENT
	const geom::coord & residue::get_carbon_beta_coord() const {
		return carbon_beta_coord;
	}

	/// \brief TODOCUMENT
	const geom::rotation & residue::get_frame() const {
		return frame;
	}

	bool operator==(const residue &,
	                const residue &);

	const chain_label & get_chain_label(const residue &);
	const residue_name & get_pdb_residue_name(const residue &);

	const int & pdb_number(const residue &);
	int pdb_number_or_value_if_null(const residue &,
	                                const int &);
	bool has_pdb_insert(const residue &);
	bool has_pdb_insert_or_value_if_null(const residue &,
	                                     const bool &);
	const char & pdb_insert(const residue &);
	char pdb_insert_or_value_if_null(const residue &,
	                                 const char &);
	char pdb_insert_or_value_if_null_or_absent(const residue &,
	                                           const char &);

	bool residue_matches_residue_id(const residue &,
	                                const residue_id &);

	char get_amino_acid_letter(const residue &);

	std::string get_amino_acid_code(const residue &);

	std::string get_amino_acid_name(const residue &);

	std::string ssap_legacy_alignment_left_side_string(const cath::residue &);

	std::string ssap_legacy_alignment_left_side_gap_string();

	std::string ssap_legacy_alignment_right_side_string(const cath::residue &);

	std::string ssap_legacy_alignment_right_side_gap_string();

	std::string get_pdb_residue_id_string(const residue &);

	int get_accessi_of_residue(const residue &);

	std::ostream & operator<<(std::ostream &,
	                          const residue &);

	geom::rotation construct_residue_frame(const geom::coord &,
	                                       const geom::coord &,
	                                       const geom::coord &);

	residue combine_residues_from_dssp_and_pdb(const residue &,
	                                           const residue &,
	                                           const file::dssp_skip_angle_skipping & = file::dssp_skip_angle_skipping::BREAK_ANGLES);

	bool is_null_residue(const residue &);

	void wipe_secondary_structure(residue &);
} // namespace cath
#endif
