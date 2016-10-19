/// \file
/// \brief The pdb_residue class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_PDB_PDB_RESIDUE_H
#define _CATH_TOOLS_SOURCE_FILE_PDB_PDB_RESIDUE_H

#include <boost/operators.hpp>

#include "common/type_aliases.h"
#include "file/file_type_aliases.h"
#include "structure/chain_label.h"
#include "structure/residue_name.h"
#include "structure/structure_type_aliases.h"

#include <vector>

namespace cath { class amino_acid; }
namespace cath { namespace file { class pdb_atom; } }
namespace cath { namespace geom { class coord; } }
namespace cath { namespace geom { class rotation; } }
namespace cath { class residue; }

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		///
		/// The amino acid cannot be stored here, rather than in the pdb_atom entries because
		/// some residues (eg residues 22 and 25 in 1cbnA) have different amino acid values in
		/// the same residue. Storing that the residue level would mean that information could
		/// not be reproduced when outputting the pdb.
		class pdb_residue final : boost::additive<pdb_residue, geom::coord> {
		private:
			/// \brief TODOCUMENT
			chain_label  the_chain_label;

			/// \brief TODOCUMENT
			residue_name the_residue_name;

			/// \brief TODOCUMENT
			pdb_atom_vec atoms;

		public:
			/// \brief TODOCUMENT
			using const_iterator = pdb_atom_vec::const_iterator;

			pdb_residue(const chain_label &,
			            const residue_name &,
			            const pdb_atom_vec &);

			pdb_residue(const chain_label &,
			            const residue_name &,
			            pdb_atom_vec &&);

			chain_label get_chain_label() const;
			residue_name get_residue_name() const;
//			coord get_carbon_alpha_coord() const;
// 			std::string get_amino_acid() const;
			size_t get_num_atoms() const;
			const pdb_atom & get_atom_cref_of_index(const size_t &) const;

			void set_chain_label(const chain_label &);

			void rotate(const geom::rotation &);
			void operator+=(const geom::coord &);
			void operator-=(const geom::coord &);

			const_iterator begin() const;
			const_iterator end() const;
		};

		const amino_acid & get_amino_acid(const pdb_residue &);
		char_opt get_amino_acid_letter(const pdb_residue &);
		std::string get_amino_acid_code(const pdb_residue &);
		std::string get_amino_acid_name(const pdb_residue &);

		bool has_atom_of_id_of_residue(const pdb_residue &,
		                               const std::string &);
		bool has_nitrogen_coord_of_residue(const pdb_residue &);
		bool has_carbon_alpha_coord_of_residue(const pdb_residue &);
		bool has_carbon_coord_of_residue(const pdb_residue &);
		bool has_carbon_beta_coord_of_residue(const pdb_residue &);
		bool has_oxygen_coord_of_residue(const pdb_residue &);
		bool is_backbone_complete(const pdb_residue &);

		geom::coord get_atom_of_id_of_residue(const pdb_residue &,
		                                      const std::string &);
		geom::coord get_nitrogen_coord_of_residue(const pdb_residue &);
		geom::coord get_carbon_alpha_coord_of_residue(const pdb_residue &);
		geom::coord get_carbon_coord_of_residue(const pdb_residue &);
		geom::coord get_carbon_beta_coord_of_residue(const pdb_residue &);
		geom::coord get_oxygen_coord_of_residue(const pdb_residue &);

		geom::coord fake_carbon_beta_coord_of_residue(const pdb_residue &);
		geom::coord get_or_predict_carbon_beta_coord_of_residue(const pdb_residue &);
		geom::rotation get_ssap_frame_of_residue(const pdb_residue &);

		std::ostream & write_pdb_file_entry(std::ostream &,
		                                    const pdb_residue &);

		geom::angle_angle_pair<double> get_psi_of_this_and_phi_of_next(const pdb_residue &,
		                                                               const pdb_residue &);

		residue build_residue_of_pdb_residue(const pdb_residue &,
		                                     const geom::doub_angle &,
		                                     const geom::doub_angle &,
		                                     const size_t &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_residue_vec &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_residue &);
	} // namespace file
} // namespace cath

#endif
