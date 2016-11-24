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

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/operators.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "file/file_type_aliases.h"
#include "file/pdb/pdb_atom.h"
#include "structure/chain_label.h"
#include "structure/residue_name.h"
#include "structure/structure_type_aliases.h"

#include <vector>

using namespace cath::common::literals;

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
			bool empty() const;
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

		const geom::coord & get_atom_of_id_of_residue(const pdb_residue &,
		                                              const std::string &);
		const geom::coord & get_nitrogen_coord_of_residue(const pdb_residue &);
		const geom::coord & get_carbon_alpha_coord_of_residue(const pdb_residue &);
		const geom::coord & get_carbon_coord_of_residue(const pdb_residue &);
		const geom::coord & get_carbon_beta_coord_of_residue(const pdb_residue &);
		const geom::coord & get_oxygen_coord_of_residue(const pdb_residue &);

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

		bool dssp_might_skip_residue(const pdb_residue &);

		bool dssp_will_skip_residue(const pdb_residue &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_residue_vec &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_residue &);

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline bool has_atom_of_id_of_residue(const pdb_residue &arg_pdb_residue, ///< The residue from which to extract the atom coordinates
		                                      const std::string &arg_pdb_atom_id  ///< TODOCUMENT
		                                      ) {
			return boost::algorithm::any_of(
				arg_pdb_residue,
				[&] (const pdb_atom &x) { return x.get_element_type() == arg_pdb_atom_id; }
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline bool has_nitrogen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                          ) {
			return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_NITROGEN );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline bool has_carbon_alpha_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                              ) {
			return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON_ALPHA );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline bool has_carbon_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                        ) {
			return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline bool has_carbon_beta_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                             ) {
			return has_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON_BETA );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline bool has_oxygen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                        ) {
			return has_atom_of_id_of_residue(arg_pdb_residue, pdb_atom::PDB_ID_OXYGEN);
		}

		/// \brief Get an ato of the specified type from the specified residue
		///
		/// This is organised so that, in the cases where there are multiple atoms of the same type,
		/// it will replicate DSSP's implicit policy :
		///
		///  * Prefer an atom that meets alt_locn_is_dssp_accepted() (ie altlocn of ' ' or 'A'), then
		///  * Prefer an atom nearer the back
		///
		/// \relates pdb_residue
		inline const geom::coord & get_atom_of_id_of_residue(const pdb_residue &arg_pdb_residue, ///< The residue from which to extract the atom coordinates
		                                                     const std::string &arg_pdb_atom_id  ///< The type of atom to retrieve (eg pdb_atom::PDB_ID_CARBON_ALPHA)
		                                                     ) {
			size_opt best_seen_atom_ctr;
			for (const size_t &atom_ctr : boost::irange( 0_z, arg_pdb_residue.get_num_atoms() ) | boost::adaptors::reversed) {
				const pdb_atom &the_atom = arg_pdb_residue.get_atom_cref_of_index( atom_ctr );
				if ( the_atom.get_element_type() == arg_pdb_atom_id ) {
					if ( alt_locn_is_dssp_accepted( the_atom ) ) {
						return the_atom.get_coord();
					}
					if ( ! best_seen_atom_ctr ) {
						best_seen_atom_ctr = atom_ctr;
					}
				}
			}
			if ( best_seen_atom_ctr ) {
				return arg_pdb_residue.get_atom_cref_of_index( *best_seen_atom_ctr ).get_coord();
			}
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
				"Cannot find atom of type "
				+ arg_pdb_atom_id
				+ " within the "
				+ boost::lexical_cast<std::string>( arg_pdb_residue.get_num_atoms() )
				+ " atom(s) of residue "
				+ boost::lexical_cast<std::string>( arg_pdb_residue.get_residue_name() )
			));
			return geom::coord::ORIGIN_COORD; // To appease Eclipse's syntax parser
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline const geom::coord & get_nitrogen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                                         ) {
			return get_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_NITROGEN );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline const geom::coord & get_carbon_alpha_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                                             ) {
			return get_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON_ALPHA );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline const geom::coord & get_carbon_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                                       ) {
			return get_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline const geom::coord & get_carbon_beta_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                                            ) {
			return get_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_CARBON_BETA );
		}

		/// \brief TODOCUMENT
		///
		/// \relates pdb_residue
		inline const geom::coord & get_oxygen_coord_of_residue(const pdb_residue &arg_pdb_residue ///< The residue from which to extract the atom coordinates
		                                                       ) {
			return get_atom_of_id_of_residue( arg_pdb_residue, pdb_atom::PDB_ID_OXYGEN );
		}

	} // namespace file
} // namespace cath

#endif
