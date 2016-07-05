/// \file
/// \brief The pdb class header

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

#ifndef PDB_H_INCLUDED
#define PDB_H_INCLUDED

#include <boost/operators.hpp>

#include "file/file_type_aliases.h"
#include "file/pdb/pdb_base.h"
#include "structure/structure_type_aliases.h"

#include <iostream>
#include <vector>

namespace cath { namespace chop { class domain; } }
namespace cath { namespace file { class pdb_residue; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { class protein; }

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		///
		/// \todo Change to do reading and writing via streams
		class pdb final : public pdb_base, private boost::additive<pdb, geom::coord> {
		private:
			/// \brief TODOCUMENT
			pdb_residue_vec pdb_residues;

			virtual void do_read_file(const boost::filesystem::path &) override final;
			virtual void do_append_to_file(const boost::filesystem::path &) const override final;
			virtual void do_set_chain_label(const chain_label &) override final;
			virtual residue_name_vec do_get_residue_names_of_first_chain__backbone_unchecked() const override final;
			virtual geom::coord do_get_residue_ca_coord_of_index__backbone_unchecked(const size_t &) const override final;
			virtual size_t do_get_num_atoms() const override final;

			virtual void do_rotate(const geom::rotation &) override final;
			virtual void do_add(const geom::coord &) override final;
			virtual void do_subtract(const geom::coord &) override final;

		public:
			virtual ~pdb() noexcept = default;

			size_t get_num_backbone_complete_residues() const;
			size_t get_index_of_backbone_complete_index(const size_t &) const;

			size_t get_num_residues() const;
			const pdb_residue & get_residue_cref_of_index__backbone_unchecked(const size_t &) const;
			void set_residues(const pdb_residue_vec &);

			const pdb_residue & get_residue_cref_of_backbone_complete_index(const size_t &) const;
			geom::coord get_residue_ca_coord_of_backbone_complete_index(const size_t &) const;
			residue_name_vec get_backbone_complete_residue_names_of_first_chain(const bool & = true) const;

			/// \brief TODOCUMENT
			using const_iterator = pdb_residue_vec::const_iterator;

			const_iterator begin() const;
			const_iterator end() const;
		};

		pdb read_pdb_file(const boost::filesystem::path &);

		pdb read_domain_from_pdb_file(const boost::filesystem::path &,
		                              const chop::domain &);

		pdb read_pdb_file(std::istream &);
		std::istream & read_pdb_file(std::istream &,
		                             pdb &);
		pdb_list read_end_separated_pdb_files(std::istream &);
		std::ostream & write_pdb_file(std::ostream &,
		                              const pdb &);

		amino_acid_vec get_amino_acid_list(const pdb &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb &);

//		const pdb_residue & get_residue_ref_of_index__offset_1(const pdb &,
//		                                                       const size_t &);

		geom::doub_angle_doub_angle_pair_vec get_phi_and_psi_angles(const pdb &);

		pdb backbone_complete_subset_of_pdb(const pdb &,
		                                    std::ostream & = std::cerr);

		protein build_protein_of_pdb(const pdb &,
		                             std::ostream & = std::cerr);

		protein build_protein_of_pdb_and_name(const pdb &,
		                                      const std::string &,
		                                      std::ostream & = std::cerr);
	}
}

#endif
