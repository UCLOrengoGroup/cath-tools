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

#ifndef _CATH_TOOLS_SOURCE_FILE_PDB_PDB_H
#define _CATH_TOOLS_SOURCE_FILE_PDB_PDB_H

#include <boost/operators.hpp>
#include <boost/optional.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "file/file_type_aliases.hpp"
#include "file/pdb/pdb_base.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "structure/structure_type_aliases.hpp"

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
			size_t get_num_backbone_complete_residues() const;
			size_t get_index_of_backbone_complete_index(const size_t &) const;

			size_t get_num_residues() const;
			const pdb_residue & get_residue_cref_of_index__backbone_unchecked(const size_t &) const;
			void set_residues(const pdb_residue_vec &);
			void set_residues(pdb_residue_vec &&);

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

		geom::doub_angle_doub_angle_pair_vec get_phi_and_psi_angles(const pdb &);

		pdb backbone_complete_subset_of_pdb(const pdb &,
		                                    const ostream_ref_opt & = boost::none,
		                                    const bool & = false);

		protein build_protein_of_pdb(const pdb &,
		                             const ostream_ref_opt & = boost::none);

		protein build_protein_of_pdb_and_name(const pdb &,
		                                      const std::string &,
		                                      const ostream_ref_opt & = boost::none);

		size_set get_protein_res_indices_that_dssp_might_skip(const pdb &,
		                                                      const ostream_ref_opt & = boost::none);

		/// \brief Get the number of residues held in this pdb
		inline size_t pdb::get_num_residues() const {
			return pdb_residues.size();
		}

		/// \brief Get the residue of the specified index
		///        (with no checking for which residues are backbone-complete)
		inline const pdb_residue & pdb::get_residue_cref_of_index__backbone_unchecked(const size_t &arg_index ///< The index of the residue to retun
		                                                                              ) const {
#ifndef NDEBUG
			if ( arg_index >= get_num_residues() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Unable to get_residue_ca_coord_of_index__backbone_unchecked() for index >= number of residues"
				));
			}
#endif
			return pdb_residues[ arg_index ];
		}

		/// \brief Standard const begin method for the range of residues
		inline auto pdb::begin() const -> const_iterator {
			return common::cbegin( pdb_residues );
		}

		/// \brief Standard const end method for the range of residues
		inline auto pdb::end() const -> const_iterator {
			return common::cend( pdb_residues );
		}

	} // namespace file
} // namespace cath

#endif
