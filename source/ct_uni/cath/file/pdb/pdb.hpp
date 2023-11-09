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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_HPP

#include <filesystem>
#include <iostream>
#include <optional>
#include <vector>

#include <boost/operators.hpp>

#include "cath/chopping/region/region.hpp"
#include "cath/chopping/region/regions_limiter.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/config.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/file/pdb/dssp_skip_policy.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/file/pdb/pdb_write_mode.hpp"
#include "cath/structure/structure_type_aliases.hpp"

// clang-format off
namespace cath { class protein; }
namespace cath::file { class pdb_list; }
namespace cath::file { class pdb_residue; }
namespace cath::file { struct protein_info; }
// clang-format on

namespace cath::file {

	/// \brief TODOCUMENT
	///
	/// \todo Change to do reading and writing via streams
	class pdb final : private boost::additive<pdb, geom::coord> {
	private:
		/// \brief TODOCUMENT
		pdb_residue_vec pdb_residues;

		/// \brief The residues that appeared after a TER record in their respective chains
		pdb_residue_vec post_ter_residues;

	public:
		void read_file(const ::std::filesystem::path &);
		void append_to_file(const ::std::filesystem::path &) const;
		pdb & set_chain_label(const chain_label &);
		[[nodiscard]] residue_id_vec get_residue_ids_of_first_chain__backbone_unchecked() const;
		[[nodiscard]] geom::coord    get_residue_ca_coord_of_index__backbone_unchecked( const size_t & ) const;
		[[nodiscard]] size_t         get_num_atoms() const;

		pdb & rotate(const geom::rotation &);
		pdb & operator+=(const geom::coord &);
		pdb & operator-=(const geom::coord &);

		[[nodiscard]] bool               empty() const;
		[[nodiscard]] size_t             get_num_residues() const;
		[[nodiscard]] const pdb_residue &get_residue_of_index__backbone_unchecked( const size_t & ) const;
		pdb & set_residues(pdb_residue_vec);
		pdb & set_post_ter_residues(pdb_residue_vec);

		[[nodiscard]] const pdb_residue_vec &get_post_ter_residues() const;

		/// \brief TODOCUMENT
		using const_iterator = pdb_residue_vec::const_iterator;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	backbone_complete_indices get_backbone_complete_indices(const pdb &);

	size_t get_num_backbone_complete_residues(const pdb &);

	size_t get_index_of_backbone_complete_index(const pdb &,
	                                            const size_t &);

	const pdb_residue & get_residue_of_backbone_complete_index(const pdb &,
	                                                           const size_t &);

	const pdb_residue & get_residue_of_backbone_complete_index(const pdb &,
	                                                           const backbone_complete_indices &,
	                                                           const size_t &);

	geom::coord get_residue_ca_coord_of_backbone_complete_index(const pdb &,
	                                                            const size_t &);

	geom::coord get_residue_ca_coord_of_backbone_complete_index(const pdb &,
	                                                            const backbone_complete_indices &,
	                                                            const size_t &);

	size_t get_num_region_limited_backbone_complete_residues(const pdb &,
	                                                         const chop::region_vec_opt & = ::std::nullopt);

	size_t get_index_of_region_limited_backbone_complete_index(const pdb &,
	                                                           const size_t &,
	                                                           const chop::region_vec_opt & = ::std::nullopt);

	const pdb_residue & get_residue_of_region_limited_backbone_complete_index(const pdb &,
	                                                                          const size_t &,
	                                                                          const chop::region_vec_opt & = ::std::nullopt);
	geom::coord get_residue_ca_coord_of_region_limited_backbone_complete_index(const pdb &,
	                                                                           const size_t &,
	                                                                           const chop::region_vec_opt & = ::std::nullopt);


	residue_id_vec get_backbone_complete_residue_ids_of_first_chain(const pdb &,
	                                                                const bool & = true);
	residue_id_vec get_backbone_complete_residue_ids(const pdb &);

	pdb read_pdb_file(const ::std::filesystem::path &);

	pdb read_pdb_file(std::istream &);
	std::istream & read_pdb_file(std::istream &,
	                             pdb &);
	pdb read_pdb(const std::string &);
	pdb_list read_end_separated_pdb_files(std::istream &);

	std::string to_pdb_file_string(const pdb &,
	                               const chop::region_vec_opt & = ::std::nullopt,
	                               const pdb_write_mode & = pdb_write_mode::ONLY_OR_LAST_PDB);
	std::ostream & write_pdb_file(std::ostream &,
	                              const pdb &,
	                              const chop::region_vec_opt & = ::std::nullopt,
	                              const pdb_write_mode & = pdb_write_mode::ONLY_OR_LAST_PDB);
	void write_pdb_file(const ::std::filesystem::path &,
	                    const pdb &,
	                    const chop::region_vec_opt & = ::std::nullopt,
	                    const pdb_write_mode & = pdb_write_mode::ONLY_OR_LAST_PDB);

	std::string pdb_file_to_string(const pdb &,
	                               const chop::region_vec_opt & = ::std::nullopt,
	                               const pdb_write_mode & = pdb_write_mode::ONLY_OR_LAST_PDB);

	amino_acid_vec get_amino_acid_list(const pdb &);

	size_vec indices_of_residues_following_chain_breaks(const pdb &);

	bool has_multiple_chain_labels(const pdb &);

	std::ostream & operator<<(std::ostream &,
	                          const pdb &);

	geom::doub_angle_doub_angle_pair_vec get_phi_and_psi_angles(const pdb &,
	                                                            // const size_vec &,
	                                                            const dssp_skip_angle_skipping &);

	pdb_size_vec_pair backbone_complete_subset_of_pdb(const pdb &,
	                                                  const ostream_ref_opt & = ::std::nullopt,
	                                                  const dssp_skip_res_skipping & = dssp_skip_res_skipping::DONT_SKIP);

	std::pair<protein, protein_info> build_protein_of_pdb(const pdb &,
	                                                      const ostream_ref_opt & = ::std::nullopt,
	                                                      const dssp_skip_policy & = dssp_skip_policy::DONT_SKIP__DONT_BREAK_ANGLES);

	protein build_protein_of_pdb_and_name(const pdb &,
	                                      const name_set &,
	                                      const ostream_ref_opt & = ::std::nullopt);

	size_set get_protein_res_indices_that_dssp_might_skip(const pdb &,
	                                                      const ostream_ref_opt & = ::std::nullopt);

	pdb get_regions_limited_pdb(const chop::region_vec_opt &,
	                            const pdb &);

	pdb backbone_complete_region_limited_subset_of_pdb(const pdb &,
	                                                   const chop::region_vec_opt &,
	                                                   const ostream_ref_opt & = ::std::nullopt);

	/// \brief Get whether this PDB is empty of residues
	inline bool pdb::empty() const {
		return pdb_residues.empty();
	}


	/// \brief Get the number of residues held in this pdb
	inline size_t pdb::get_num_residues() const {
		return pdb_residues.size();
	}

	/// \brief Get the residue of the specified index
	///        (with no checking for which residues are backbone-complete)
	inline const pdb_residue & pdb::get_residue_of_index__backbone_unchecked(const size_t &prm_index ///< The index of the residue to retun
	                                                                         ) const {
		if constexpr ( common::IS_IN_DEBUG_MODE ) {
			if ( prm_index >= get_num_residues() ) {
				BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
				  "Unable to get_residue_ca_coord_of_index__backbone_unchecked() for index >= number of residues" ) );
			}
		}
		return pdb_residues[ prm_index ];
	}

	/// \brief Standard const begin method for the range of residues
	inline auto pdb::begin() const -> const_iterator {
		return ::std::cbegin( pdb_residues );
	}

	/// \brief Standard const end method for the range of residues
	inline auto pdb::end() const -> const_iterator {
		return ::std::cend( pdb_residues );
	}

	/// \brief Find the index of preceding residue in the same chain of the specified PDB
	///        as the residue at the specified index, or return nullopt if there is none
	///
	/// \relates pdb
	inline size_opt index_of_preceding_residue_in_same_chain(const pdb    &prm_pdb,  ///< The PDB containing the residues in question
	                                                         const size_t &prm_index ///< The index of the query residue
	                                                         ) {
		const auto &chain = get_chain_label( prm_pdb.get_residue_of_index__backbone_unchecked( prm_index ) );
		for (const size_t &index : common::indices( prm_index ) | boost::adaptors::reversed ) {
			if ( chain == get_chain_label( prm_pdb.get_residue_of_index__backbone_unchecked( index ) ) ) {
				return index;
			}
		}
		return ::std::nullopt;
	}

} // namespace cath::file

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_HPP
