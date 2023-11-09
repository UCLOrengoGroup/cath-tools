/// \file
/// \brief The pdb_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_LIST_HPP

#include <iostream>
#include <optional>
#include <vector>

#include <boost/operators.hpp>
#include <boost/range.hpp>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/structure/structure_type_aliases.hpp"

// clang-format off
namespace cath { class protein_list; }
namespace cath::file { class name_set_list; }
namespace cath::file { class pdb; }
// clang-format on

namespace cath::file {

	/// \brief TODOCUMENT
	class pdb_list final {
	private:
		/// \brief TODOCUMENT
		pdb_vec pdbs;

	public:
		pdb_list() = default;
		explicit pdb_list(pdb_vec);

		void push_back( const pdb & );
		void push_back( pdb && );

		void reserve(const size_t &);

		[[nodiscard]] size_t size() const;
		[[nodiscard]] bool   empty() const;

		pdb & operator[](const size_t &);
		const pdb & operator[](const size_t &) const;

		// Provide iterators to make this into a range
		using const_iterator = std::vector<pdb>::const_iterator;
		using iterator       = const_iterator;
//		iterator begin();
//		iterator end();
		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	/// \brief Ctor from a vector<pdb>
	///
	/// \param prm_pdbs The pdbs from which this pdb_list should be constructed
	inline pdb_list::pdb_list( pdb_vec prm_pdbs ) : pdbs{ std::move( prm_pdbs ) } {
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_pdb TODOCUMENT
	inline void pdb_list::push_back( const pdb &prm_pdb ) {
		pdbs.push_back( prm_pdb );
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_pdb TODOCUMENT
	inline void pdb_list::push_back( pdb &&prm_pdb ) {
		pdbs.push_back( ::std::move( prm_pdb ) );
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_size TODOCUMENT
	inline void pdb_list::reserve( const size_t &prm_size ) {
		pdbs.reserve( prm_size );
	}

	/// \brief TODOCUMENT
	inline size_t pdb_list::size() const {
		return pdbs.size();
	}

	/// \brief TODOCUMENT
	inline bool pdb_list::empty() const {
		return pdbs.empty();
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_index TODOCUMENT
	inline pdb &pdb_list::operator[]( const size_t &prm_index ) {
		return pdbs[ prm_index ];
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_index TODOCUMENT
	inline const pdb &pdb_list::operator[]( const size_t &prm_index ) const {
		return pdbs[ prm_index ];
	}

	/// \brief TODOCUMENT
	inline auto pdb_list::begin() const -> const_iterator {
		return cbegin( pdbs );
	}

	/// \brief TODOCUMENT
	inline auto pdb_list::end() const -> const_iterator {
		return cend( pdbs );
	}

	pdb_list read_pdb_files(const path_vec &);

	pdb_list pdb_list_of_backbone_complete_subset_pdbs(const pdb_list &,
	                                                   const ostream_ref_opt & = ::std::nullopt);

	pdb_list pdb_list_of_backbone_complete_region_limited_subset_pdbs(const pdb_list &,
	                                                                  const chop::region_vec_opt_vec &,
	                                                                  const ostream_ref_opt & = ::std::nullopt);

	protein_list build_protein_list_of_pdb_list( const pdb_list &, const ostream_ref_opt & = ::std::nullopt );

	protein_list build_protein_list_of_pdb_list_and_names(const pdb_list &,
	                                                      const name_set_list &);

	amino_acid_vec_vec get_amino_acid_lists(const pdb_list &);

	residue_id_vec_vec get_residue_ids_of_first_chains__backbone_unchecked(const pdb_list &);

	backbone_complete_indices_vec get_backbone_complete_indices(const pdb_list &);

	residue_id_vec_vec get_backbone_complete_residue_ids(const pdb_list &);

	residue_id_vec_vec get_backbone_complete_residue_ids_of_first_chains(const pdb_list &);

} // namespace cath::file

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PDB_LIST_HPP
