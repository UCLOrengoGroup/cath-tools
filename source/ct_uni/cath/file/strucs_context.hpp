/// \file
/// \brief The strucs_context class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_STRUCS_CONTEXT_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_STRUCS_CONTEXT_HPP

#include <optional>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"                  // for str_vec
#include "cath/file/name_set/name_set_list.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_list.hpp"                    // for pdb_list

namespace cath {
	namespace file {

		/// \brief Store the context of structures: the PDBs, IDs and any regions
		class strucs_context final {
		private:
			/// \brief The PDBs of the structures
			file::pdb_list           pdbs;

			/// \brief The IDs of the structures
			name_set_list            name_sets;

			/// \brief The key regions of each structure to which this refers (or nullopt where it refers to all of a structure)
			chop::region_vec_opt_vec regions;

		public:
			explicit strucs_context(file::pdb_list);
			strucs_context(file::pdb_list,
			               name_set_list);
			strucs_context(file::pdb_list,
			               name_set_list,
			               chop::region_vec_opt_vec);

			const file::pdb_list & get_pdbs() const;
			const name_set_list & get_name_sets() const;
			const chop::region_vec_opt_vec & get_regions() const;

			strucs_context & set_pdbs(const file::pdb_list &);
		};

		size_t size(const strucs_context &);

		/// \brief Ctor for strucs_context
		inline strucs_context::strucs_context(pdb_list prm_pdbs ///< The PDBs of the structures
		                                      ) : strucs_context{
		                                          	prm_pdbs,
		                                          	name_set_list           ( prm_pdbs.size() ),
		                                          	chop::region_vec_opt_vec( prm_pdbs.size() )
		                                          } {
		}

		/// \brief Ctor for strucs_context
		inline strucs_context::strucs_context(pdb_list      prm_pdbs,     ///< The PDBs of the structures
		                                      name_set_list prm_name_sets ///< The IDs of the structures
		                                      ) : strucs_context{
		                                          	prm_pdbs,
		                                          	prm_name_sets,
		                                          	chop::region_vec_opt_vec( prm_pdbs.size() )
		                                          } {
		}

		/// \brief Ctor for strucs_context
		inline strucs_context::strucs_context(pdb_list                 prm_pdbs,      ///< The PDBs of the structures
		                                      name_set_list            prm_name_sets, ///< The IDs of the structures
		                                      chop::region_vec_opt_vec prm_regions    ///< The key regions of each structure to which this refers (or nullopt where it refers to all of a structure)
		                                      ) : pdbs      { std::move( prm_pdbs      ) },
		                                          name_sets { std::move( prm_name_sets ) },
		                                          regions   { std::move( prm_regions   ) } {
			if ( pdbs.size() != name_sets.size() || pdbs.size() != regions.size() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot construct a strucs_context with inconsistent numbers of PDBs / names / regions"));
			}
		}

		/// \brief Getter for the PDBs of the structures
		inline const file::pdb_list & strucs_context::get_pdbs() const {
			return pdbs;
		}

		/// \brief Getter for the IDs of the structures
		inline const name_set_list & strucs_context::get_name_sets() const {
			return name_sets;
		}

		/// \brief Getter for the key regions of each structure to which this refers (or nullopt where it refers to all of a structure)
		inline const chop::region_vec_opt_vec & strucs_context::get_regions() const {
			return regions;
		}

		/// \brief Set the PDBs to the specified PDBs
		inline strucs_context & strucs_context::set_pdbs(const file::pdb_list &prm_pdbs ///< The PDBs to set
		                                                 ) {
			if ( prm_pdbs.size() != size( *this ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot set an inconsistent number of PDBs in a strucs_context"));
			}
			pdbs = prm_pdbs;
			return *this;
		}

		/// \brief Get the number of entries in the specified strucs_context
		inline size_t size(const strucs_context &prm_strucs_context ///< The strucs_context to query
		                   ) {
			return prm_strucs_context.get_name_sets().size();
		}

		strucs_context strucs_context_of_backbone_complete_subset_pdbs(const strucs_context &,
		                                                               const ostream_ref_opt & = ::std::nullopt);

		strucs_context strucs_context_of_backbone_complete_region_limited_subset_pdbs(const strucs_context &,
		                                                                              const ostream_ref_opt & = ::std::nullopt);

		void restrict_pdbs(strucs_context &);

		strucs_context restrict_pdbs_copy(strucs_context);

		size_t get_num_regions_set(const strucs_context &);

		std::string to_string(const strucs_context &);

		file::pdb_list get_restricted_pdbs(const strucs_context &);

		protein_list build_protein_list(const strucs_context &);

		chop::domain_opt get_domain_opt_of_index(const strucs_context &,
		                                         const size_t &);

		void non_crypto_hash(size_t &,
		                     const strucs_context &);

		size_t non_crypto_hash_copy(size_t,
		                            const strucs_context &);

	} // namespace file
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_STRUCS_CONTEXT_HPP
