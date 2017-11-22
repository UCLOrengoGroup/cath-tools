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

#ifndef _CATH_TOOLS_SOURCE_FILE_STRUCS_CONTEXT_H
#define _CATH_TOOLS_SOURCE_FILE_STRUCS_CONTEXT_H

#include "chopping/chopping_type_aliases.hpp"
#include "chopping/region/region.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/type_aliases.hpp"                  // for str_vec
#include "file/name_set/name_set_list.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_list.hpp"                    // for pdb_list

namespace cath {
	namespace file {

		/// \brief Store the context of structures: the PDBs, IDs and any regions
		class strucs_context final {
		private:
			/// \brief The PDBs of the structures
			file::pdb_list           pdbs;

			/// \brief The IDs of the structures
			name_set_list            name_sets;

			/// \brief The key regions of each structure to which this refers (or none where it refers to all of a structure)
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
		inline strucs_context::strucs_context(pdb_list arg_pdbs ///< The PDBs of the structures
		                                      ) : strucs_context{
		                                          	arg_pdbs,
		                                          	name_set_list           ( arg_pdbs.size() ),
		                                          	chop::region_vec_opt_vec( arg_pdbs.size() )
		                                          } {
		}

		/// \brief Ctor for strucs_context
		inline strucs_context::strucs_context(pdb_list      arg_pdbs,     ///< The PDBs of the structures
		                                      name_set_list arg_name_sets ///< The IDs of the structures
		                                      ) : strucs_context{
		                                          	arg_pdbs,
		                                          	arg_name_sets,
		                                          	chop::region_vec_opt_vec( arg_pdbs.size() )
		                                          } {
		}

		/// \brief Ctor for strucs_context
		inline strucs_context::strucs_context(pdb_list                 arg_pdbs,      ///< The PDBs of the structures
		                                      name_set_list            arg_name_sets, ///< The IDs of the structures
		                                      chop::region_vec_opt_vec arg_regions    ///< The key regions of each structure to which this refers (or none where it refers to all of a structure)
		                                      ) : pdbs      { std::move( arg_pdbs      ) },
		                                          name_sets { std::move( arg_name_sets ) },
		                                          regions   { std::move( arg_regions   ) } {
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

		/// \brief Getter for the key regions of each structure to which this refers (or none where it refers to all of a structure)
		inline const chop::region_vec_opt_vec & strucs_context::get_regions() const {
			return regions;
		}

		/// \brief Set the PDBs to the specified PDBs
		inline strucs_context & strucs_context::set_pdbs(const file::pdb_list &arg_pdbs ///< The PDBs to set
		                                                 ) {
			if ( arg_pdbs.size() != size( *this ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot set an inconsistent number of PDBs in a strucs_context"));
			}
			pdbs = arg_pdbs;
			return *this;
		}

		/// \brief Get the number of entries in the specified strucs_context
		inline size_t size(const strucs_context &arg_strucs_context ///< The strucs_context to query
		                   ) {
			return arg_strucs_context.get_name_sets().size();
		}

		strucs_context strucs_context_of_backbone_complete_subset_pdbs(const strucs_context &,
		                                                               const ostream_ref_opt & = boost::none);

		strucs_context strucs_context_of_backbone_complete_region_limited_subset_pdbs(const strucs_context &,
		                                                                              const ostream_ref_opt & = boost::none);

		void restrict_pdbs(strucs_context &);

		strucs_context restrict_pdbs_copy(strucs_context);

		size_t get_num_regions_set(const strucs_context &);

		std::string to_string(const strucs_context &);

		file::pdb_list get_restricted_pdbs(const strucs_context &);

		protein_list build_protein_list(const strucs_context &);

	} // namespace file
} // namespace cath

#endif
