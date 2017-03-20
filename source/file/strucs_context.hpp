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
#include "common/type_aliases.hpp"                  // for str_vec
#include "exception/invalid_argument_exception.hpp"
#include "file/pdb/pdb_list.hpp"                    // for pdb_list

namespace cath {
	namespace file {

		/// \brief Store the context of structures: the PDBs, IDs and any regions
		class strucs_context final {
		private:
			/// \brief The PDBs of the structures
			file::pdb_list           pdbs;

			/// \brief The IDs of the structures
			str_vec                  names;

			/// \brief The key regions of each structure to which this refers (or none where it refers to all of a structure)
			chop::region_vec_opt_vec regions;

		public:
			strucs_context(const file::pdb_list,
			               const str_vec,
			               const chop::region_vec_opt_vec);

			const file::pdb_list & get_pdbs() const;
			const str_vec & get_names() const;
			const chop::region_vec_opt_vec & get_regions() const;

			strucs_context & set_pdbs(const file::pdb_list &);
		};

		size_t get_num_entries(const strucs_context &);

		/// \brief Ctor for strucs_context
		inline strucs_context::strucs_context(const pdb_list                 arg_pdbs,   ///< The PDBs of the structures
		                                      const str_vec                  arg_names,  ///< The IDs of the structures
		                                      const chop::region_vec_opt_vec arg_regions ///< The key regions of each structure to which this refers (or none where it refers to all of a structure)
		                                      ) : pdbs    { std::move( arg_pdbs    ) },
		                                          names   { std::move( arg_names   ) },
		                                          regions { std::move( arg_regions ) } {
			if ( pdbs.size() != names.size() || pdbs.size() != regions.size() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot construct a strucs_context with inconsistent numbers of PDBs / names / regions"));
			}
		}

		/// \brief Getter for the PDBs of the structures
		inline const file::pdb_list & strucs_context::get_pdbs() const {
			return pdbs;
		}

		/// \brief Getter for the IDs of the structures
		inline const str_vec & strucs_context::get_names() const {
			return names;
		}

		/// \brief Getter for the key regions of each structure to which this refers (or none where it refers to all of a structure)
		inline const chop::region_vec_opt_vec & strucs_context::get_regions() const {
			return regions;
		}

		/// \brief Set the PDBs to the specified PDBs
		inline strucs_context & strucs_context::set_pdbs(const file::pdb_list &arg_pdbs ///< The PDBs to set
		                                                 ) {
			if ( arg_pdbs.size() != get_num_entries( *this ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot set an inconsistent number of PDBs in a strucs_context"));
			}
			pdbs = arg_pdbs;
			return *this;
		}

		/// \brief Get the number of entries in the specified strucs_context
		inline size_t get_num_entries(const strucs_context &arg_strucs_context ///< The strucs_context to query
		                              ) {
			return arg_strucs_context.get_names().size();
		}

	} // namespace file
} // namespace cath

#endif
