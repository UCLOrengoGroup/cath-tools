/// \file
/// \brief The alignment_context class header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_ALIGNMENT_CONTEXT_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_ALIGNMENT_CONTEXT_H

#include "alignment/alignment.hpp"                  // for alignment
#include "chopping/chopping_type_aliases.hpp"
#include "common/type_aliases.hpp"                  // for str_vec
#include "file/pdb/pdb_list.hpp"                    // for pdb_list

namespace cath { namespace sup { class superposition; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	namespace align {

		/// \brief Store an alignment along with the context of the PDBs and ids of the structures
		///
		/// ATM, this is little more than a tuple<pdb_list, str_vec, alignment>
		class alignment_context final {
		private:
			/// \brief TODOCUMENT
			file::pdb_list           pdbs;

			/// \brief TODOCUMENT
			str_vec                  names;

			/// \brief TODOCUMENT
			alignment                the_alignment;

			/// \brief For each PDB: the regions of the PDB to which the alignment refers,
			///        or none if it refers to all of it
			chop::region_vec_opt_vec regions;

		public:
			alignment_context(const file::pdb_list &,
			                  const str_vec &,
			                  const alignment &,
			                  const chop::region_vec_opt_vec &);

			const file::pdb_list & get_pdbs() const;
			const str_vec & get_names() const;
			const alignment & get_alignment() const;
			const chop::region_vec_opt_vec & get_regions() const;
		};

		sup::superposition_context make_superposition_context(const alignment_context &,
		                                                      const sup::superposition &);

	} // namespace align
} // namespace cath

#endif
