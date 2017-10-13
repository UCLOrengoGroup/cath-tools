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

#include "alignment/alignment.hpp" // for alignment
#include "file/strucs_context.hpp"

namespace cath { namespace file { class name_set_list; } }
namespace cath { namespace sup { class superposition; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	namespace align {

		/// \brief Store an alignment along with the context of the PDBs and ids of the actual structures being superposed
		class alignment_context final {
		private:
			/// \brief TODOCUMENT
			alignment            the_alignment;

			/// \brief TODOCUMENT
			file::strucs_context context;

		public:
			alignment_context(alignment,
			                  file::strucs_context);

			alignment_context(alignment,
			                  const file::pdb_list &,
			                  const file::name_set_list &,
			                  const chop::region_vec_opt_vec &);

			const alignment & get_alignment() const;
			const file::strucs_context & get_strucs_context() const;
		};

		alignment_context make_restricted_alignment_context(alignment,
		                                                    file::strucs_context);

		const file::pdb_list & get_pdbs(const alignment_context &);
		const file::name_set_list & get_name_sets(const alignment_context &);
		const chop::region_vec_opt_vec & get_regions(const alignment_context &);

		file::pdb_list get_restricted_pdbs(const alignment_context &);

		sup::superposition_context make_superposition_context(const alignment_context &,
		                                                      const sup::superposition &);

	} // namespace align
} // namespace cath

#endif
