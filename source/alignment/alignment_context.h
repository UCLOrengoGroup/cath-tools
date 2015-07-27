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

#ifndef ALIGNMENT_CONTEXT_H_INCLUDED
#define ALIGNMENT_CONTEXT_H_INCLUDED


#include "alignment/alignment.h"                  // for alignment
#include "common/type_aliases.h"                  // for str_vec
#include "file/pdb/pdb_list.h"                    // for pdb_list
#include "superposition/superposition_context.h"

namespace cath { namespace sup { class superposition; } }

namespace cath {
	namespace align {

		/// \brief Store a superposition along with the the context of the pdbs and ids of the structures
		///
		/// ATM, this is little more than a tuple<pdb_list, str_vec, alignment>
		class alignment_context final {
		private:
			/// \brief TODOCUMENT
			file::pdb_list pdbs;

			/// \brief TODOCUMENT
			str_vec        names;

			/// \brief TODOCUMENT
			alignment      the_alignment;

		public:
			alignment_context(const file::pdb_list &,
			                  const str_vec &,
			                  const alignment &);

			const file::pdb_list & get_pdbs()      const;
			const str_vec        & get_names()     const;
			const alignment      & get_alignment() const;
		};

		sup::superposition_context make_superposition_context(const alignment_context &,
		                                                      const sup::superposition &);

	}
}

#endif
