/// \file
/// \brief The superposition_context class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef SUPERPOSITION_CONTEXT_H_INCLUDED
#define SUPERPOSITION_CONTEXT_H_INCLUDED

#include <boost/optional.hpp>

#include "alignment/alignment.h"
#include "common/type_aliases.h"
#include "file/pdb/pdb_list.h"
#include "superposition/superposition.h"

namespace cath {
	namespace align {
		class alignment_context;
	}
}

namespace cath {
	namespace sup {

		/// \brief Store a superposition along with the the context of the pdbs and ids of the actual structures being superposed
		///
		/// ATM, this is little more than a tuple<pdb_list, str_vec, superposition> with nice names and extra functionality
		/// of optionally storing an alignment
		class superposition_context final {
		private:
			file::pdb_list pdbs;
			str_vec        names;
			superposition  the_superposition;
			boost::optional<align::alignment> any_alignment;

		public:
			superposition_context(const file::pdb_list &,
			                      const str_vec &,
			                      const superposition &);
			superposition_context(const file::pdb_list &,
			                      const str_vec &,
			                      const superposition &,
			                      const align::alignment &);

			const file::pdb_list &   get_pdbs_cref()          const;
			const str_vec &          get_names_cref()         const;
			const superposition &    get_superposition_cref() const;
			bool has_alignment() const;
			const align::alignment & get_alignment_cref() const;
		};

		align::alignment_context make_alignment_context(const superposition_context &);
	}
}

#endif