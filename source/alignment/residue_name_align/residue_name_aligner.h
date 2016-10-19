/// \file
/// \brief The residue_name_aligner class header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_RESIDUE_NAME_ALIGN_RESIDUE_NAME_ALIGNER_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_RESIDUE_NAME_ALIGN_RESIDUE_NAME_ALIGNER_H

#include "alignment/align_type_aliases.h"
#include "common/type_aliases.h"
#include "structure/protein/protein_list.h"
#include "structure/structure_type_aliases.h"

#include <iosfwd>
#include <set>
#include <utility>
#include <vector>

namespace cath {
	namespace align {
		class alignment;
		class residue_scorer;
	} // namespace align
	class protein_list;
} // namespace cath

namespace cath {
	namespace align {

		/// \brief A class with a public static method to construct an alignment by just pairing up residues with the same name
		///
		/// \todo Extend this code to work with more than two residue lists.
		///       Look at residue_name_aligner_test.cpp for relevant notes.
		class residue_name_aligner final {
		private:
			residue_name_aligner() = delete;

		public:
			static alignment residue_name_align(const residue_name_vec_vec &);
		};

		alignment residue_name_align_and_residue_score(const residue_name_vec_vec &,
		                                               const residue_scorer &,
		                                               const protein_list &);

		alignment residue_name_align_and_residue_score_if_multi(const residue_name_vec_vec &,
		                                                        const residue_scorer &,
		                                                        const protein_list &);

	} // namespace align
} // namespace cath

#endif
