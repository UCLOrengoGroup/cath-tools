/// \file
/// \brief The alignment action header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_ALIGNMENT_ACTION_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_ALIGNMENT_ACTION_H

#include "alignment/align_type_aliases.h"

namespace cath { namespace align { class alignment; } }
namespace cath { class protein_list; }

namespace cath {
	namespace align {
		namespace detail {
			using aln_ent_ind_tup      = std::tuple<const alignment &, const size_t &, const size_t &>;

			using aln_ent_ind_tup_pair = std::pair<aln_ent_ind_tup, aln_ent_ind_tup>;

			enum class glued_row_type {
				FROM_A,
				FROM_B,
				FROM_BOTH
			};

			void append_glued_row(alignment &,
			                      const aln_ent_ind_tup_pair &,
			                      const glued_row_type &);
		} // namespace detail

		alignment glue_two_alignments(const alignment &,
		                              const size_t &,
		                              const alignment &,
		                              const size_t & = 0);

		alignment build_alignment_from_parts(const size_size_alignment_tuple_vec &,
		                                     const protein_list &);

	} // namespace align
} // namespace cath

#endif
