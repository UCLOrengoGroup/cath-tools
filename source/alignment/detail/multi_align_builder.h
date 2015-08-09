/// \file
/// \brief The multi_align_builder class header

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

#ifndef MULTI_ALIGN_BUILDER_H_INCLUDED
#define MULTI_ALIGN_BUILDER_H_INCLUDED

#include "alignment/align_type_aliases.h"
#include "alignment/detail/multi_align_group.h"
#include "alignment/refiner/alignment_refiner.h"

#include <vector>

namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { class pdb_list; } }

namespace cath {
	namespace align {
		namespace detail {

			/// \brief Build up a full alignment by gluing together alignments based on a spanning tree
			///
			/// This is an implementation class for build_alignment_from_parts()
			class multi_align_builder final {
			private:
				/// \brief A vector of the group that are being built
				///
				/// Note that this may contain inactive groups - remnants of data that
				/// has now been merged into another group. Use group_index_of_entry
				/// to find active groups
				std::vector<multi_align_group> groups;

				/// \brief A mapping from each entry's index to the index of the group
				///        in which it now resides
				size_vec                       group_index_of_entry;

				/// \brief An alignment_refiner with which to refine the alignments being built
				alignment_refiner              the_refiner;

				size_t find_group_of_entry(const size_t &) const;
				void update_group_index_of_entry(const size_t &);

			public:
				multi_align_builder(const size_t &);

				size_set get_active_groups() const;
				const multi_align_group & get_group_of_index(const size_t &) const;

				void add_alignment_branch(const size_t &,
				                          const size_t &,
				                          const alignment &,
				                          const protein_list &);

				alignment get_alignment() const;
			};

			std::ostream & operator<<(std::ostream &,
			                          const multi_align_builder &);

			void add_alignment_branch(multi_align_builder &,
			                          const size_size_alignment_tuple &,
			                          const protein_list &);
		}
	}
}

#endif
