/// \file
/// \brief The merge class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_MERGE_H
#define _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_MERGE_H

#include <boost/filesystem/path.hpp> // ***** Only required for merge_vec function *****

#include "clustagglom/clustagglom_type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief Represent a merge between to items (in the clustering sense)
		///
		/// TODO: Prevent node_a     == node_b?
		/// TODO: Prevent merge_node <= node_a?
		/// TODO: Prevent merge_node <= node_b?
		struct merge final {
			/// \brief The index of the first item to merge
			item_idx node_a;

			/// \brief The index of the second item to merge
			item_idx node_b;

			/// \brief The new number assigned for the newly created merge
			///
			/// This must be strictly greater than node_a or node_b
			item_idx merge_node;

			/// \brief The dissimilarity associated with this
			strength dissim;

			merge(const item_idx &,
			      const item_idx &,
			      const item_idx &,
			      const strength &);
		};

		/// \brief Ctor like aggregate-initialization
		inline merge::merge(const item_idx &arg_node_a,     ///< The index of the first item to merge
		                    const item_idx &arg_node_b,     ///< The index of the second item to merge
		                    const item_idx &arg_merge_node, ///< The new number assigned for the newly created merge
		                    const strength &arg_dissim      ///< The dissimilarity associated with this
		                    ) : node_a    { arg_node_a     },
		                        node_b    { arg_node_b     },
		                        merge_node{ arg_merge_node },
		                        dissim    { arg_dissim     } {
		}

		/// \brief Standard equality operator for merge
		///
		/// \relates merge
		///
		/// Note that this requires the merge_node to be the same but it may be
		/// worth having some operation that ignores that one field
		inline bool operator==(const merge &arg_lhs, ///< The first  merge to compare
		                       const merge &arg_rhs  ///< The second merge to compare
		                       ) {
			return (
				arg_lhs.node_a     == arg_rhs.node_a
				&&
				arg_lhs.node_b     == arg_rhs.node_b
				&&
				arg_lhs.merge_node == arg_rhs.merge_node
				&&
				arg_lhs.dissim     == arg_rhs.dissim
			);
		}

		std::string to_string(const merge &);
		std::string to_string(const merge_vec &);
		void write_merge_list(std::ostream &,
		                      const merge_vec &);
		void write_merge_list(const boost::filesystem::path &,
		                      const merge_vec &);
		merge_vec read_merge_list(std::istream &);
		merge_vec read_merge_list(const boost::filesystem::path &);

		std::ostream & operator<<(std::ostream &,
		                          const merge &);

	} // namespace clust
} // namespace cath

#endif
