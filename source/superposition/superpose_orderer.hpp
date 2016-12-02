/// \file
/// \brief The superpose_orderer class header

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

#ifndef _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSE_ORDERER_H
#define _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSE_ORDERER_H

#include "common/type_aliases.hpp"

#include <deque>
#include <vector>

namespace superpose_orderer_test_suite { struct half_matrix_indexing; }
namespace superpose_orderer_test_suite { struct half_matrix_indexing_throws_for_invalid_indices; }

namespace cath {
	namespace sup {

		/// \brief Find a maximum spanning tree over the graph of scores specified to give a sensible strategy.
		///
		/// This is useful for superposing multiple structures using pairwise superpositions.
		///
		/// The class itself just stores the half-matrix of scores. The action happens in  get_order_of_pairs_to_superpose().
		/// \todo Give this class a better name.
		class superpose_orderer final {
		public:
			using size_type = doub_vec::size_type;

		private:
			friend struct superpose_orderer_test_suite::half_matrix_indexing;
			friend struct superpose_orderer_test_suite::half_matrix_indexing_throws_for_invalid_indices;

			bool_deq  is_scored;
			doub_vec  scores;
			size_type num_items;

			static size_type half_matrix_index_of_indices(const size_type &,
			                                              const size_type &);
			void check_index_pair(const size_type &,
			                      const size_type &) const;

		public:
			explicit superpose_orderer(const size_type &);

			size_type get_num_items() const;
			bool has_score(const size_type &,
			               const size_type &) const;
			double get_score(const size_type &,
			                 const size_type &) const;

			void set_score(const size_type &,
			               const size_type &,
			               const double &);
		};

		size_size_pair_vec get_spanning_tree(const superpose_orderer &);

		size_size_pair_vec get_spanning_tree_ordered_by_desc_score(const superpose_orderer &);

		size_size_pair_vec order_spanning_tree_by_desc_score(const superpose_orderer &,
		                                                     const size_size_pair_vec &);

		superpose_orderer make_superpose_orderer(const size_size_pair_doub_map &);

		size_size_pair_vec get_spanning_tree_ordered_by_desc_score(const size_size_pair_doub_map &);
	} // namespace sup
} // namespace cath

#endif
