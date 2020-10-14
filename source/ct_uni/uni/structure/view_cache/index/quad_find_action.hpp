/// \file
/// \brief The quad_find_action class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_INDEX_QUAD_FIND_ACTION_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_INDEX_QUAD_FIND_ACTION_H

#include "structure/protein/protein.hpp"
#include "structure/view_cache/index/view_cache_index_entry.hpp"

namespace cath {
	namespace index {

		/// \brief Keep a running total score for all residue quads seen
		///
		/// This is the standard action to pass to view_cache_index::perform_action_on_all_match_at_leaves() and it just
		/// adds the score of each residue quad it's shown to total_score
		///
		/// \todo Is it necessary for this to hold references to the two proteins?
		class quad_find_action final {
		private:
			/// \brief A running total of the scores for all residue quads for which this action has been performed
			double total_score = 0.0;

			/// \brief A const reference to the first protein being compared
			const protein &protein_a;

			/// \brief A const reference to the second protein being compared
			const protein &protein_b;

		public:
			quad_find_action(const protein &,
			                 const protein &);
			
			const double & get_total_score() const;

			void operator()(const view_cache_index_entry &,
			                const view_cache_index_entry &);

			void operator()(const size_size_pair &,
			                const size_size_pair &);
		};

		/// \brief Add the score arising from the two from/to residue pairs (represented by the two specified view_cache_index_entry objects)
		///         to the running total_score
		inline void quad_find_action::operator()(const view_cache_index_entry &prm_cache_a, ///< A view_cache_index_entry representing the first  from/to residue pair
		                                         const view_cache_index_entry &prm_cache_b  ///< A view_cache_index_entry representing the second from/to residue pair
		                                         ) {
			const double sq_dist = detail::squared_distance( prm_cache_a, prm_cache_b );
			const float_score_type score = 1.0 - ( sqrt( sq_dist ) / 7.0 );
//			const double score   = ( 10.0 / ( sq_dist + 10.0 ) );
			total_score += score;
		}

		/// \brief Add the score arising from the two from/to residue pairs (represented by the two specified from/to index pairs)
		///         to the running total_score
		inline void quad_find_action::operator()(const size_size_pair &prm_indices_a, ///< The indices of the first  from/to residue pair in their protein
		                                         const size_size_pair &prm_indices_b  ///< The indices of the second from/to residue pair in their protein
		                                         ) {
			const double sq_dist = detail::squared_distance( prm_indices_a, prm_indices_b, protein_a, protein_b );
			const float_score_type score = 1.0 - ( sqrt( sq_dist ) / 7.0 );
//			const double score   = ( 10.0 / ( sq_dist + 10.0 ) );
			total_score += score;
		}

	} // namespace index
} // namespace cath

#endif
