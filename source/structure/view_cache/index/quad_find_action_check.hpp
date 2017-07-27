/// \file
/// \brief The quad_find_action_check class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_INDEX_QUAD_FIND_ACTION_CHECK_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_INDEX_QUAD_FIND_ACTION_CHECK_H

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/tools/old/impl.hpp>

#include "common/algorithm/contains.hpp"
#include "exception/out_of_range_exception.hpp"
#include "structure/protein/protein.hpp"
#include "structure/view_cache/index/detail/vcie_match_criteria.hpp"
#include "structure/view_cache/index/view_cache_index_entry.hpp"

namespace cath {
	namespace index {

		/// \brief Keep a running total score for all residue quads seen
		///
		/// This is the standard action to pass to view_cache_index::perform_action_on_all_match_at_leaves() and it just
		/// adds the score of each residue quad it's shown to total_score
		///
		/// \todo Is it necessary for this to hold references to the two proteins?
		class quad_find_action_check final {
		private:
			using idx_idx_idx_idx_tuple     = std::tuple<detail::index_type, detail::index_type, detail::index_type, detail::index_type>;
			using idx_idx_idx_idx_tuple_set = std::set<idx_idx_idx_idx_tuple>;

			/// \brief A running total of the scores for all residue quads for which this action has been performed
			double total_score = 0.0;

			/// \brief A const reference to the first protein being compared
			const protein &protein_a;

			/// \brief A const reference to the second protein being compared
			const protein &protein_b;

			/// \brief A const reference to the criteria
			const detail::vcie_match_criteria &the_criteria;

			/// \brief TODOCUMENT
			idx_idx_idx_idx_tuple_set previously_seen_indices;

		public:
			quad_find_action_check(const protein &,
			                       const protein &,
			                       const detail::vcie_match_criteria &);
			
			const double & get_total_score() const;

			void operator()(const view_cache_index_entry &,
			                const view_cache_index_entry &);

			// void operator()(const size_size_pair &,
			//                 const size_size_pair &);
		};

		/// \brief Add the score arising from the two from/to residue pairs (represented by the two specified view_cache_index_entry objects)
		///         to the running total_score
		inline void quad_find_action_check::operator()(const view_cache_index_entry &arg_cache_a, ///< A view_cache_index_entry representing the first  from/to residue pair
		                                               const view_cache_index_entry &arg_cache_b  ///< A view_cache_index_entry representing the second from/to residue pair
		                                               ) {
			const auto &a_from_index = arg_cache_a.get_from_index();
			const auto &a_to_index   = arg_cache_a.get_to_index();
			const auto &b_from_index = arg_cache_b.get_from_index();
			const auto &b_to_index   = arg_cache_b.get_to_index();

			const auto  a_indices    = size_size_pair{ a_from_index, a_to_index };
			const auto  b_indices    = size_size_pair{ b_from_index, b_to_index };
			if ( ! the_criteria( a_indices, b_indices, protein_a, protein_b ) ) {
				BOOST_THROW_EXCEPTION(common::out_of_range_exception(
					"The match_at_nodes scan has accepted an entry pair that is then rejected by standard comparisons "
					+ boost::lexical_cast<std::string>( arg_cache_a )
					+ " and "
					+ boost::lexical_cast<std::string>( arg_cache_b )
				));
			}

			const auto the_indices = std::make_tuple( a_from_index, a_to_index, b_from_index, b_to_index );
			if ( common::contains( previously_seen_indices, the_indices ) ) {
				BOOST_THROW_EXCEPTION(common::out_of_range_exception("Actioning indices that have already been seen"));
			}
			previously_seen_indices.insert( the_indices );


			const double sq_dist       = detail::squared_distance( arg_cache_a, arg_cache_b );
			const double sq_dist_check = detail::squared_distance( a_indices, b_indices, protein_a, protein_b );
			if ( ! boost::test_tools::check_is_close( sq_dist, sq_dist_check, boost::math::fpc::percent_tolerance( 0.001 ) ) ) {
				BOOST_THROW_EXCEPTION(common::out_of_range_exception(
					"Squared distance is different (from vcies: "
					+ std::to_string( sq_dist )
					+ "; recalculated: "
					+ std::to_string( sq_dist_check )
					+ ")"
				));
			}
			const float_score_type score = static_cast<float_score_type>( 1.0 ) - ( sqrt( sq_dist ) / static_cast<float_score_type>( 7.0 ) );
//			const double score   = ( 10.0 / ( sq_dist + 10.0 ) );
			total_score += score;
		}

// 		/// \brief Add the score arising from the two from/to residue pairs (represented by the two specified from/to index pairs)
// 		///         to the running total_score
// 		inline void quad_find_action_check::operator()(const size_size_pair &arg_indices_a, ///< The indices of the first  from/to residue pair in their protein
// 		                                               const size_size_pair &arg_indices_b  ///< The indices of the second from/to residue pair in their protein
// 		                                               ) {
// 			const double sq_dist = detail::squared_distance( arg_indices_a, arg_indices_b, protein_a, protein_b );
// 			const float_score_type score = static_cast<float_score_type>( 1.0 ) - ( sqrt( sq_dist ) / static_cast<float_score_type>( 7.0 ) );
// //			const double score   = ( 10.0 / ( sq_dist + 10.0 ) );
// 			total_score += score;
// 		}

	} // namespace index
} // namespace cath

#endif
