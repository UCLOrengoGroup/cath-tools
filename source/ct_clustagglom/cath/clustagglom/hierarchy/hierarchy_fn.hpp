/// \file
/// \brief The hierarchy functions header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_FN_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_FN_HPP

#include <functional>

#include "cath/clustagglom/hierarchy.hpp"
#include "cath/clustagglom/hierarchy/hierarchy_group.hpp"
#include "cath/clustagglom/hierarchy/hierarchy_value.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"

namespace cath::clust::detail {

	/// \brief Provide depth-first traversal of a hierarchy
	class depth_first_hierachy_traverser final {
	private:
		/// \brief The initial depth from which traversal should start
		static constexpr size_t   INIT_DEPTH = 0;

		/// \brief The initial group index from which traversal of each layer should start
		static constexpr size_t   INIT_GROUP = 0;

		/// \brief The index from which the counters should start (eg like C in CATH starts at 1)
		static constexpr item_idx CTR_INIT   = 1;

		/// \brief The counters (like CATHSOLID) during the traversal
		item_vec ctrs;

		/// \brief Recursively traverse the specified hierarchy from the specified depth and group
		///        calling the specified function at each entry
		template <typename Fn>
		void traverse_recurse(const hierarchy  &prm_hierarchy, ///< The hierarchy to traverse
		                      const size_t     &prm_depth,     ///< The depth from which the traversal should be conducted
		                      const size_t     &prm_group,     ///< The index of the group at the specified depth from which the traversal should be conducted
		                      Fn              &&prm_fn         ///< The function to call at each entry
		                      ) {
			const hierarchy_group &entry = prm_hierarchy[ prm_depth ][ prm_group ];

			// Add a new layer of counter
			ctrs.push_back( CTR_INIT );

			// Loop through the values in the group
			for (const size_t &value_idx : common::indices( entry.size() ) ) {
				const hierarchy_value &value = entry[ value_idx ];

				// If this value refers to an entry, temporarily fill out the counters and call the function
				if ( value.get_type() == hierarchy_ref::ENTRY ) {
					ctrs.resize( prm_hierarchy.size(), CTR_INIT );
					::std::invoke( prm_fn, ctrs, value.get_index() );
					ctrs.resize( prm_depth + 1 );
				}
				// Otherwise recurse into traversing within the child group
				else {
					traverse_recurse(
						prm_hierarchy,
						prm_depth + 1,
						value.get_index(),
						prm_fn
					);
				}

				// Increment the last counter
				++ctrs.back();
			}

			// Remove the last counter
			ctrs.pop_back();
		}

	public:

		/// \brief Perform a depth-first traversal through the specified hierarchy and call
		///        the specified function with the (CATHSOLID-like) counters and index for each leaf node
		template <typename Fn>
		void traverse(const hierarchy  &prm_hierarchy, ///< The hierarchy to traverse
		              Fn              &&prm_fn         ///< The function to apply to each node
		              ) {
			ctrs.clear();
			ctrs.reserve( prm_hierarchy.size() );
			traverse_recurse( prm_hierarchy, INIT_DEPTH, INIT_GROUP, std::forward<Fn>( prm_fn ) );
		}

	};

	/// \brief Perform a depth-first traversal through the specified hierarchy and call
	///        the specified function with the (CATHSOLID-like) counters and index for each leaf node
	template <typename Fn>
	void depth_first_traverse_hierachy(const hierarchy  &prm_hierarchy, ///< The hierarchy to traverse
	                                   Fn              &&prm_fn         ///< The function to apply to each node
	                                   ) {
		depth_first_hierachy_traverser traverser;
		traverser.traverse( prm_hierarchy, std::forward<Fn>( prm_fn ) );
	}

} // namespace cath::clust::detail

#endif // CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_HIERARCHY_HIERARCHY_FN_HPP
