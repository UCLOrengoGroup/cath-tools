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

#ifndef _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_HIERARCHY_HIERARCHY_FN_HPP
#define _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_HIERARCHY_HIERARCHY_FN_HPP

#include "clustagglom/hierarchy.hpp"
#include "clustagglom/hierarchy/hierarchy_group.hpp"
#include "clustagglom/hierarchy/hierarchy_value.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/indices.hpp"

namespace cath {
	namespace clust {
		namespace detail {

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
				void traverse_recurse(const hierarchy  &arg_hierarchy, ///< The hierarchy to traverse
				                      const size_t     &arg_depth,     ///< The depth from which the traversal should be conducted
				                      const size_t     &arg_group,     ///< The index of the group at the specified depth from which the traversal should be conducted
				                      Fn              &&arg_fn         ///< The function to call at each entry
				                      ) {
					const hierarchy_group &entry = arg_hierarchy[ arg_depth ][ arg_group ];

					// Add a new layer of counter
					ctrs.push_back( CTR_INIT );

					// Loop through the values in the group
					for (const size_t &value_idx : common::indices( entry.size() ) ) {
						const hierarchy_value &value = entry[ value_idx ];

						// If this value refers to an entry, temporarily fill out the counters and call the function
						if ( value.get_type() == hierarchy_ref::ENTRY ) {
							ctrs.resize( arg_hierarchy.size(), CTR_INIT );
							common::invoke( arg_fn, ctrs, value.get_index() );
							ctrs.resize( arg_depth + 1 );
						}
						// Otherwise recurse into traversing within the child group
						else {
							traverse_recurse(
								arg_hierarchy,
								arg_depth + 1,
								value.get_index(),
								arg_fn
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
				void traverse(const hierarchy  &arg_hierarchy, ///< The hierarchy to traverse
				              Fn              &&arg_fn         ///< The function to apply to each node
				              ) {
					ctrs.clear();
					ctrs.reserve( arg_hierarchy.size() );
					traverse_recurse( arg_hierarchy, INIT_DEPTH, INIT_GROUP, std::forward<Fn>( arg_fn ) );
				}

			};

			/// \brief Perform a depth-first traversal through the specified hierarchy and call
			///        the specified function with the (CATHSOLID-like) counters and index for each leaf node
			template <typename Fn>
			void depth_first_traverse_hierachy(const hierarchy  &arg_hierarchy, ///< The hierarchy to traverse
			                                   Fn              &&arg_fn         ///< The function to apply to each node
			                                   ) {
				depth_first_hierachy_traverser traverser;
				traverser.traverse( arg_hierarchy, std::forward<Fn>( arg_fn ) );
			}

		} // namespace detail

	} // namespace clust
} // namespace cath

#endif
