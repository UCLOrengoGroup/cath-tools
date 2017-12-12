/// \file
/// \brief The hierarchy class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_HIERARCHY_H
#define _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_HIERARCHY_H

#include <boost/filesystem/path.hpp>

#include "clustagglom/clustagglom_type_aliases.hpp"
#include "clustagglom/hierarchy/hierarchy_layer.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/cpp17/invoke.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"

using namespace cath::common::literals;

namespace cath { namespace common { class id_of_str_bidirnl; } }

namespace cath {
	namespace clust {

		/// \brief A hierarchy (like CATH-SOLID) for a bunch of entries
		///
		/// This stores the structure in a bunch of layers which each contain a bunch of
		/// groups which each contain a bunch of values (rather than storing CATH-SOLID-like indices).
		///
		/// The root level of the hierarchy is first.
		///
		/// \invariant There must be at least one heirarchy_layer
		///
		/// \invariant The first heirarchy_layer must always contain exactly one hierarchy_group
		///            (like the virtual root node before the C in the CATH hierarchy)
		///
		/// \invariant The values in the groups of the last hierarchy_layer must all refer to entries
		///
		/// \invariant The indices in hierarchy_ref::CLUSTER must always be less than the number of groups
		///            in the following layer.
		///
		/// \invariant The indices associated with the hierarchy_ref::CLUSTER should be a permutation
		///            of the indices of the groups in the next layer
		///
		/// At each depth in the hierarchy, a hierarchy value can either refer to an entry
		/// or to a group at the next deepest level. (Though all entries in the last level must
		/// of course refer to entries).
		///
		/// \todo Consider adding a parameter to the ctor to specify whether to add preceding layer to group
		///       (and throw if not but top layer has != 1 groups (allow completely empty?))
		struct hierarchy final {
		private:
			/// \brief The layers in the hierarchy
			hierarchy_layer_vec layers;

		public:
			/// \brief A const_iterator type alias as part of making this a range over hierarchy_layers
			using const_iterator = hierarchy_layer_vec::const_iterator;

			hierarchy();
			explicit hierarchy(hierarchy_layer_vec);

			size_t size() const;

			hierarchy_layer & operator[](const size_t &);
			const hierarchy_layer & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Construct one layer with one empty group
		inline hierarchy::hierarchy() : layers{ { hierarchy_layer{ { hierarchy_group{} } } } } {
		}

		/// \brief Ctor from the specified vector of hierarchy_layers
		///
		/// \pre `! arg_layers.empty()`
		///
		/// \pre `arg_layers[ 0 ].size() == 1`
		inline hierarchy::hierarchy(hierarchy_layer_vec arg_layers ///< The hierarchy_layers from which to construct this hierarchy
		                            ) : layers{ std::move( arg_layers ) } {
			if ( layers.empty() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create a hierarchy with no layers"));
			}
			if ( common::front( layers ).size() != 1 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("First layer in hierarchy must contain exactly one group"));
			}
		}

		/// \brief Return the number of hierarchy_layers (ie the depth of the hierarchy)
		inline size_t hierarchy::size() const {
			return layers.size();
		}

		/// \brief Get the hierarchy_layer associated with the sequence with the specified index
		inline hierarchy_layer & hierarchy::operator[](const size_t &arg_index ///< The index of the hierarchy_layer to retrieve
		                                               ) {
			return layers[ arg_index ];
		}

		/// \brief Get the hierarchy_layer associated with the sequence with the specified index
		inline const hierarchy_layer & hierarchy::operator[](const size_t &arg_index ///< The index of the hierarchy_layer to retrieve
		                                                     ) const {
			return layers[ arg_index ];
		}

		/// \brief Standard const begin() method, as part of making this into a range over the hierarchy_layer entries
		inline auto hierarchy::begin() const -> const_iterator {
			return common::cbegin( layers );
		}

		/// \brief Standard const end() method, as part of making this into a range over the hierarchy_layer entries
		inline auto hierarchy::end() const -> const_iterator {
			return common::cend( layers );
		}

		namespace detail {

			item_idx first_index(const hierarchy &,
			                     const hierarchy_value &,
			                     const size_t &);

		} // namespace detail

		std::string to_string(const hierarchy &);

		std::ostream & operator<<(std::ostream &,
		                          const hierarchy &);

		void sort_hierarchy(hierarchy &,
		                    const size_vec &);

		hierarchy sort_hierarchy_copy(hierarchy,
		                              const size_vec &);

		std::ostream & write_cluster(std::ostream &,
		                             const hierarchy &,
		                             const common::id_of_str_bidirnl &);

		void write_cluster(const boost::filesystem::path &,
		                   const hierarchy &,
		                   const common::id_of_str_bidirnl &);

		hierarchy make_hierarchy_from_reversed_layers(hierarchy_layer_vec);

		/// \brief Make a hierarchy by adding a root layer (including the specified leaf indices)
		///        to the specified reverse-order layers (ie deepest first) and then reversing them
		///
		/// \pre `reverse(arg_layers)` must still make a valid hierarchy
		///
		/// \relates hierarchy
		template <typename Rng>
		hierarchy make_hierarchy_from_reversed_without_root_layer(hierarchy_layer_vec  arg_layers,      ///< The reverse-order layers from which to make the hierarchy (ie deepest first)
		                                                          const Rng           &arg_leaf_indices ///< The indices of any leaf items that must be added as separate items in the final layer's group
		                                                          ) {
			// Create a single group for the final layer
			hierarchy_group new_group;

			// If there are previous layers then add each of the last layer's groups as children of this group
			if ( ! arg_layers.empty() ) {
				for (const auto &x : common::indices( arg_layers.back().size() ) ) {
					new_group.emplace_back( hierarchy_ref::CLUSTER, x );
				}
			}

			// Add any specified leaf nodes to the group
			for (const auto &leaf_index : arg_leaf_indices) {
				new_group.emplace_back( hierarchy_ref::ENTRY, leaf_index );
			}

			// Add a final with the group
			hierarchy_layer layer;
			layer.emplace_back( std::move( new_group ) );
			arg_layers.push_back( std::move( layer ) );

			// Make a hierarchy from the (reversed) layers
			return make_hierarchy_from_reversed_layers( std::move( arg_layers ) );
		}

	} // namespace clust
} // namespace cath

#endif
