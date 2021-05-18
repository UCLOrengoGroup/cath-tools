/// \file
/// \brief The links class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINKS_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINKS_HPP

#include <algorithm>
#include <filesystem>
#include <functional>

#include "cath/clustagglom/link.hpp"
#include "cath/clustagglom/link_dirn.hpp"
#include "cath/clustagglom/link_list.hpp"
#include "cath/common/algorithm/contains.hpp"
#include "cath/common/cpp17/apply.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath { namespace common { class id_of_str_bidirnl; } }

namespace cath {
	namespace clust {

		namespace detail {

			/// \brief The number of previously seen scores for some particular entry
			enum class num_prev_scores : bool {
				ZERO, ///< Previously seen zero times before
				ONE,  ///< Previously seen one  time  before
			};

			/// \brief Type alias for a pair of num_prev_scores and strength
			using num_prev_scores_strength_pair     = std::pair<num_prev_scores, strength>;

			/// \brief Type alias for a vector of num_prev_scores_strength_pair values
			using num_prev_scores_strength_pair_vec = std::vector<num_prev_scores_strength_pair>;

		} // namespace detail


		/// \brief A sparse matrix of links between items
		///
		/// TODO: Consider preventing self-linking?
		class links final {
		private:
			/// \brief The links
			link_list_vec the_link_lists;

			/// \brief A data structure that's used for merging and is kept
			///        around between merges to avoid it being reallocated
			detail::num_prev_scores_strength_pair_vec update_dists_data;

		public:
			/// \brief A const_iterator type alias as part of making this a range over links
			using const_iterator = link_list_vec_citr;

			links() = default;

			[[nodiscard]] bool   empty() const;
			[[nodiscard]] size_t size() const;

			links & add_link_symmetrically(const size_t &,
			                               const size_t &,
			                               const strength &);

			template <typename Fn>
			links & merge(const item_idx &,
			              const item_idx &,
			              const item_idx &,
			              Fn &&);

			/// TODOSOON: ***** SHOULD THIS BE HERE???? *******
			/// TODOSOON: ***** SHOULD THIS BE HERE???? *******
			/// TODOSOON: ***** SHOULD THIS BE HERE???? *******
			/// TODOSOON: ***** SHOULD THIS BE HERE???? *******
			/// TODOSOON: ***** SHOULD THIS BE HERE???? *******
			link_list & operator[](const size_t &);
			const link_list & operator[](const size_t &) const;

			[[nodiscard]] const_iterator begin() const;
			[[nodiscard]] const_iterator end() const;
		};

		/// \brief Return whether this collection of links is empty
		inline bool links::empty() const {
			return the_link_lists.empty();
		}

		/// \brief Return the number of lists of links
		inline size_t links::size() const {
			return the_link_lists.size();
		}

		/// \brief Non-const-overload of standard subscript operator
		inline link_list & links::operator[](const size_t &prm_index ///< The index of the list of links to access
		                                     ) {
			return the_link_lists[ prm_index ];
		}

		/// \brief Const-overload of standard subscript operator
		inline const link_list & links::operator[](const size_t &prm_index ///< The index of the list of links to access
		                                           ) const {
			return the_link_lists[ prm_index ];
		}

		/// \brief Symmetrically add a link of the specified dissimilarity between the specified items
		inline links & links::add_link_symmetrically(const size_t   &prm_index_a,      ///< The first  item to be linked
		                                             const size_t   &prm_index_b,      ///< The second item to be linked
		                                             const strength &prm_dissimilarity ///< The dissimilarity between the two items
		                                             ) {
			if ( prm_index_a == prm_index_b ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot add a self-link to links"));
			}
			if ( the_link_lists.size() < std::max( prm_index_a, prm_index_b ) + 1 ) {
				the_link_lists.resize( std::max( prm_index_a, prm_index_b ) + 1 );
			}
			the_link_lists[ prm_index_a ].emplace_back( prm_index_b, prm_dissimilarity );
			the_link_lists[ prm_index_b ].emplace_back( prm_index_a, prm_dissimilarity );
			return *this;
		}

		/// \brief Merge the two specified clusters into a new cluster of the specified label
		///        using the specified function to determine the dissimilarities between the
		///        new cluster and others based on the previous cluster's dissimilarities to
		///        the cluster in question
		template <typename Fn>
		links & links::merge(const item_idx &prm_a,         ///< The first  cluster to merge
		                     const item_idx &prm_b,         ///< The second cluster to merge
		                     const item_idx &prm_new_label, ///< The label of the new cluster to form
		                     Fn            &&prm_fn         ///< The function from the index of the target cluster and the two dissimilarities with that cluster from the mergee clusters. The target cluster may no-longer exist. May return none for no link.
		                                                    ///< !!!!! At the moment, this function is only called where both links exists - OK for complete-linkage but not others
		                     ) {
			// Ensure that there is enough space in the_link_lists
			the_link_lists.resize( prm_new_label + 1 );

			// Update distances...
			update_dists_data.assign(
				prm_new_label + 1,
				detail::num_prev_scores_strength_pair{ detail::num_prev_scores::ZERO, 0.0 }
			);
			// Copy all the dissimilarities from the first links_list into update_dists_data
			for (const link &x : the_link_lists[ prm_a ] ) {
				update_dists_data[ x.node ] = { detail::num_prev_scores::ONE, x.dissim };
			}

			// Swap the already-allocated memory from the_link_lists[ prm_a ] with the_link_lists[ prm_new_label ] and clear it
			std::swap( the_link_lists[ prm_a ], the_link_lists[ prm_new_label ] );
			the_link_lists[ prm_new_label ].clear();

			// Update all the links between the new cluster and each target cluster
			for (const link &x : the_link_lists[ prm_b ] ) {
				const auto &update_dists_val = update_dists_data[ x.node ];
				if ( update_dists_val.first == detail::num_prev_scores::ONE ) {
					const strength_opt dissim_opt = ::std::invoke( prm_fn, x.node, update_dists_val.second, x.dissim );
					if ( dissim_opt ) {
						the_link_lists[ prm_new_label ].emplace_back( x.node,        *dissim_opt );
						the_link_lists[ x.node        ].emplace_back( prm_new_label, *dissim_opt );
					}
				}
			}

			// Free up memory that's no longer required
			the_link_lists[ prm_a ].clear();
			the_link_lists[ prm_a ].shrink_to_fit();
			the_link_lists[ prm_b ].clear();
			the_link_lists[ prm_b ].shrink_to_fit();

			return *this;
		}

		/// \brief Standard begin() method to allow iteration over the lists of links
		inline auto links::begin() const -> const_iterator {
			return ::std::cbegin( the_link_lists );
		}

		/// \brief Standard end() method to allow iteration over the lists of links
		inline auto links::end() const -> const_iterator {
			return ::std::cend  ( the_link_lists );
		}

		/// \brief Symmetrically add the specified link tuple to the specified links
		///
		/// This forwards the parts of the tuple to the add_link_symmetrically() method of links
		inline void add_link_symmetrically(links                        &prm_links, ///< The links to which the new link should be added
		                                   const item_item_strength_tpl &prm_link   ///< The link specified as a tuple
		                                   ) {
			// Pass the parts of the tuple as arguments to the add_link_symmetrically() member function of prm_links
			common::apply(
				[&] (auto...args) {
					prm_links.add_link_symmetrically(
						std::forward< decltype( args ) >( args )...
					);
				},
				prm_link
			);
		}

		links make_links(const item_item_strength_tpl_vec &);

		size_size_pair_vec get_spanning_tree_of_subset(const links &,
		                                               const size_set &);

		std::string to_string(const links &);

		void write_ordered_links(std::ostream &,
		                         const links &,
		                         const common::id_of_str_bidirnl &,
		                         const size_vec &);

		void write_ordered_links(const ::std::filesystem::path &,
		                         const links &,
		                         const common::id_of_str_bidirnl &,
		                         const size_vec &);

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINKS_HPP
