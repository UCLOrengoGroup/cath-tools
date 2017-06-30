/// \file
/// \brief The full_hit_prune_builder class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_DETAIL_FULL_HIT_PRUNE_BUILDER_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_DETAIL_FULL_HIT_PRUNE_BUILDER_H

#include <boost/range/algorithm/equal.hpp>

#include "common/boost_addenda/tribool/tribool.hpp"
#include "resolve_hits/first_hit_is_better.hpp"
#include "resolve_hits/first_hit_is_better.hpp"
#include "resolve_hits/full_hit.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"
#include "resolve_hits/seg_dupl_hit_policy.hpp"

#include <unordered_map>

namespace cath {

	namespace common {

			/// \brief Attempt to do a *reasonably* efficient "move" of a deque into a vector
			///
			/// Note: this doesn't offer the strong exception safety guarantee; it just
			///       moves elements across without checking if T's move operations are noexcept
			template <typename T>
			std::vector<T> make_vector_of_rvalue_deque(std::deque<T> &&arg_deque ///< An rvalue deque<T> to "move" into a vector<T>
			                                           ) {
				std::vector<T> results;
				results.reserve( arg_deque.size() );
				for (const T &element : arg_deque) {
					// Better to push_back() or emplace_back() if given an rvalue T?
					results.push_back( std::move( element ) );
				}
				arg_deque.clear();
				return results;
			}

	}

	namespace rslv {
		namespace detail {

			/// \brief Hasher for reference_wrapper<full_hit>
			///        which just hashes the full_hit's segments
			struct full_hit_ref_segs_hasher final {
				/// \brief The function operator that performs the hash on the full_hit's segments
				size_t operator()(const std::reference_wrapper<full_hit> &arg_value
				                  ) const {
					return calc_hash( arg_value.get().get_segments() );
				}
			};

			/// \brief Equality function object for reference_wrapper<full_hit>
			///        which just compares the full_hit's segments
			struct full_hit_ref_segs_equal_to final {
				/// \brief The function operator that performs compares two full_hits' segments
				size_t operator()(const std::reference_wrapper<full_hit> &arg_lhs, ///< The first  reference_wrapper<full_hit> to compare
				                  const std::reference_wrapper<full_hit> &arg_rhs  ///< The second reference_wrapper<full_hit> to compare
				                  ) const {
					return boost::range::equal(
						arg_lhs.get().get_segments(),
						arg_rhs.get().get_segments()
					);
				}
			};

			/// \brief Build a full_hit_vec of hits whilst pruning hits that have a better hit (or not worse) with *identical* segments
			///
			/// This does a better job than calc_hit_prune_builder
			class full_hit_prune_builder final {
			private:
				/// \brief Whether the strictly-worse hits should be preserved or pruned
				seg_dupl_hit_policy policy;

				/// \brief The full_hit_vec to be built up
				///
				/// Using deque rather than vector to ensure that the references stored in the unordered_map
				/// below aren't invalidated as this grows
				std::deque<full_hit> hits;

				/// \brief Type alias for an unordered_map from reference_wrapper<full_hit> to size_t
				///        that uses full_hit_ref_segs_hasher and full_hit_ref_segs_equal_to to make
				///        hashing and equality-comparison be performed on the full_hits' segments
				using full_hit_to_index_uomap = std::unordered_map<
					std::reference_wrapper<full_hit>,
					size_t,
					full_hit_ref_segs_hasher,
					full_hit_ref_segs_equal_to
				>;

				/// \brief A lookup from a full_hit_ref to an index the index of an example full_hit with those segments
				full_hit_to_index_uomap index_of_full_hit_ref;

			public:
				/// \brief Ctor to populate require_strictly_worse_hits
				inline explicit full_hit_prune_builder(const seg_dupl_hit_policy &arg_policy ///< Whether the strictly-worse hits should be preserved or pruned
				                                       ) : policy { arg_policy } {
				}

				/// \brief Reserve space for the specified number of hits
				inline void reserve(const size_t &arg_capacity ///< The number of hits for which space should be reserved
				                    ) {
					index_of_full_hit_ref.reserve( arg_capacity );
				}


				/// \brief Return the number of hits stored thus far (after any pruning)
				inline size_t size() const {
					return hits.size();
				}

				/// \brief Whether there haven't been any hits built so far
				inline bool empty() const {
					return hits.empty();
				}

				/// \brief Add a hit (or skip it if it's found to be worse than existing hit)
				inline void add_hit(full_hit arg_full_hit ///< The hit to be added
				                    ) {
					if ( policy == seg_dupl_hit_policy::PRESERVE ) {
						hits.push_back( std::move( arg_full_hit ) );
						return;
					}
					const auto itr = index_of_full_hit_ref.find( arg_full_hit );
					if ( itr == common::cend( index_of_full_hit_ref ) ) {
						const auto hits_size_before = hits.size();
						hits.push_back( std::move( arg_full_hit ) );
						index_of_full_hit_ref.emplace( hits.back(), hits_size_before );
					}
					else {
						const auto &comp_hit = hits[ itr->second ];
						const auto  result   = first_hit_is_better( arg_full_hit, comp_hit );
						if ( common::is_true( result ) ) {
							hits[ itr->second ] = std::move( arg_full_hit );
						}
						else {
						}
					}
				}

				/// \brief Get the built hits
				inline full_hit_list get_built_hits() {
					full_hit_list result{ common::make_vector_of_rvalue_deque( std::move( hits ) ) };
					index_of_full_hit_ref.clear();
					return result;
				}
			};

		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif
