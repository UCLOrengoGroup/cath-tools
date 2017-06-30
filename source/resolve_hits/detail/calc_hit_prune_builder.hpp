/// \file
/// \brief The calc_hit_prune_builder class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_DETAIL_CALC_HIT_PRUNE_BUILDER_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_DETAIL_CALC_HIT_PRUNE_BUILDER_H

#include "common/boost_addenda/tribool/tribool.hpp"
#include "resolve_hits/calc_hit.hpp"
#include "resolve_hits/first_hit_is_better.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"
#include "resolve_hits/seg_dupl_hit_policy.hpp"

#include <unordered_map>

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Build a calc_hit_vec of hits whilst pruning some (but not reliably all)
			/// of the hits that have a better hit with identical segments
			///
			/// This works by storing an unordered_map of the hash of the segments signature
			/// to an index of an example calc_hit with that signature
			///
			/// It's a bit unusual to store a hash value (which the unordered_map will hash again).
			/// This means different signatures may collide but that should be relatively rare and is tolerable.
			/// And it allows the unordered_map to be much faster because it's just storing pairs of ints rather
			/// than having to store a copy of each hit's seq_seg_run (including a vector).
			///
			/// full_hit_prune_builder does a better job using a reference_wrapper
			class calc_hit_prune_builder final {
			private:
				/// \brief Whether the strictly-worse hits should be preserved or pruned
				seg_dupl_hit_policy policy;

				/// \brief The calc_hit_vec to be built up
				calc_hit_vec hits;

				/// \brief A lookup from a hash of the segments signature to an index of an example calc_hit
				///        with that segments signature
				std::unordered_map<size_t, size_t> index_of_signature_hash;

			public:
				/// \brief Ctor from require_strictly_worse_hits
				inline explicit calc_hit_prune_builder(const seg_dupl_hit_policy &arg_policy ///< Whether the strictly-worse hits should be preserved or pruned
				                                       ) : policy { arg_policy } {
				}

				/// \brief Reserve space for the specified number of hits
				inline void reserve(const size_t &arg_capacity ///< The number of hits for which space should be reserved
				                    ) {
					hits.reserve( arg_capacity );
					index_of_signature_hash.reserve( arg_capacity );
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
				inline void add_hit(calc_hit             arg_calc_hit, ///< The hit to be added
				                    const full_hit_list &arg_full_hits ///< The full_hit_list from which the calc_hit comes (for allowing comparison with first_hit_is_better() )
				                    ) {
					if ( policy == seg_dupl_hit_policy::PRESERVE ) {
						hits.push_back( std::move( arg_calc_hit ) );
						return;
					}

					const size_t hash_value = calc_hash( arg_calc_hit.get_segments() );
					const auto   itr        = index_of_signature_hash.find( hash_value );
					if ( itr == common::cend( index_of_signature_hash ) ) {
						index_of_signature_hash.emplace( hash_value, hits.size() );
						hits.push_back( std::move( arg_calc_hit ) );
					}
					else {
						const auto &comp_hit = hits[ itr->second ];
						const auto  result   = first_hit_is_better( arg_calc_hit, comp_hit, arg_full_hits );
						if ( common::is_true( result ) ) {
							hits[ itr->second ] = std::move( arg_calc_hit );
						}
						else if ( boost::logic::indeterminate( result ) ) {
							hits.push_back( std::move( arg_calc_hit ) );
						}
					}
				}

				/// \brief Get a non-const reference to the built hits
				inline calc_hit_vec & get_built_hits() {
					return hits;
				}
			};

		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif
