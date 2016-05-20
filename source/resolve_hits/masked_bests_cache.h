/// \file
/// \brief The masked_bests_cache class header

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

#ifndef MASKED_BESTS_CACHE_H_INCLUDED
#define MASKED_BESTS_CACHE_H_INCLUDED

// #include <boost/log/trivial.hpp>
// #include <boost/range/algorithm/adjacent_find.hpp>
// #include <boost/range/algorithm/sort.hpp>

// #include "common/c++14/cbegin_cend.h"
// #include "common/chrono/duration_to_seconds_string.h"
#include "resolve_hits/hit.h"
#include "resolve_hits/hit_seg.h"
#include "resolve_hits/scored_arch_proxy.h"

#include <iostream> // ***** TEMPORARY *****
#include <unordered_map>

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief TODOCUMENT
			struct hit_seg_vec_hasher final {

				/// \brief Hash function for hit_vec
				size_t operator()(const hit_seg_vec &arg_hit_seg_vec ///< The hit_seg_vec to hash
				                  ) const {
					size_t seed = 0;
					for (const hit_seg &the_hit : arg_hit_seg_vec) {
						for (const auto &index : { the_hit.get_start_arrow().get_index(),
						                           the_hit.get_stop_arrow ().get_index(), } ) {
							using index_type = std::decay_t< decltype( index ) >;
							seed ^= std::hash< index_type >{}( index  ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
						}
					}
					return seed;
				}
			};

			// This removes zero-length regions
			// /// \brief TODOCUMENT
			// hit_seg_vec get_unmasked_regions_before_arrow(const hit_vec   &arg_hits,
			//                                               const res_arrow &arg_arrow
			//                                               ) {
			// 	const auto hit_segs = get_start_sorted_hit_segs( arg_hits );
			// 	hit_seg_vec results;
			// 	auto prev_stop = start_arrow();
			// 	for (const hit_seg &the_hit_seg : hit_segs) {
			// 		if ( arg_arrow <= the_hit_seg.get_start_arrow()  ) {
			// 			break;
			// 		}
			// 		if ( prev_stop <  the_hit_seg.get_start_arrow() ) {
			// 			results.emplace_back( prev_stop, the_hit_seg.get_start_arrow() );
			// 		}
			// 		prev_stop = the_hit_seg.get_stop_arrow();
			// 	}
			// 	if ( prev_stop < arg_arrow ) {
			// 		results.emplace_back( prev_stop, arg_arrow );
			// 	}
			// 	return results;
			// }

			/// \brief TODOCUMENT
			///
			/// This removes zero-length regions this means that different hit_vecs might
			/// collapse to the same thing
			inline hit_seg_vec get_unmasked_regions_before_arrow(const hit_vec   &arg_hits,
			                                                     const res_arrow &arg_arrow
			                                                     ) {
				const auto hit_segs = get_start_sorted_hit_segs( arg_hits );
				hit_seg_vec results;
				auto prev_stop = start_arrow();
				for (const hit_seg &the_hit_seg : hit_segs) {
					if ( arg_arrow <= the_hit_seg.get_start_arrow()  ) {
						break;
					}
					if ( prev_stop <  the_hit_seg.get_start_arrow() ) {
						results.emplace_back( prev_stop, the_hit_seg.get_start_arrow() );
					}
					prev_stop = the_hit_seg.get_stop_arrow();
				}
				if ( prev_stop < arg_arrow ) {
					results.emplace_back( prev_stop, arg_arrow );
				}
				return results;
			}
		}

		/// \brief TODOCUMENT
		class masked_bests_cache final {
		private:
			/// \brief TODOCUMENT
			std::unordered_map<hit_seg_vec, scored_arch_proxy, detail::hit_seg_vec_hasher> store;

		public:
			const scored_arch_proxy & get_best_for_unmasked(const hit_seg_vec &) const;
			void store_best_for_unmasked(const hit_seg_vec &&,
			                             const scored_arch_proxy &);
		};

		/// \brief TODOCUMENT
		inline const scored_arch_proxy & masked_bests_cache::get_best_for_unmasked(const hit_seg_vec &arg_unmasked ///< TODOCUMENT
		                                                                           ) const {
			// std::cerr << "Reading :";
			// for (const hit_seg &the_seg : arg_unmasked) {
			// 	std::cerr << " " << get_start_res_index( the_seg ) << "-" << get_stop_res_index( the_seg );
			// }
			// std::cerr << "\n";
			// cerr << "Writing : " <<  << "\n";
			return store.at( arg_unmasked );
		}

		/// \brief TODOCUMENT
		inline void masked_bests_cache::store_best_for_unmasked(const hit_seg_vec       &&arg_unmasked,              ///< TODOCUMENT
		                                                        const scored_arch_proxy  &arg_best_scored_arch_proxy ///< TODOCUMENT
		                                                        ) {
			// std::cerr << "Writing :";
			// for (const hit_seg &the_seg : arg_unmasked) {
			// 	std::cerr << " " << get_start_res_index( the_seg ) << "-" << get_stop_res_index( the_seg );
			// }
			// std::cerr << "\n";
			// cerr << "Writing : " <<  << "\n";
			store.emplace(
				std::move( arg_unmasked ),
				std::move( arg_best_scored_arch_proxy )
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates masked_bests_cache
		inline const scored_arch_proxy & get_best_for_masks_up_to_arrow(const masked_bests_cache &arg_masked_bests_cache, ///< TODOCUMENT
		                                                                const hit_vec            &arg_mask_hits,          ///< TODOCUMENT
		                                                                const res_arrow          &arg_stop_arrow          ///< TODOCUMENT
		                                                                ) {
			return arg_masked_bests_cache.get_best_for_unmasked(
				detail::get_unmasked_regions_before_arrow( arg_mask_hits, arg_stop_arrow )
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates masked_bests_cache
		///
		///     Store       : score
		///     Store       : arch - discontigs
		///     Check       : nothing in arch - discontigs should overrun current_arrow
		///     Check       : nothing previously stored there
		///     This needs to remove mask from arch
		inline void store_best_for_masks_up_to_arrow(masked_bests_cache      &arg_masked_bests_cache,     ///< TODOCUMENT
		                                             const scored_arch_proxy &arg_best_scored_arch_proxy, ///< TODOCUMENT
		                                             const hit_vec           &arg_mask_hits,              ///< TODOCUMENT
		                                             const res_arrow         &arg_stop_arrow              ///< TODOCUMENT
		                                             ) {
			// auto &&masked_scored_arch =  arg_best_scored_arch_proxy - arg_mask_hits;

// #ifndef NDEBUG
			
// #endif

// #ifndef NDEBUG
// 			for (const hit &the_hit : masked_scored_arch.get_arch() ) {
// 				if ( the_hit.get_stop_arrow() > arg_stop_arrow ) {
// 					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot store masked arch that overruns the stop arrow"));
// 				}
// 			}
// #endif

			arg_masked_bests_cache.store_best_for_unmasked(
				detail::get_unmasked_regions_before_arrow( arg_mask_hits, arg_stop_arrow ),
				// std::forward< decltype( masked_scored_arch ) >( masked_scored_arch )
				arg_best_scored_arch_proxy
			);
		}

	}
}

#endif
