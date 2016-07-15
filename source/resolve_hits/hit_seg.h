/// \file
/// \brief The hit_seg class header

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

#ifndef HIT_SEG_H_INCLUDED
#define HIT_SEG_H_INCLUDED

#include <boost/core/ignore_unused.hpp>
#include <boost/operators.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "exception/invalid_argument_exception.h"
#include "resolve_hits/res_arrow.h"

namespace cath {
	namespace rslv {

		/// \brief Represent a single segment of a sequence hit (ie domain)
		///
		/// This can also be used to represent a fragment
		///
		/// Like all cath-resolve-hits code, this assumes simple residue numbering
		/// and is hence unsuitable for use with raw PDB residue numbers
		class hit_seg final : private boost::equality_comparable<hit_seg> {
		private:
			/// \brief The start boundary of the segment
			res_arrow start;

			/// \brief The stop boundary of the segment
			res_arrow stop;

		public:
			constexpr hit_seg(const res_arrow &,
			                  const res_arrow &);

			constexpr const res_arrow & get_start_arrow() const;
			constexpr const res_arrow & get_stop_arrow() const;

			hit_seg & set_start_arrow(const res_arrow &);
			hit_seg & set_stop_arrow(const res_arrow &);

			static auto get_hit_seg_start_less() {
				return [] (const hit_seg &x, const hit_seg &y) {
					return ( x.get_start_arrow() < y.get_start_arrow() );
				};
			}
		};

		/// \todo Come GCC >= 5 (with relaxed constexpr), add this constexpr sanity_check()
		//
		// template <typename... Ts> constexpr void constexpr_ignore_unused(Ts &&...) {}
		//
		///// \brief Sanity check this hit_seg and throw an exception if a problem is detected
		//constexpr void hit_seg::sanity_check() const {
		//	const int dummy_var = ( start < stop ) ? 0
		//	                                       : throw std::invalid_argument( "Cannot create hit_seg with start residue before the stop residue" );
		//	constexpr_ignore_unused( dummy_var );
		//}

		/// \brief Ctor hit_seg from a start and stop
		constexpr hit_seg::hit_seg(const res_arrow &arg_start, ///< The residue boundary of the segment's start
		                           const res_arrow &arg_stop   ///< The residue boundary of the segment's stop
		                           ) : start( arg_start ),
		                               stop ( arg_stop  ) {
			/// \todo GCC >= 5 (with relaxed constexpr), add this call to the constexpr sanity_check()
			//sanity_check();
		}

		/// \brief Getter for the start boundary
		constexpr const res_arrow & hit_seg::get_start_arrow() const {
			return start;
		}

		/// \brief Getter for the stop boundary
		constexpr const res_arrow & hit_seg::get_stop_arrow() const {
			return stop;
		}

		/// \brief Setter for the start boundary
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr (and remove inline)
		inline hit_seg & hit_seg::set_start_arrow(const res_arrow &arg_start ///< The start boundary to set
		                                          ) {
			start = arg_start;
			return *this;
		}

		/// \brief Setter for the stop boundary
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr (and remove inline)
		inline hit_seg & hit_seg::set_stop_arrow(const res_arrow &arg_stop ///< The stop boundary to set
		                                         ) {
			stop = arg_stop;
			return *this;
		}

		/// \brief Get the start residue index of the specified segment
		///
		/// \relates hit_seg
		constexpr const residx_t & get_start_res_index(const hit_seg &arg_hit_seg ///< The hit_seg to query
		                                               ) {
			return arg_hit_seg.get_start_arrow().res_after();
		}

		/// \brief Get the stop residue index of the specified segment
		///
		/// \relates hit_seg
		constexpr residx_t get_stop_res_index(const hit_seg &arg_hit_seg ///< The hit_seg to query
		                                      ) {
			return arg_hit_seg.get_stop_arrow().res_before();
		}

		/// \brief Get the length of the specified hit_seg
		///
		/// \relates hit_seg
		constexpr size_t get_length(const hit_seg &arg_hit_seg ///< The hit_seg to query
		                            ) {
			return static_cast<size_t>(
				arg_hit_seg.get_stop_arrow ().get_index()
				-
				arg_hit_seg.get_start_arrow().get_index()
			);
		}

		/// \brief Build a hit_seg of the specified start/stop residue indices
		///
		/// \relates hit_seg
		constexpr hit_seg hit_seg_of_res_idcs(const residx_t &arg_start_res_idx, ///< The segment's start residue index
		                                      const residx_t &arg_stop_res_idx   ///< The segment's stop  residue index
		                                      ) {
			return {
				arrow_before_res( arg_start_res_idx ),
				arrow_after_res ( arg_stop_res_idx  )
			};
		}

		/// \brief Build a hit_seg from a pair of start/stop residue indices
		///
		/// \relates hit_seg
		constexpr hit_seg hit_seg_of_res_idx_pair(const residx_residx_pair &arg_res_idx_pair ///< The segments start/stop residue indices
		                                          ) {
			return {
				arrow_before_res( arg_res_idx_pair.first  ),
				arrow_after_res ( arg_res_idx_pair.second )
			};
		}

		/// \brief Return whether the two specified hit_segs are identical
		///
		/// \relates hit_seg
		constexpr bool operator==(const hit_seg &arg_hit_seg_a, ///< The first  hit_seg to compare
		                          const hit_seg &arg_hit_seg_b  ///< The second hit_seg to compare
		                          ) {
			return (
				arg_hit_seg_a.get_start_arrow() == arg_hit_seg_b.get_start_arrow()
				&&
				arg_hit_seg_a.get_stop_arrow()  == arg_hit_seg_b.get_stop_arrow()
			);
		}

		/// \brief Return whether the two specified hit_segs overlap
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \relates hit_seg
		constexpr bool hit_segs_overlap(const hit_seg &arg_hit_seg_a, ///< The first  hit_seg to query
		                                const hit_seg &arg_hit_seg_b  ///< The second hit_seg to query
		                                ) {
			return (
				arg_hit_seg_a.get_start_arrow() <  arg_hit_seg_b.get_stop_arrow()
				&&
				arg_hit_seg_b.get_start_arrow() <  arg_hit_seg_a.get_stop_arrow()
			);
		}

		/// \brief In-place sort the specified segments by their starts
		///
		/// \relates hit_seg
		inline void start_sort_hit_segs(hit_seg_vec &arg_hit_segs ///< The segments to sort
		                                ) {
			boost::range::sort(
				arg_hit_segs,
				hit_seg::get_hit_seg_start_less()
			);
		}

		/// \brief Return a copy of the specified segments, sorted by their starts
		///
		/// \relates hit_seg
		inline hit_seg_vec start_sort_hit_segs_copy(hit_seg_vec arg_hit_segs ///< The segments to sort
		                                            ) {
			start_sort_hit_segs( arg_hit_segs );
			return arg_hit_segs;
		}

		hit_seg_vec make_fragments_of_segments(hit_seg_vec);
		bool segments_are_start_sorted_and_non_overlapping(const hit_seg_vec &);
		hit_seg_vec make_fragments_of_start_sorted_segments(const hit_seg_vec &);
		std::string to_string(const hit_seg &);
		std::ostream & operator<<(std::ostream &,
		                          const hit_seg &);

	}
}

#endif
