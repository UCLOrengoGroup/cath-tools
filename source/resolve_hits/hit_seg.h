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

#include "exception/invalid_argument_exception.h"
#include "resolve_hits/res_arrow.h"

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		///
		class hit_seg final {
		private:
			/// \brief TODOCUMENT
			res_arrow start;

			/// \brief TODOCUMENT
			res_arrow stop;

			//constexpr void sanity_check() const;

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

		/// \brief TODOCUMENT
		template <typename... Ts>
		constexpr void constexpr_ignore_unused(Ts &&...) {}

		/// \todo Come GCC >= 5 (with relaxed constexpr), add this constexpr sanity_check()
		//
		///// \brief TODOCUMENT
		//constexpr void hit_seg::sanity_check() const {
		//	const int dummy_var = ( start < stop ) ? 0
		//	                                       : throw std::invalid_argument( "Cannot create hit_seg with start residue before the stop residue" );
		//	constexpr_ignore_unused( dummy_var );
		//}

		/// \brief TODOCUMENT
		constexpr hit_seg::hit_seg(const res_arrow &arg_start, ///< TODOCUMENT
		                           const res_arrow &arg_stop   ///< TODOCUMENT
		                           ) : start( arg_start ),
		                               stop ( arg_stop  ) {
			/// \todo GCC >= 5 (with relaxed constexpr), add this call to the constexpr sanity_check()
			//sanity_check();
		}

		/// \brief TODOCUMENT
		constexpr const res_arrow & hit_seg::get_start_arrow() const {
			return start;
		}

		/// \brief TODOCUMENT
		constexpr const res_arrow & hit_seg::get_stop_arrow() const {
			return stop;
		}

		/// \brief TODOCUMENT
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr (and remove inline)
		inline hit_seg & hit_seg::set_start_arrow(const res_arrow &arg_start ///< TODOCUMENT
		                                          ) {
			start = arg_start;
			return *this;
		}

		/// \brief TODOCUMENT
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr (and remove inline)
		inline hit_seg & hit_seg::set_stop_arrow(const res_arrow &arg_stop ///< TODOCUMENT
		                                         ) {
			stop = arg_stop;
			return *this;
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit_seg
		constexpr const residx_t & get_start_res_index(const hit_seg &arg_hit_seg ///< TODOCUMENT
		                                               ) {
			return arg_hit_seg.get_start_arrow().res_after();
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit_seg
		constexpr residx_t get_stop_res_index(const hit_seg &arg_hit_seg ///< TODOCUMENT
		                                      ) {
			return arg_hit_seg.get_stop_arrow().res_before();
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit_seg
		constexpr size_t get_length(const hit_seg &arg_hit_seg ///< TODOCUMENT
		                            ) {
			return static_cast<size_t>(
				arg_hit_seg.get_stop_arrow ().get_index()
				-
				arg_hit_seg.get_start_arrow().get_index()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit_seg
		constexpr hit_seg hit_seg_of_res_idcs(const residx_t &arg_start_res_idx, ///< TODOCUMENT
		                                      const residx_t &arg_stop_res_idx   ///< TODOCUMENT
		                                      ) {
			return {
				arrow_before_res( arg_start_res_idx ),
				arrow_after_res ( arg_stop_res_idx  )
			};
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit_seg
		constexpr hit_seg hit_seg_of_res_idx_pair(const residx_residx_pair &arg_res_idx_pair ///< TODOCUMENT
		                                          ) {
			return {
				arrow_before_res( arg_res_idx_pair.first  ),
				arrow_after_res ( arg_res_idx_pair.second )
			};
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit_seg
		constexpr bool operator==(const hit_seg &arg_hit_seg_a, ///< TODOCUMENT
		                          const hit_seg &arg_hit_seg_b  ///< TODOCUMENT
		                          ) {
			return (
				arg_hit_seg_a.get_start_arrow() == arg_hit_seg_b.get_start_arrow()
				&&
				arg_hit_seg_a.get_stop_arrow()  == arg_hit_seg_b.get_stop_arrow()
			);
		}

		/// \brief TODOCUMENT
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \relates hit_seg
		constexpr bool hit_segs_overlap(const hit_seg &arg_hit_seg_a, ///< TODOCUMENT
		                                const hit_seg &arg_hit_seg_b  ///< TODOCUMENT
		                                ) {
			return (
				arg_hit_seg_a.get_start_arrow() <  arg_hit_seg_b.get_stop_arrow()
				&&
				arg_hit_seg_b.get_start_arrow() <  arg_hit_seg_a.get_stop_arrow()
			);
		}

		void start_sort_hit_segs(hit_seg_vec &);
		hit_seg_vec start_sort_hit_segs_copy(hit_seg_vec);
		hit_seg_vec make_fragments_of_segments(hit_seg_vec);
		bool segments_are_start_sorted_and_non_overlapping(const hit_seg_vec &);
		hit_seg_vec make_fragments_of_start_sorted_segments(const hit_seg_vec &);
		std::string to_string(const hit_seg &);
		std::ostream & operator<<(std::ostream &,
		                          const hit_seg &);

		

	}
}

#endif
