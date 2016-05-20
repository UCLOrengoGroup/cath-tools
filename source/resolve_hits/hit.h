/// \file
/// \brief The hit class header

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

#ifndef HIT_H_INCLUDED
#define HIT_H_INCLUDED

#include <boost/range/irange.hpp>

#include "common/algorithm/transform_build.h"
#include "common/c++14/cbegin_cend.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "resolve_hits/hit_seg.h"
#include "resolve_hits/res_arrow.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

using namespace cath::common::literals;

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		///
		/// Invariants:
		///  * all segments must have the same residue_locating
		///     (ie whether they locate their residues by names and/or indices)
		///
		/// \todo Change internals to use arrows rather than residx_t values and then
		//        make get_*_res_index* methods into non-member functions
		///      
		class hit final {
		private:
			/// \brief TODOCUMENT
			res_arrow start_arrow;

			/// \brief TODOCUMENT
			res_arrow stop_arrow;

			/// \brief TODOCUMENT
			resscr_t score;

			/// \brief TODOCUMENT
			std::string label;

			/// \brief TODOCUMENT
			hit_seg_vec fragments;

			void sanity_check() const;

			// using const_iterator = residx_residx_pair_vec::const_iterator;

			// const_iterator begin() const;
			// const_iterator end() const;

		public:
			hit(const res_arrow &,
			    const res_arrow &,
			    const resscr_t &,
			    const std::string &);

			hit(const hit_seg_vec &,
			    const resscr_t &,
			    const std::string &);

			bool is_discontig() const;
			size_t get_num_segments() const; // return fragments.size() + 1;
			const res_arrow & get_start_arrow_of_segment(const size_t &) const;
			const res_arrow & get_stop_arrow_of_segment(const size_t &) const;

			const res_arrow   & get_start_arrow() const;
			const res_arrow   & get_stop_arrow()  const;
			const resscr_t    & get_score()       const;
			const std::string & get_label()       const;

			static auto get_hit_start_less() {
				return [] (const hit &x, const hit &y) {
					return ( x.get_start_arrow() < y.get_start_arrow() );
				};
			}

			static auto get_hit_stop_less() {
				return [] (const hit &x, const hit &y) {
					return ( x.get_stop_arrow() < y.get_stop_arrow() );
				};
			}
		};

		using hit_vec = std::vector<hit>;

		/// \brief TODOCUMENT
		inline hit_vec operator+(hit_vec    arg_hit_vec, ///< TODOCUMENT
		                         const hit &arg_hit      ///< TODOCUMENT
		                         ) {
			arg_hit_vec.push_back( arg_hit );
			return arg_hit_vec;
		}

		enum class hit_output_format {
			CLASS,
			JON
		};

		std::string get_segments_string(const hit &);
		std::string to_string(const hit &,
		                      const hit_output_format & = hit_output_format::CLASS);
		std::ostream & operator<<(std::ostream &,
		                          const hit &);
		bool operator==(const hit &,
		                const hit &);
		bool any_interaction(const hit &,
		                     const hit &);

		/// \brief TODOCUMENT
		inline void hit::sanity_check() const {
			if ( stop_arrow < start_arrow ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("start index must not be greater than the stop index"));
			}
			if ( score <= 0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("score must be strictly greater than 0 else the algorithm doesn't work (because there's no way to no how to trade scores off against empty space)"));
			}
		}

		// /// \brief TODOCUMENT
		// inline auto hit::begin() const -> const_iterator {
		// 	return common::cbegin( fragments );
		// }

		// /// \brief TODOCUMENT
		// inline auto hit::end() const -> const_iterator {
		// 	return common::cend  ( fragments );
		// }

		/// \brief Ctor for contiguous hit
		inline hit::hit(const res_arrow   &arg_start_arrow, ///< TODOCUMENT
		                const res_arrow   &arg_stop_arrow,  ///< TODOCUMENT
		                const resscr_t    &arg_score,       ///< TODOCUMENT
		                const std::string &arg_label        ///< TODOCUMENT
		                ) : start_arrow ( arg_start_arrow ),
		                    stop_arrow  ( arg_stop_arrow  ),
		                    score       ( arg_score       ),
		                    label       ( arg_label       ) {
			sanity_check();
		}

		inline hit::hit(const hit_seg_vec &arg_segments,
		                const resscr_t    &arg_score,
		                const std::string &arg_label
		                ) : start_arrow    ( arg_segments.front().get_start_arrow()      ),
		                    stop_arrow     ( arg_segments.back ().get_stop_arrow ()      ),
		                    score          ( arg_score                                   ),
		                    label          ( arg_label                                   ),
		                    fragments      ( make_fragments_of_segments( arg_segments )  ) {
			sanity_check();
		}

		/// \brief TODOCUMENT
		inline bool hit::is_discontig() const {
			return ! fragments.empty();
		}

		/// \brief TODOCUMENT
		inline size_t hit::get_num_segments() const {
			return fragments.size() + 1;
		}

		/// \brief TODOCUMENT
		inline const res_arrow & hit::get_start_arrow_of_segment(const size_t &arg_segment_index ///< TODOCUMENT
		                                                         ) const {
			return ( arg_segment_index > 0                ) ? fragments[ arg_segment_index - 1 ].get_stop_arrow()
			                                                : start_arrow;
		}

		/// \brief TODOCUMENT
		inline const res_arrow & hit::get_stop_arrow_of_segment(const size_t &arg_segment_index ///< TODOCUMENT
		                                                        ) const {
			return ( arg_segment_index < fragments.size() ) ? fragments[ arg_segment_index     ].get_start_arrow()
			                                                : stop_arrow;
		}
		
		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline const hit_seg get_hit_seg_of_segment(const hit    &arg_hit,    ///< TODOCUMENT
		                                            const size_t &arg_seg_idx ///< TODOCUMENT
		                                            ) {
			return {
				arg_hit.get_start_arrow_of_segment( arg_seg_idx ),
				arg_hit.get_stop_arrow_of_segment ( arg_seg_idx )
			};
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline const hit_seg_vec get_hit_segs(const hit &arg_hit ///< TODOCUMENT
		                                      ) {
			return common::transform_build<hit_seg_vec>(
				boost::irange( 0_z, arg_hit.get_num_segments() ),
				[&] (const size_t &x) {
					return get_hit_seg_of_segment( arg_hit, x );
				}
			);
		}

		/// \brief TODOCUMENT
		template <typename T, typename U>
		void append(std::vector<T> &arg_vec, ///< TODOCUMENT
		            const U        &arg_rng  ///< TODOCUMENT
		            ) {
			arg_vec.insert(
				common::cend  ( arg_vec ),
				common::cbegin( arg_rng ),
				common::cend  ( arg_rng )
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline const hit_seg_vec get_hit_segs(const hit_vec &arg_hit_vec ///< TODOCUMENT
		                                      ) {
			hit_seg_vec results;
			for (const hit &the_hit : arg_hit_vec) {
				append( results, get_hit_segs( the_hit ) );
			}
			return results;
		}

		/// \brief TODOCUMENT
		inline void start_sort_hit_segs(hit_seg_vec &arg_hit_segs ///< TODOCUMENT
		                                ) {
			boost::range::sort(
				arg_hit_segs,
				[] (const hit_seg &x, const hit_seg &y) { return ( x.get_start_arrow() < y.get_start_arrow() ); }
			);
		}

		/// \brief TODOCUMENT
		inline hit_seg_vec start_sort_hit_segs_copy(hit_seg_vec arg_hit_segs ///< TODOCUMENT
		                                            ) {
			start_sort_hit_segs( arg_hit_segs );
			return arg_hit_segs;
		}

		/// \brief TODOCUMENT
		inline hit_seg_vec get_start_sorted_hit_segs(const hit_vec &arg_hit_vec ///< TODOCUMENT
		                                             ) {
			return start_sort_hit_segs_copy( get_hit_segs( arg_hit_vec ) );
		}

		/// \brief TODOCUMENT
		inline const res_arrow & hit::get_start_arrow() const {
			return start_arrow;
		}

		/// \brief TODOCUMENT
		inline const res_arrow & hit::get_stop_arrow()  const {
			return stop_arrow;
		}

		/// \brief TODOCUMENT
		inline const resscr_t & hit::get_score() const {
			return score;
		}

		/// \brief TODOCUMENT
		inline const std::string & hit::get_label() const {
			return label;
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline const residx_t & get_start_res_index_of_segment(const hit    &arg_hit,          ///< TODOCUMENT
		                                                       const size_t &arg_segment_index ///< TODOCUMENT
		                                                       ) {
			return arg_hit.get_start_arrow_of_segment( arg_segment_index ).res_after();
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline residx_t get_stop_res_index_of_segment(const hit    &arg_hit,          ///< TODOCUMENT
		                                              const size_t &arg_segment_index ///< TODOCUMENT
		                                              ) {
			return arg_hit.get_stop_arrow_of_segment( arg_segment_index ).res_before();
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline const residx_t & get_start_res_index(const hit &arg_hit ///< TODOCUMENT
		                                            ) {
			return arg_hit.get_start_arrow().res_after();
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline residx_t get_stop_res_index(const hit &arg_hit ///< TODOCUMENT
		                                   ) {
			return arg_hit.get_stop_arrow().res_before();
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline res_arrow get_stop_of_first_segment(const hit &arg_hit ///< TODOCUMENT
		                                           ) {
			if ( ! arg_hit.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_stop_of_first_segment of contiguous hit"));
			}
			return arg_hit.get_stop_arrow_of_segment( 0 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline res_arrow get_start_of_last_segment(const hit &arg_hit ///< TODOCUMENT
		                                           ) {
			if ( ! arg_hit.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_start_of_last_segment of contiguous hit"));
			}
			return arg_hit.get_start_arrow_of_segment( arg_hit.get_num_segments() - 1 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline hit make_hit_from_res_indices(const residx_t    &arg_start_res_idx, ///< TODOCUMENT
		                                     const residx_t    &arg_stop_res_idx,  ///< TODOCUMENT
		                                     const resscr_t    &arg_score,         ///< TODOCUMENT
		                                     const std::string &arg_label          ///< TODOCUMENT
		                                     ) {
			return {
				arrow_before_res( arg_start_res_idx ),
				arrow_after_res ( arg_stop_res_idx  ),
				arg_score,
				arg_label
			};
		}

		/// \brief TODOCUMENT
		///
		/// \relates hit
		inline hit make_hit_from_res_indices(const residx_residx_pair_vec &arg_residue_index_segments, ///< TODOCUMENT
		                                     const resscr_t               &arg_score,                  ///< TODOCUMENT
		                                     const std::string            &arg_label                   ///< TODOCUMENT
		                                     ) {
			return {
				common::transform_build<hit_seg_vec>(
					arg_residue_index_segments,
					hit_seg_of_res_idx_pair
				),
				arg_score,
				arg_label
			};
		}

		/// \brief Return whether the either of the two specified hits overlaps, interleaves or straddles the other
		///
		/// \relates hit
		inline bool any_interaction(const hit &arg_hit_a, ///< TODOCUMENT
		                            const hit &arg_hit_b  ///< TODOCUMENT
		                            ) {
			return (
				arg_hit_a.get_start_arrow() < arg_hit_b.get_stop_arrow()
				&&
				arg_hit_b.get_start_arrow() < arg_hit_a.get_stop_arrow()
			);
		}
		
		/// \brief TODOCUMENT
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \todo This can be made more efficient by iterating over the two segment
		///       lists simultaneously
		///
		/// \relates hit
		inline bool hits_overlap(const hit &arg_hit_a, ///< TODOCUMENT
		                         const hit &arg_hit_b  ///< TODOCUMENT
		                         ) {
			if ( ! any_interaction( arg_hit_a, arg_hit_b ) ) {
				return false;
			}
			for (const auto &seg_ctr_a : boost::irange( 0_z, arg_hit_a.get_num_segments() ) ) {
				for (const auto &seg_ctr_b : boost::irange( 0_z, arg_hit_b.get_num_segments() ) ) {
					const bool seg_overlap = hit_segs_overlap(
						get_hit_seg_of_segment( arg_hit_a, seg_ctr_a ),
						get_hit_seg_of_segment( arg_hit_b, seg_ctr_b )
					);
					if ( seg_overlap ) {
						return true;
					}
				}
			}
			return false;
		}

		/// \brief TODOCUMENT
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \relates hit
		inline bool hit_overlaps_with_any_of_hits(const hit     &arg_hit_lhs, ///< TODOCUMENT
		                                          const hit_vec &arg_hit_vec  ///< TODOCUMENT
		                                          ) {
			for (const hit &hit_rhs : arg_hit_vec) {
				if ( hits_overlap( arg_hit_lhs, hit_rhs ) ) {
					return true;
				}
			}
			return false;
		}

		inline bool second_right_intersperses_first(const hit &arg_hit_a, ///< TODOCUMENT
		                                            const hit &arg_hit_b  ///< TODOCUMENT
		                                            ) {
			if ( ! arg_hit_a.is_discontig() || ! arg_hit_b.is_discontig() ) {
				return false;
			}
			const bool ends_are_ok = (
				arg_hit_a.get_start_arrow() < arg_hit_b.get_start_arrow()
				&&
				arg_hit_a.get_stop_arrow()  < arg_hit_b.get_stop_arrow()
				&&
				arg_hit_b.get_start_arrow() < arg_hit_a.get_stop_arrow()
			);
			if ( ! ends_are_ok ) {
				return false;
			}

			return ( ! hits_overlap( arg_hit_a, arg_hit_b ) );
		}



	}
}

#endif
