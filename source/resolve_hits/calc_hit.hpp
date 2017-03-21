/// \file
/// \brief The calc_hit class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_CALC_HIT_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_CALC_HIT_H

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "common/algorithm/append.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "resolve_hits/hit_seg.hpp"
#include "resolve_hits/res_arrow.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

using namespace cath::common::literals;

namespace cath {
	namespace rslv {

		class calc_hit;
		inline res_arrow get_stop_of_first_segment(const calc_hit &);
		inline res_arrow get_start_of_last_segment(const calc_hit &);

		/// \brief Represent a single calc_hit (ie one domain) with a score, label and one or more segments
		///
		/// Like all cath-resolve-hits code, this assumes simple residue numbering
		/// and is hence unsuitable for use with raw PDB residue numbers.
		///
		/// This is structured to be a reasonably compact 40 bytes:
		///  * 4 x  4-byte numbers
		///  * 1 x 24-byte vector
		/// ...so that it can be processed efficiently in a vector
		class calc_hit final {
		private:
			/// \brief The boundary at the start of the first segment
			res_arrow start_arrow;

			/// \brief The boundary at the end of the last segment
			res_arrow stop_arrow;

			/// \brief The score associated with this calc_hit
			///
			/// This must be greater than 0.0
			resscr_t score;

			/// \brief The index of the label for this calc_hit (in some corresponding list labels)
			hitidx_t label_idx;

			/// \brief The (possibly empty) list of the boundaries associated with any gaps between this calc_hit's segments
			hit_seg_vec fragments;

			void sanity_check() const;

		public:
			calc_hit(res_arrow,
			         res_arrow,
			         const resscr_t &,
			         const hitidx_t &);

			calc_hit(const hit_seg_vec &,
			         const resscr_t &,
			         const hitidx_t &);

			calc_hit(res_arrow,
			         res_arrow,
			         hit_seg_vec,
			         const resscr_t &,
			         const hitidx_t &);

			bool is_discontig() const;
			size_t get_num_segments() const;
			const res_arrow & get_start_arrow_of_segment(const size_t &) const;
			const res_arrow & get_stop_arrow_of_segment(const size_t &) const;

			const res_arrow   & get_start_arrow() const;
			const res_arrow   & get_stop_arrow()  const;
			const resscr_t    & get_score()       const;
			const hitidx_t    & get_label_idx() const;

			static auto get_hit_start_less() {
				return [] (const calc_hit &x, const calc_hit &y) {
					return ( x.get_start_arrow() < y.get_start_arrow() );
				};
			}

			static auto get_hit_stop_less() {
				return [] (const calc_hit &x, const calc_hit &y) {
					return ( x.get_stop_arrow() < y.get_stop_arrow() );
				};
			}

			static auto get_hit_first_seg_stop_less() {
				return [] (const calc_hit &x, const calc_hit &y) {
					return ( get_stop_of_first_segment( x ) < get_stop_of_first_segment( y ) );
				};
			}

			static auto get_hit_last_seg_start_less() {
				return [] (const calc_hit &x, const calc_hit &y) {
					return ( get_start_of_last_segment( x ) < get_start_of_last_segment( y ) );
				};
			}
		};

		/// \brief Addition operator to add a calc_hit to a copy of calc_hit_vec
		///
		/// \relates calc_hit_vec
		///
		/// \relatesalso calc_hit
		inline calc_hit_vec operator+(calc_hit_vec    arg_hit_vec, ///< The calc_hit_vec from which to take a copy to which the calc_hit should then be added
		                              const calc_hit &arg_hit      ///< The calc_hit to add
		                              ) {
			arg_hit_vec.push_back( arg_hit );
			return arg_hit_vec;
		}

		std::string get_segments_string(const calc_hit &);
		std::string to_string(const calc_hit &);
		std::ostream & operator<<(std::ostream &,
		                          const calc_hit &);
		bool operator==(const calc_hit &,
		                const calc_hit &);
		bool any_interaction(const calc_hit &,
		                     const calc_hit &);

		/// \brief Sanity check that the calc_hit is sensible and throw an exception if not
		inline void calc_hit::sanity_check() const {
			if ( stop_arrow < start_arrow ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Start index must not be greater than the stop index"));
			}
			if ( score <= 0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
					"Hit cannot be processed because its score of "
					+ ::std::to_string( get_score() )
					+ " isn't greater than 0, which is required for the algorithm to work (because otherwise there's no way to know how to trade scores off against empty space)"
				));
			}
		}

		/// \brief Ctor for contiguous calc_hit
		inline calc_hit::calc_hit(res_arrow        arg_start_arrow, ///< The start boundary of the continuous calc_hit
		                          res_arrow        arg_stop_arrow,  ///< The end boundary of the continuous calc_hit
		                          const resscr_t  &arg_score,       ///< The score associated with the calc_hit
		                          const hitidx_t  &arg_label_idx    ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
		                          ) : start_arrow ( std::move( arg_start_arrow ) ),
		                              stop_arrow  ( std::move( arg_stop_arrow  ) ),
		                              score       ( arg_score                    ),
		                              label_idx   ( arg_label_idx                ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous calc_hit from segments
		inline calc_hit::calc_hit(const hit_seg_vec &arg_segments, ///< The segments of the calc_hit
		                          const resscr_t    &arg_score,    ///< The score associated with the calc_hit
		                          const hitidx_t    &arg_label_idx ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
		                          ) : start_arrow ( arg_segments.front().get_start_arrow()      ),
		                              stop_arrow  ( arg_segments.back ().get_stop_arrow ()      ),
		                              score       ( arg_score                                   ),
		                              label_idx   ( arg_label_idx                               ),
		                              fragments   ( make_fragments_of_segments( arg_segments )  ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous calc_hit from start, stop and fragments
		inline calc_hit::calc_hit(res_arrow arg_start_arrow,       ///< The boundary at the start of the first segment
		                          res_arrow arg_stop_arrow,        ///< The boundary at the end of the last segment
		                          hit_seg_vec arg_fragments,       ///< The (possibly empty) list of the boundaries associated with any gaps between this calc_hit's segments
		                          const resscr_t    &arg_score,    ///< The score associated with the calc_hit
		                          const hitidx_t    &arg_label_idx ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
		                          ) : start_arrow (std::move( arg_start_arrow        )),
		                              stop_arrow  (std::move( arg_stop_arrow         )),
		                              score       ( arg_score              ),
		                              label_idx   ( arg_label_idx          ),
		                              fragments   (std::move( arg_fragments          )) {
			sanity_check();
		}

		/// \brief Return whether this calc_hit is discontiguous
		inline bool calc_hit::is_discontig() const {
			return ! fragments.empty();
		}

		/// \brief Return the number of segments in this calc_hit
		inline size_t calc_hit::get_num_segments() const {
			return fragments.size() + 1;
		}

		/// \brief Get the start boundary of the segment with the specified index
		inline const res_arrow & calc_hit::get_start_arrow_of_segment(const size_t &arg_segment_index ///< The index of the segment whose start arrow should be returned
		                                                              ) const {
			return ( arg_segment_index > 0                ) ? fragments[ arg_segment_index - 1 ].get_stop_arrow()
			                                                : start_arrow;
		}

		/// \brief Get the stop boundary of the segment with the specified index
		inline const res_arrow & calc_hit::get_stop_arrow_of_segment(const size_t &arg_segment_index ///< The index of the segment whose stop arrow should be returned
		                                                             ) const {
			return ( arg_segment_index < fragments.size() ) ? fragments[ arg_segment_index     ].get_start_arrow()
			                                                : stop_arrow;
		}

		/// \brief Get the length of the specified calc_hit's segment corresponding to the specified index
		///
		/// \relates calc_hit
		inline size_t get_length_of_hit_seg(const calc_hit &arg_hit,    ///< The calc_hit to query
		                                    const size_t   &arg_seg_idx ///< The index of the segment who length should be returned
		                                    ) {
			return static_cast<size_t>(
				arg_hit.get_stop_arrow_of_segment ( arg_seg_idx )
				-
				arg_hit.get_start_arrow_of_segment( arg_seg_idx )
			);
		}
		
		/// \brief Get the specified calc_hit's segment corresponding to the specified index
		///
		/// \relates calc_hit
		inline hit_seg get_hit_seg_of_seg_idx(const calc_hit &arg_hit,    ///< The calc_hit to query
		                                      const size_t   &arg_seg_idx ///< The index of the segment to return
		                                      ) {
			return {
				arg_hit.get_start_arrow_of_segment( arg_seg_idx ),
				arg_hit.get_stop_arrow_of_segment ( arg_seg_idx )
			};
		}

		/// \brief Get a vector of the segments in this calc_hit
		///
		/// \relates calc_hit
		inline const hit_seg_vec get_hit_segs(const calc_hit &arg_hit ///< The calc_hit to query
		                                      ) {
			return common::transform_build<hit_seg_vec>(
				boost::irange( 0_z, arg_hit.get_num_segments() ),
				[&] (const size_t &x) {
					return get_hit_seg_of_seg_idx( arg_hit, x );
				}
			);
		}

		/// \brief Get the (possibly-repeated, non-sorted) segments from the specified hits
		///
		/// \relates calc_hit
		inline const hit_seg_vec get_hit_segs(const calc_hit_vec &arg_hit_vec ///< The hits whose segments should be returned
		                                      ) {
			hit_seg_vec results;
			for (const calc_hit &the_hit : arg_hit_vec) {
				common::append( results, get_hit_segs( the_hit ) );
			}
			return results;
		}

		/// \brief Get a vector of the specified hits' segments, sorted by their starts
		inline hit_seg_vec get_start_sorted_hit_segs(const calc_hit_vec &arg_hit_vec ///< The vector of hits to query
		                                             ) {
			return start_sort_hit_segs_copy( get_hit_segs( arg_hit_vec ) );
		}

		/// \brief Get the (first) start of this calc_hit
		inline const res_arrow & calc_hit::get_start_arrow() const {
			return start_arrow;
		}

		/// \brief Get the (last) stop of this calc_hit
		inline const res_arrow & calc_hit::get_stop_arrow()  const {
			return stop_arrow;
		}

		/// \brief Get the score associated with this calc_hit
		inline const resscr_t & calc_hit::get_score() const {
			return score;
		}

		/// \brief Get the index of the label associated with this calc_hit (where the index refers to some other list of hits' labels)
		inline const hitidx_t & calc_hit::get_label_idx() const {
			return label_idx;
		}

		/// \brief Get the start residue index of the segment of specified index in the specified calc_hit
		///
		/// \relates calc_hit
		inline const residx_t & get_start_res_index_of_segment(const calc_hit &arg_hit,          ///< The calc_hit to query
		                                                       const size_t   &arg_segment_index ///< The index of the segment to query
		                                                       ) {
			return arg_hit.get_start_arrow_of_segment( arg_segment_index ).res_after();
		}

		/// \brief Get the stop residue index of the segment of specified index in the specified calc_hit
		///
		/// \relates calc_hit
		inline residx_t get_stop_res_index_of_segment(const calc_hit &arg_hit,          ///< The calc_hit to query
		                                              const size_t   &arg_segment_index ///< The index of the segment to query
		                                              ) {
			return arg_hit.get_stop_arrow_of_segment( arg_segment_index ).res_before();
		}

		/// \brief Get the start residue index of the specified calc_hit
		///
		/// \relates calc_hit
		inline const residx_t & get_start_res_index(const calc_hit &arg_hit ///< The calc_hit to query
		                                            ) {
			return arg_hit.get_start_arrow().res_after();
		}

		/// \brief Get the stop residue index of the specified calc_hit
		///
		/// \relates calc_hit
		inline residx_t get_stop_res_index(const calc_hit &arg_hit ///< The calc_hit to query
		                                   ) {
			return arg_hit.get_stop_arrow().res_before();
		}

		/// \brief Get the stop of the first segment in the specified calc_hit
		///
		/// \pre `arg_hit.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates calc_hit
		inline res_arrow get_stop_of_first_segment(const calc_hit &arg_hit ///< The calc_hit to query
		                                           ) {
			if ( ! arg_hit.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_stop_of_first_segment of contiguous calc_hit"));
			}
			return arg_hit.get_stop_arrow_of_segment( 0 );
		}

		/// \brief Get the start of the last segment in the specified calc_hit
		///
		/// \pre `arg_hit.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates calc_hit
		inline res_arrow get_start_of_last_segment(const calc_hit &arg_hit ///< The calc_hit to query
		                                           ) {
			if ( ! arg_hit.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_start_of_last_segment of contiguous calc_hit"));
			}
			return arg_hit.get_start_arrow_of_segment( arg_hit.get_num_segments() - 1 );
		}

		/// \brief Get the total length of the specified calc_hit (ie the sum of its segments' lengths)
		///
		/// \relates calc_hit
		inline size_t get_total_length(const calc_hit &arg_hit ///< The calc_hit to query
		                               ) {
			return boost::accumulate(
				boost::irange( 0_z, arg_hit.get_num_segments() )
					| boost::adaptors::transformed( [&] (const size_t &x) {
						return get_length_of_hit_seg( arg_hit, x );
					} ),
				0_z
			);
		}

		/// \brief Make a continuous calc_hit from the residue indices
		///
		/// \relates calc_hit
		inline calc_hit make_hit_from_res_indices(const residx_t &arg_start_res_idx, ///< The start residue index
		                                          const residx_t &arg_stop_res_idx,  ///< The stop residue index
		                                          const resscr_t &arg_score,         ///< The score associated with the calc_hit
		                                          const hitidx_t &arg_label_idx      ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
		                                          ) {
			return {
				arrow_before_res( arg_start_res_idx ),
				arrow_after_res ( arg_stop_res_idx  ),
				arg_score,
				arg_label_idx
			};
		}

		/// \brief Make a calc_hit
		///
		/// \relates calc_hit
		inline calc_hit make_hit_from_res_indices(const residx_residx_pair_vec &arg_residue_index_segments, ///< The residue index start/stop pairs of the calc_hit's segments
		                                          const resscr_t               &arg_score,                  ///< The score associated with the calc_hit
		                                          const hitidx_t               &arg_label_idx               ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
		                                          ) {
			return {
				common::transform_build<hit_seg_vec>(
					arg_residue_index_segments,
					hit_seg_of_res_idx_pair
				),
				arg_score,
				arg_label_idx
			};
		}

		/// \brief Return whether the either of the two specified hits overlaps, interleaves or straddles the other
		///
		/// \relates calc_hit
		inline bool any_interaction(const calc_hit &arg_hit_a, ///< The first  calc_hit to query
		                            const calc_hit &arg_hit_b  ///< The second calc_hit to query
		                            ) {
			return (
				arg_hit_a.get_start_arrow() < arg_hit_b.get_stop_arrow()
				&&
				arg_hit_b.get_start_arrow() < arg_hit_a.get_stop_arrow()
			);
		}
		
		/// \brief Return whether the two specified hits overlap with each other
		///
		/// This requires there to be a genuine overlap of segments, not just that one
		/// calc_hit interleaves or straddles the other
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \todo This can be made more efficient by iterating over the two segment
		///       lists simultaneously
		///
		/// \relates calc_hit
		inline bool hits_overlap(const calc_hit &arg_hit_a, ///< The first  calc_hit to query
		                         const calc_hit &arg_hit_b  ///< The second calc_hit to query
		                         ) {
			if ( ! any_interaction( arg_hit_a, arg_hit_b ) ) {
				return false;
			}
			for (const auto &seg_ctr_a : boost::irange( 0_z, arg_hit_a.get_num_segments() ) ) {
				for (const auto &seg_ctr_b : boost::irange( 0_z, arg_hit_b.get_num_segments() ) ) {
					const bool seg_overlap = hit_segs_overlap(
						get_hit_seg_of_seg_idx( arg_hit_a, seg_ctr_a ),
						get_hit_seg_of_seg_idx( arg_hit_b, seg_ctr_b )
					);
					if ( seg_overlap ) {
						return true;
					}
				}
			}
			return false;
		}

		/// \brief Return whether the specified calc_hit overlaps with any of the hits in the specified list of hits
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \relates calc_hit
		inline bool hit_overlaps_with_any_of_hits(const calc_hit &arg_hit_lhs, ///< The calc_hit to query
		                                          const calc_hit_vec  &arg_hit_vec  ///< The list of hits to query
		                                          ) {
			for (const calc_hit &hit_rhs : arg_hit_vec) {
				if ( hits_overlap( arg_hit_lhs, hit_rhs ) ) {
					return true;
				}
			}
			return false;
		}

		/// \brief Return whether the second calc_hit right-intersperses the first
		///
		/// This means that both are discontiguous hits that don't actually
		/// overlap each other but the second's start is within the first
		/// and the first's stop in within the second. Like this:
		///
		/// ~~~~~
		///   ***     ***
		///       ***     ***
		/// ~~~~~
		///
		/// \relates calc_hit
		inline bool second_right_intersperses_first(const calc_hit &arg_hit_a, ///< The first  calc_hit to query
		                                            const calc_hit &arg_hit_b  ///< The second calc_hit to query
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

		/// \brief Return whether the second calc_hit right-intersperses or inside-intersperses the first
		///
		/// This means that both are discontiguous hits that don't actually
		/// overlap each other but the second's start is within the first
		///
		/// ~~~~~
		///   ***     ***     ***
		///       ***     ***
		/// ~~~~~
		///
		/// \relates calc_hit
		inline bool second_right_or_inside_intersperses_first(const calc_hit &arg_hit_a, ///< The first  calc_hit to query
		                                                      const calc_hit &arg_hit_b  ///< The second calc_hit to query
		                                                      ) {
			if ( ! arg_hit_a.is_discontig() || ! arg_hit_b.is_discontig() ) {
				return false;
			}
			const bool ends_are_ok = (
				arg_hit_a.get_start_arrow() < arg_hit_b.get_start_arrow()
				&&
				arg_hit_b.get_start_arrow() < arg_hit_a.get_stop_arrow()
			);
			if ( ! ends_are_ok ) {
				return false;
			}

			return ( ! hits_overlap( arg_hit_a, arg_hit_b ) );
		}

	} // namespace rslv
} // namespace cath

#endif
