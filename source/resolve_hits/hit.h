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

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "common/algorithm/transform_build.h"
#include "common/cpp14/cbegin_cend.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "resolve_hits/hit_seg.h"
#include "resolve_hits/res_arrow.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

using namespace cath::common::literals;

namespace cath {
	namespace rslv {

		class hit;
		inline res_arrow get_stop_of_first_segment(const hit &);
		inline res_arrow get_start_of_last_segment(const hit &);

		/// \brief Represent a single hit (ie one domain) with a score, label and one or more segments
		///
		/// Invariants:
		///  * all segments must have the same residue_locating
		///    (ie whether they locate their residues by names and/or indices)
		///
		/// \todo Change internals to use arrows rather than residx_t values and then
		///       make get_*_res_index* methods into non-member functions
		///
		/// Like all cath-resolve-hits code, this assumes simple residue numbering
		/// and is hence unsuitable for use with raw PDB residue numbers.
		///
		/// This is structured to be a reasonably compact 40 bytes:
		///  * 4 x  4-byte numbers
		///  * 1 x 24-byte vector
		/// ...so that it can be processed effieciently in a vector
		class hit final {
		private:
			/// \brief The boundary at the start of the first segment
			res_arrow start_arrow;

			/// \brief The boundary at the end of the last segment
			res_arrow stop_arrow;

			/// \brief The score associated with this hit
			///
			/// This must be greater than 0.0
			resscr_t score;

			/// \brief The index of the label for this hit (in some corresponding list labels)
			hitidx_t label_idx;

			/// \brief The (possibly empty) list of the boundaries associated with any gaps between this hit's segments
			hit_seg_vec fragments;

			void sanity_check() const;

		public:
			hit(const res_arrow &,
			    const res_arrow &,
			    const resscr_t &,
			    const hitidx_t &);

			hit(const hit_seg_vec &,
			    const resscr_t &,
			    const hitidx_t &);

			hit(const res_arrow &,
			    const res_arrow &,
			    const hit_seg_vec &,
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
			const std::string & get_label(const str_vec &) const;

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

			static auto get_hit_first_seg_stop_less() {
				return [] (const hit &x, const hit &y) {
					return ( get_stop_of_first_segment( x ) < get_stop_of_first_segment( y ) );
				};
			}

			static auto get_hit_last_seg_start_less() {
				return [] (const hit &x, const hit &y) {
					return ( get_start_of_last_segment( x ) < get_start_of_last_segment( y ) );
				};
			}
		};

		/// \brief Addition operator to add a hit to a copy of hit_vec
		///
		/// \relates hit_vec
		///
		/// \relatesalso hit
		inline hit_vec operator+(hit_vec    arg_hit_vec, ///< The hit_vec from which to take a copy to which the hit should then be added
		                         const hit &arg_hit      ///< The hit to add
		                         ) {
			arg_hit_vec.push_back( arg_hit );
			return arg_hit_vec;
		}

		/// \brief Represent the different formats in which hits can be output
		enum class hit_output_format {
			CLASS, ///< Default format for outputting C++ classes
			JON    ///< Format in which hits should be output for Jon's use
		};

		std::string get_segments_string(const hit &);
		std::string to_string(const hit &,
		                      const str_vec &,
		                      const hit_output_format & = hit_output_format::CLASS,
		                      const std::string & = std::string{});

		bool operator==(const hit &,
		                const hit &);
		bool any_interaction(const hit &,
		                     const hit &);

		/// \brief Sanity check that the hit is sensible and throw an exception if not
		inline void hit::sanity_check() const {
			if ( stop_arrow < start_arrow ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Start index must not be greater than the stop index"));
			}
			if ( score <= 0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("score must be strictly greater than 0 else the algorithm doesn't work (because there's no way to no how to trade scores off against empty space)"));
			}
		}

		/// \brief Ctor for contiguous hit
		inline hit::hit(const res_arrow &arg_start_arrow, ///< The start boundary of the continuous hit
		                const res_arrow &arg_stop_arrow,  ///< The end boundary of the continuous hit
		                const resscr_t  &arg_score,       ///< The score associated with the hit
		                const hitidx_t  &arg_label_idx    ///< The index of the label associated with the hit (in some other list of hits' labels)
		                ) : start_arrow ( arg_start_arrow ),
		                    stop_arrow  ( arg_stop_arrow  ),
		                    score       ( arg_score       ),
		                    label_idx   ( arg_label_idx       ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous hit from segments
		inline hit::hit(const hit_seg_vec &arg_segments, ///< The segments of the hit
		                const resscr_t    &arg_score,    ///< The score associated with the hit
		                const hitidx_t    &arg_label_idx ///< The index of the label associated with the hit (in some other list of hits' labels)
		                ) : start_arrow ( arg_segments.front().get_start_arrow()      ),
		                    stop_arrow  ( arg_segments.back ().get_stop_arrow ()      ),
		                    score       ( arg_score                                   ),
		                    label_idx   ( arg_label_idx                               ),
		                    fragments   ( make_fragments_of_segments( arg_segments )  ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous hit from start, stop and fragments
		inline hit::hit(const res_arrow   &arg_start_arrow, ///< The boundary at the start of the first segment
		                const res_arrow   &arg_stop_arrow,  ///< The boundary at the end of the last segment
		                const hit_seg_vec &arg_fragments,   ///< The (possibly empty) list of the boundaries associated with any gaps between this hit's segments
		                const resscr_t    &arg_score,       ///< The score associated with the hit
		                const hitidx_t    &arg_label_idx    ///< The index of the label associated with the hit (in some other list of hits' labels)
		                ) : start_arrow ( arg_start_arrow        ),
		                    stop_arrow  ( arg_stop_arrow         ),
		                    score       ( arg_score              ),
		                    label_idx   ( arg_label_idx          ),
		                    fragments   ( arg_fragments          ) {
			sanity_check();
		}

		/// \brief Return whether this hit is discontiguous
		inline bool hit::is_discontig() const {
			return ! fragments.empty();
		}

		/// \brief Return the number of segments in this hit
		inline size_t hit::get_num_segments() const {
			return fragments.size() + 1;
		}

		/// \brief Get the start boundary of the segment with the specified index
		inline const res_arrow & hit::get_start_arrow_of_segment(const size_t &arg_segment_index ///< The index of the segment whose start arrow should be returned
		                                                         ) const {
			return ( arg_segment_index > 0                ) ? fragments[ arg_segment_index - 1 ].get_stop_arrow()
			                                                : start_arrow;
		}

		/// \brief Get the stop boundary of the segment with the specified index
		inline const res_arrow & hit::get_stop_arrow_of_segment(const size_t &arg_segment_index ///< The index of the segment whose stop arrow should be returned
		                                                        ) const {
			return ( arg_segment_index < fragments.size() ) ? fragments[ arg_segment_index     ].get_start_arrow()
			                                                : stop_arrow;
		}

		/// \brief Get the length of the specified hit's segment corresponding to the specified index
		///
		/// \relates hit
		inline size_t get_length_of_hit_seg(const hit    &arg_hit,    ///< The hit to query
		                                    const size_t &arg_seg_idx ///< The index of the segment who length should be returned
		                                    ) {
			return static_cast<size_t>(
				arg_hit.get_stop_arrow_of_segment ( arg_seg_idx ).get_index()
				-
				arg_hit.get_start_arrow_of_segment( arg_seg_idx ).get_index()
			);
		}
		
		/// \brief Get the specified hit's segment corresponding to the specified index
		///
		/// \relates hit
		inline hit_seg get_hit_seg_of_segment(const hit    &arg_hit,    ///< The hit to query
		                                      const size_t &arg_seg_idx ///< The index of the segment to return
		                                      ) {
			return {
				arg_hit.get_start_arrow_of_segment( arg_seg_idx ),
				arg_hit.get_stop_arrow_of_segment ( arg_seg_idx )
			};
		}

		/// \brief Get a vector of the segments in this hit
		///
		/// \relates hit
		inline const hit_seg_vec get_hit_segs(const hit &arg_hit ///< The hit to query
		                                      ) {
			return common::transform_build<hit_seg_vec>(
				boost::irange( 0_z, arg_hit.get_num_segments() ),
				[&] (const size_t &x) {
					return get_hit_seg_of_segment( arg_hit, x );
				}
			);
		}

		/// \brief Append a range of values to a vector of values
		template <typename T, typename U>
		void append(std::vector<T> &arg_vec, ///< The vector to which the values should be appended
		            const U        &arg_rng  ///< The range of values to append
		            ) {
			arg_vec.insert(
				common::cend  ( arg_vec ),
				common::cbegin( arg_rng ),
				common::cend  ( arg_rng )
			);
		}

		/// \brief Get the (possibly-repeated, non-sorted) segments from the specified hits
		///
		/// \relates hit
		inline const hit_seg_vec get_hit_segs(const hit_vec &arg_hit_vec ///< The hits whose segments should be returned
		                                      ) {
			hit_seg_vec results;
			for (const hit &the_hit : arg_hit_vec) {
				append( results, get_hit_segs( the_hit ) );
			}
			return results;
		}

		/// \brief Get a vector of the specified hits' segments, sorted by their starts
		inline hit_seg_vec get_start_sorted_hit_segs(const hit_vec &arg_hit_vec ///< The vector of hits to query
		                                             ) {
			return start_sort_hit_segs_copy( get_hit_segs( arg_hit_vec ) );
		}

		/// \brief Get the (first) start of this hit
		inline const res_arrow & hit::get_start_arrow() const {
			return start_arrow;
		}

		/// \brief Get the (last) stop of this hit
		inline const res_arrow & hit::get_stop_arrow()  const {
			return stop_arrow;
		}

		/// \brief Get the score associated with this hit
		inline const resscr_t & hit::get_score() const {
			return score;
		}

		/// \brief Get the index of the label associated with this hit (where the index refers to some other list of hits' labels)
		inline const hitidx_t & hit::get_label_idx() const {
			return label_idx;
		}

		/// \brief Get the label associated wit this hit from the specified, corresponding list of hits' labels
		inline const std::string & hit::get_label(const str_vec &arg_hit_labels ///< The corresponding list of hits' labels
		                                          ) const {
			return arg_hit_labels[ label_idx ];
		}

		/// \brief Get the start residue index of the segment of specified index in the specified hit
		///
		/// \relates hit
		inline const residx_t & get_start_res_index_of_segment(const hit    &arg_hit,          ///< The hit to query
		                                                       const size_t &arg_segment_index ///< The index of the segment to query
		                                                       ) {
			return arg_hit.get_start_arrow_of_segment( arg_segment_index ).res_after();
		}

		/// \brief Get the stop residue index of the segment of specified index in the specified hit
		///
		/// \relates hit
		inline residx_t get_stop_res_index_of_segment(const hit    &arg_hit,          ///< The hit to query
		                                              const size_t &arg_segment_index ///< The index of the segment to query
		                                              ) {
			return arg_hit.get_stop_arrow_of_segment( arg_segment_index ).res_before();
		}

		/// \brief Get the start residue index of the specified hit
		///
		/// \relates hit
		inline const residx_t & get_start_res_index(const hit &arg_hit ///< The hit to query
		                                            ) {
			return arg_hit.get_start_arrow().res_after();
		}

		/// \brief Get the stop residue index of the specified hit
		///
		/// \relates hit
		inline residx_t get_stop_res_index(const hit &arg_hit ///< The hit to query
		                                   ) {
			return arg_hit.get_stop_arrow().res_before();
		}

		/// \brief Get the stop of the first segment in the specified hit
		///
		/// \pre `arg_hit.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates hit
		inline res_arrow get_stop_of_first_segment(const hit &arg_hit ///< The hit to query
		                                           ) {
			if ( ! arg_hit.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_stop_of_first_segment of contiguous hit"));
			}
			return arg_hit.get_stop_arrow_of_segment( 0 );
		}

		/// \brief Get the start of the last segment in the specified hit
		///
		/// \pre `arg_hit.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates hit
		inline res_arrow get_start_of_last_segment(const hit &arg_hit ///< The hit to query
		                                           ) {
			if ( ! arg_hit.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_start_of_last_segment of contiguous hit"));
			}
			return arg_hit.get_start_arrow_of_segment( arg_hit.get_num_segments() - 1 );
		}

		/// \brief Get the total length of the specified hit (ie the sum of its segments' lengths)
		///
		/// \relates hit
		inline size_t get_total_length(const hit &arg_hit ///< The hit to query
		                               ) {
			return boost::accumulate(
				boost::irange( 0_z, arg_hit.get_num_segments() )
					| boost::adaptors::transformed( [&] (const size_t &x) {
						return get_length_of_hit_seg( arg_hit, x );
					} ),
				0_z
			);
		}

		/// \brief Make a continuous hit from the residue indices
		///
		/// \relates hit
		inline hit make_hit_from_res_indices(const residx_t &arg_start_res_idx, ///< The start residue index
		                                     const residx_t &arg_stop_res_idx,  ///< The stop residue index
		                                     const resscr_t &arg_score,         ///< The score associated with the hit
		                                     const hitidx_t &arg_label_idx      ///< The index of the label associated with the hit (in some other list of hits' labels)
		                                     ) {
			return {
				arrow_before_res( arg_start_res_idx ),
				arrow_after_res ( arg_stop_res_idx  ),
				arg_score,
				arg_label_idx
			};
		}

		/// \brief Make a hit
		///
		/// \relates hit
		inline hit make_hit_from_res_indices(const residx_residx_pair_vec &arg_residue_index_segments, ///< The residue index start/stop pairs of the hit's segments
		                                     const resscr_t               &arg_score,                  ///< The score associated with the hit
		                                     const hitidx_t               &arg_label_idx               ///< The index of the label associated with the hit (in some other list of hits' labels)
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
		/// \relates hit
		inline bool any_interaction(const hit &arg_hit_a, ///< The first  hit to query
		                            const hit &arg_hit_b  ///< The second hit to query
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
		/// hit interleaves or straddles the other
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \todo This can be made more efficient by iterating over the two segment
		///       lists simultaneously
		///
		/// \relates hit
		inline bool hits_overlap(const hit &arg_hit_a, ///< The first  hit to query
		                         const hit &arg_hit_b  ///< The second hit to query
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

		/// \brief Return whether the specified hit overlaps with any of the hits in the specified list of hits
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \relates hit
		inline bool hit_overlaps_with_any_of_hits(const hit     &arg_hit_lhs, ///< The hit to query
		                                          const hit_vec &arg_hit_vec  ///< The list of hits to query
		                                          ) {
			for (const hit &hit_rhs : arg_hit_vec) {
				if ( hits_overlap( arg_hit_lhs, hit_rhs ) ) {
					return true;
				}
			}
			return false;
		}

		/// \brief Return whether the second hit right-interspersed the first
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
		/// \relates hit
		inline bool second_right_intersperses_first(const hit &arg_hit_a, ///< The first  hit to query
		                                            const hit &arg_hit_b  ///< The second hit to query
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
