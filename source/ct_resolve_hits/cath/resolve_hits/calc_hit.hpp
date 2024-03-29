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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_CALC_HIT_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_CALC_HIT_HPP

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/numeric.hpp>

#include "cath/common/algorithm/append.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"
#include "cath/seq/seq_arrow.hpp"
#include "cath/seq/seq_seg.hpp"
#include "cath/seq/seq_seg_run.hpp"

// clang-format off
namespace cath::rslv { class calc_hit; }
// clang-format on

namespace cath::rslv {

	inline const seq::seq_arrow & get_start_arrow(const calc_hit &);
	inline const seq::seq_arrow & get_stop_arrow(const calc_hit &);
	inline seq::seq_arrow get_start_of_last_segment(const calc_hit &);
	inline seq::seq_arrow get_stop_of_first_segment(const calc_hit &);

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
		/// \brief The score associated with this calc_hit
		///
		/// This must be greater than 0.0
		resscr_t score;

		/// \brief The index of the label for this calc_hit (in some corresponding list labels)
		hitidx_t label_idx;

		/// \brief The list of the segments
		seq::seq_seg_run segments;

		void sanity_check() const;

	public:
		calc_hit(seq::seq_arrow,
		         seq::seq_arrow,
		         const resscr_t &,
		         const hitidx_t &);

		calc_hit(const seq::seq_seg_vec &,
		         const resscr_t &,
		         const hitidx_t &);

		calc_hit(seq::seq_arrow,
		         seq::seq_arrow,
		         seq::seq_seg_vec,
		         const resscr_t &,
		         const hitidx_t &);

		[[nodiscard]] const seq::seq_seg_run &get_segments() const;
		[[nodiscard]] const resscr_t &        get_score() const;
		[[nodiscard]] const hitidx_t &        get_label_idx() const;

		static auto get_hit_start_less() {
			return [] (const calc_hit &x, const calc_hit &y) {
				return ( get_start_arrow( x ) < get_start_arrow( y ) );
			};
		}

		static auto get_hit_stop_less() {
			return [] (const calc_hit &x, const calc_hit &y) {
				return ( get_stop_arrow( x ) < get_stop_arrow( y ) );
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
	inline calc_hit_vec operator+(calc_hit_vec    prm_hit_vec, ///< The calc_hit_vec from which to take a copy to which the calc_hit should then be added
	                              const calc_hit &prm_hit      ///< The calc_hit to add
	                              ) {
		prm_hit_vec.push_back( prm_hit );
		return prm_hit_vec;
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
		if ( score <= 0 ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception(
				"Hit cannot be processed because its score of "
				+ ::std::to_string( get_score() )
				+ " isn't greater than 0, which is required for the algorithm to work (because otherwise there's no way to know how to trade scores off against empty space)"
			));
		}
	}

	/// \brief Ctor for contiguous calc_hit
	inline calc_hit::calc_hit(seq::seq_arrow   prm_start_arrow, ///< The start boundary of the continuous calc_hit
	                          seq::seq_arrow   prm_stop_arrow,  ///< The end boundary of the continuous calc_hit
	                          const resscr_t  &prm_score,       ///< The score associated with the calc_hit
	                          const hitidx_t  &prm_label_idx    ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
	                          ) : score     { prm_score     },
	                              label_idx { prm_label_idx },
	                              segments  {
	                              	std::move( prm_start_arrow ),
	                              	std::move( prm_stop_arrow )
	                              } {
		sanity_check();
	}

	/// \brief Ctor for a possibly discontinuous calc_hit from segments
	inline calc_hit::calc_hit( const seq::seq_seg_vec &prm_segments, ///< The segments of the calc_hit
	                           const resscr_t &        prm_score,    ///< The score associated with the calc_hit
	                           const hitidx_t &        prm_label_idx ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
	                           ) :
	        score{ prm_score }, label_idx{ prm_label_idx }, segments{ prm_segments } {
		sanity_check();
	}

	/// \brief Ctor for a possibly discontinuous calc_hit from start, stop and fragments
	inline calc_hit::calc_hit(seq::seq_arrow     prm_start_arrow, ///< The boundary at the start of the first segment
	                          seq::seq_arrow     prm_stop_arrow,  ///< The boundary at the end of the last segment
	                          seq::seq_seg_vec   prm_fragments,   ///< The (possibly empty) list of the boundaries associated with any gaps between this calc_hit's segments
	                          const resscr_t    &prm_score,       ///< The score associated with the calc_hit
	                          const hitidx_t    &prm_label_idx    ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
	                          ) : score       { prm_score     },
	                              label_idx   { prm_label_idx },
	                              segments    {
	                              	std::move( prm_start_arrow ),
	                              	std::move( prm_stop_arrow ),
	                              	std::move( prm_fragments )
	                              } {
		sanity_check();
	}


	/// \brief Get the segments associated with this calc_hit
	inline const seq::seq_seg_run & calc_hit::get_segments() const {
		return segments;
	}

	/// \brief Get the score associated with this calc_hit
	inline const resscr_t & calc_hit::get_score() const {
		return score;
	}

	/// \brief Get the index of the label associated with this calc_hit (where the index refers to some other list of hits' labels)
	inline const hitidx_t & calc_hit::get_label_idx() const {
		return label_idx;
	}

	/// \brief Return whether this calc_hit is discontiguous
	inline bool is_discontig(const calc_hit &prm_calc_hit ///< The calc_hit to query
	                         ) {
		return prm_calc_hit.get_segments().is_discontig();
	}

	/// \brief Return the number of segments in this calc_hit
	inline size_t get_num_segments(const calc_hit &prm_calc_hit ///< The calc_hit to query
	                               ) {
		return prm_calc_hit.get_segments().get_num_segments();
	}

	/// \brief Get the start boundary of the segment with the specified index
	inline const seq::seq_arrow & get_start_arrow_of_segment(const calc_hit &prm_calc_hit,     ///< The calc_hit to query
	                                                         const size_t   &prm_segment_index ///< The index of the segment whose start arrow should be returned
	                                                         ) {
		return prm_calc_hit.get_segments().get_start_arrow_of_segment( prm_segment_index );
	}

	/// \brief Get the stop boundary of the segment with the specified index
	inline const seq::seq_arrow & get_stop_arrow_of_segment(const calc_hit &prm_calc_hit,     ///< The calc_hit to query
	                                                        const size_t   &prm_segment_index ///< The index of the segment whose stop arrow should be returned
	                                                        ) {
		return prm_calc_hit.get_segments().get_stop_arrow_of_segment( prm_segment_index );
	}

	/// \brief Get the (first) start of this calc_hit
	inline const seq::seq_arrow & get_start_arrow(const calc_hit &prm_calc_hit ///< The calc_hit to query
	                                              ) {
		return prm_calc_hit.get_segments().get_start_arrow();
	}

	/// \brief Get the (last) stop of this calc_hit
	inline const seq::seq_arrow & get_stop_arrow(const calc_hit &prm_calc_hit ///< The calc_hit to query
	                                             ) {
		return prm_calc_hit.get_segments().get_stop_arrow();
	}


	/// \brief Get the length of the specified calc_hit's segment corresponding to the specified index
	///
	/// \relates calc_hit
	inline size_t get_length_of_seq_seg(const calc_hit &prm_hit,    ///< The calc_hit to query
	                                    const size_t   &prm_seg_idx ///< The index of the segment who length should be returned
	                                    ) {
		return get_length_of_seq_seg( prm_hit.get_segments(), prm_seg_idx );
	}
	
	/// \brief Get the specified calc_hit's segment corresponding to the specified index
	///
	/// \relates calc_hit
	inline seq::seq_seg get_seq_seg_of_seg_idx(const calc_hit &prm_hit,    ///< The calc_hit to query
	                                           const size_t   &prm_seg_idx ///< The index of the segment to return
	                                           ) {
		return get_seq_seg_of_seg_idx( prm_hit.get_segments(), prm_seg_idx );
	}

	/// \brief Get a vector of the segments in this calc_hit
	///
	/// \relates calc_hit
	inline seq::seq_seg_vec get_seq_segs( const calc_hit &prm_hit ///< The calc_hit to query
	                                      ) {
		return get_seq_segs( prm_hit.get_segments() );
	}

	/// \brief Get the (possibly-repeated, non-sorted) segments from the specified hits
	///
	/// \relates calc_hit
	inline seq::seq_seg_vec get_seq_segs( const calc_hit_vec &prm_hit_vec ///< The hits whose segments should be returned
	                                      ) {
		seq::seq_seg_vec results;
		for (const calc_hit &the_hit : prm_hit_vec) {
			common::append( results, get_seq_segs( the_hit ) );
		}
		return results;
	}

	/// \brief Get a vector of the specified hits' segments, sorted by their starts
	inline seq::seq_seg_vec get_start_sorted_seq_segs(const calc_hit_vec &prm_hit_vec ///< The vector of hits to query
	                                                  ) {
		return start_sort_seq_segs_copy( get_seq_segs( prm_hit_vec ) );
	}


	/// \brief Get the start residue index of the segment of specified index in the specified calc_hit
	///
	/// \relates calc_hit
	inline const seq::residx_t & get_start_res_index_of_segment(const calc_hit &prm_hit,          ///< The calc_hit to query
	                                                            const size_t   &prm_segment_index ///< The index of the segment to query
	                                                            ) {
		return get_start_res_index_of_segment( prm_hit.get_segments(), prm_segment_index );
	}

	/// \brief Get the stop residue index of the segment of specified index in the specified calc_hit
	///
	/// \relates calc_hit
	inline seq::residx_t get_stop_res_index_of_segment(const calc_hit &prm_hit,          ///< The calc_hit to query
	                                                   const size_t   &prm_segment_index ///< The index of the segment to query
	                                                   ) {
		return get_stop_res_index_of_segment( prm_hit.get_segments(), prm_segment_index );
	}

	/// \brief Get the start residue index of the specified calc_hit
	///
	/// \relates calc_hit
	inline const seq::residx_t & get_start_res_index(const calc_hit &prm_hit ///< The calc_hit to query
	                                                 ) {
		return get_start_res_index( prm_hit.get_segments() );
	}

	/// \brief Get the stop residue index of the specified calc_hit
	///
	/// \relates calc_hit
	inline seq::residx_t get_stop_res_index(const calc_hit &prm_hit ///< The calc_hit to query
	                                        ) {
		return get_stop_res_index( prm_hit.get_segments() );
	}

	/// \brief Get the stop of the first segment in the specified calc_hit
	///
	/// \pre `prm_hit.is_discontig()` else an invalid_argument_exception will be thrown
	///
	/// \relates calc_hit
	inline seq::seq_arrow get_stop_of_first_segment(const calc_hit &prm_hit ///< The calc_hit to query
	                                                ) {
		return get_stop_of_first_segment( prm_hit.get_segments() );
	}

	/// \brief Get the start of the last segment in the specified calc_hit
	///
	/// \pre `prm_hit.is_discontig()` else an invalid_argument_exception will be thrown
	///
	/// \relates calc_hit
	inline seq::seq_arrow get_start_of_last_segment(const calc_hit &prm_hit ///< The calc_hit to query
	                                                ) {
		return get_start_of_last_segment( prm_hit.get_segments() );
	}

	/// \brief Get the total length of the specified calc_hit (ie the sum of its segments' lengths)
	///
	/// \relates calc_hit
	inline size_t get_total_length(const calc_hit &prm_hit ///< The calc_hit to query
	                               ) {
		return get_total_length( prm_hit.get_segments() );
	}

	/// \brief Make a continuous calc_hit from the residue indices
	///
	/// \relates calc_hit
	inline calc_hit make_hit_from_res_indices(const seq::residx_t &prm_start_res_idx, ///< The start residue index
	                                          const seq::residx_t &prm_stop_res_idx,  ///< The stop residue index
	                                          const resscr_t      &prm_score,         ///< The score associated with the calc_hit
	                                          const hitidx_t      &prm_label_idx      ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
	                                          ) {
		return {
			seq::arrow_before_res( prm_start_res_idx ),
			seq::arrow_after_res ( prm_stop_res_idx  ),
			prm_score,
			prm_label_idx
		};
	}

	/// \brief Make a calc_hit
	///
	/// \relates calc_hit
	inline calc_hit make_hit_from_res_indices(const seq::residx_residx_pair_vec &prm_residue_index_segments, ///< The residue index start/stop pairs of the calc_hit's segments
	                                          const resscr_t                    &prm_score,                  ///< The score associated with the calc_hit
	                                          const hitidx_t                    &prm_label_idx               ///< The index of the label associated with the calc_hit (in some other list of hits' labels)
	                                          ) {
		return {
			common::transform_build<seq::seq_seg_vec>(
				prm_residue_index_segments,
				seq::seq_seg_of_res_idx_pair
			),
			prm_score,
			prm_label_idx
		};
	}

	/// \brief Return whether the either of the two specified hits overlaps, interleaves or straddles the other
	///
	/// \relates calc_hit
	inline bool any_interaction(const calc_hit &prm_hit_a, ///< The first  calc_hit to query
	                            const calc_hit &prm_hit_b  ///< The second calc_hit to query
	                            ) {
		return any_interaction( prm_hit_a.get_segments(), prm_hit_b.get_segments() );
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
	inline bool are_overlapping(const calc_hit &prm_hit_a, ///< The first  calc_hit to query
	                            const calc_hit &prm_hit_b  ///< The second calc_hit to query
	                            ) {
		return are_overlapping( prm_hit_a.get_segments(), prm_hit_b.get_segments() );
	}

	/// \brief Return whether the specified calc_hit overlaps with any of the hits in the specified list of hits
	///
	/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
	///
	/// \relates calc_hit
	inline bool hit_overlaps_with_any_of_hits(const calc_hit     &prm_hit_lhs, ///< The calc_hit to query
	                                          const calc_hit_vec &prm_hit_vec  ///< The list of hits to query
	                                          ) {
		return ::boost::algorithm::any_of(
		  prm_hit_vec, [ & ]( const calc_hit &x ) { return are_overlapping( prm_hit_lhs, x ); } );
	}

	/// \brief Whether the segments in the first specified calc_hit never extend outside
	///        those in the second specified calc_hit
	///
	/// \relates calc_hit
	inline bool first_is_not_outside_second(const calc_hit &prm_hit_lhs, ///< The first  calc_hit to query
	                                        const calc_hit &prm_hit_rhs  ///< The second calc_hit to query
	                                        ) {
		return first_is_not_outside_second( prm_hit_lhs.get_segments(), prm_hit_rhs.get_segments() );
	}

	/// \brief Whether either of the specified calc_hit covers the other
	///
	/// \relates calc_hit
	inline bool one_covers_other(const calc_hit &prm_hit_lhs, ///< The first  calc_hit to query
	                             const calc_hit &prm_hit_rhs  ///< The second calc_hit to query
	                             ) {
		return one_covers_other( prm_hit_lhs.get_segments(), prm_hit_rhs.get_segments() );
	}

	/// \brief Whether the segments in the first specified calc_hit are shorter strictly
	///        within those  in the second specified calc_hit
	///
	/// \relates calc_hit
	inline bool first_is_shorter_and_within_second(const calc_hit &prm_hit_lhs, ///< The first  calc_hit to query
	                                               const calc_hit &prm_hit_rhs  ///< The second calc_hit to query
	                                               ) {
		return first_is_shorter_and_within_second( prm_hit_lhs.get_segments(), prm_hit_rhs.get_segments() );
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
	inline bool second_right_intersperses_first(const calc_hit &prm_hit_lhs, ///< The first  calc_hit to query
	                                            const calc_hit &prm_hit_rhs  ///< The second calc_hit to query
	                                            ) {
		return second_right_intersperses_first( prm_hit_lhs.get_segments(), prm_hit_rhs.get_segments() );
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
	inline bool second_right_or_inside_intersperses_first(const calc_hit &prm_hit_lhs, ///< The first  calc_hit to query
	                                                      const calc_hit &prm_hit_rhs  ///< The second calc_hit to query
	                                                      ) {
		return second_right_or_inside_intersperses_first( prm_hit_lhs.get_segments(), prm_hit_rhs.get_segments() );
	}

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_CALC_HIT_HPP
