/// \file
/// \brief The seq_seg_run class header

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

#ifndef _CATH_TOOLS_SOURCE_SEQ_SEQ_SEG_RUN_H
#define _CATH_TOOLS_SOURCE_SEQ_SEQ_SEG_RUN_H

#include <boost/range/irange.hpp>

#include "common/algorithm/append.hpp"
#include "common/algorithm/transform_build.hpp"
#include "seq/seq_arrow.hpp"
#include "seq/seq_seg.hpp"

namespace cath {
	namespace seq {

		class seq_seg_run;
		inline seq_arrow get_stop_of_first_segment(const seq_seg_run &);
		inline seq_arrow get_start_of_last_segment(const seq_seg_run &);

		/// \brief TODOCUMENT
		class seq_seg_run final {
		private:
			/// \brief The boundary at the start of the first segment
			seq_arrow start_arrow;

			/// \brief The boundary at the end of the last segment
			seq_arrow stop_arrow;

			/// \brief The (possibly empty) list of the boundaries associated with any gaps between this seq_seg_run's segments
			seq_seg_vec fragments;

			void sanity_check() const;

		public:
			seq_seg_run(seq_arrow,
			            seq_arrow);

			seq_seg_run(const seq_seg_vec &);

			seq_seg_run(seq_arrow,
			            seq_arrow,
			            seq_seg_vec);

			bool is_discontig() const;
			size_t get_num_segments() const;
			const seq_arrow & get_start_arrow_of_segment(const size_t &) const;
			const seq_arrow & get_stop_arrow_of_segment(const size_t &) const;

			const seq_arrow & get_start_arrow() const;
			const seq_arrow & get_stop_arrow () const;

			static auto get_seg_run_start_less() {
				return [] (const seq_seg_run &x, const seq_seg_run &y) {
					return ( x.get_start_arrow() < y.get_start_arrow() );
				};
			}

			static auto get_seg_run_stop_less() {
				return [] (const seq_seg_run &x, const seq_seg_run &y) {
					return ( x.get_stop_arrow() < y.get_stop_arrow() );
				};
			}

			static auto get_seg_run_first_seg_stop_less() {
				return [] (const seq_seg_run &x, const seq_seg_run &y) {
					return ( get_stop_of_first_segment( x ) < get_stop_of_first_segment( y ) );
				};
			}

			static auto get_seg_run_last_seg_start_less() {
				return [] (const seq_seg_run &x, const seq_seg_run &y) {
					return ( get_start_of_last_segment( x ) < get_start_of_last_segment( y ) );
				};
			}
		};

		std::string get_segments_string(const seq_seg_run &);
		std::string to_string(const seq_seg_run &);
		std::ostream & operator<<(std::ostream &,
		                          const seq_seg_run &);
		bool operator==(const seq_seg_run &,
		                const seq_seg_run &);
		bool any_interaction(const seq_seg_run &,
		                     const seq_seg_run &);

		/// \brief Sanity check that the seq_seg_run is sensible and throw an exception if not
		inline void seq_seg_run::sanity_check() const {
			if ( stop_arrow < start_arrow ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Start index must not be greater than the stop index"));
			}
		}

		/// \brief Ctor for contiguous seq_seg_run
		inline seq_seg_run::seq_seg_run(seq_arrow   arg_start_arrow, ///< The start boundary of the continuous seq_seg_run
		                                seq_arrow   arg_stop_arrow   ///< The end boundary of the continuous seq_seg_run
		                                ) : start_arrow ( std::move( arg_start_arrow ) ),
		                                    stop_arrow  ( std::move( arg_stop_arrow  ) ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous seq_seg_run from segments
		inline seq_seg_run::seq_seg_run(const seq_seg_vec &arg_segments ///< The segments of the seq_seg_run
		                                ) : start_arrow ( arg_segments.front().get_start_arrow()      ),
		                                    stop_arrow  ( arg_segments.back ().get_stop_arrow ()      ),
		                                    fragments   ( make_fragments_of_segments( arg_segments )  ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous seq_seg_run from start, stop and fragments
		inline seq_seg_run::seq_seg_run(seq_arrow     arg_start_arrow, ///< The boundary at the start of the first segment
		                                seq_arrow     arg_stop_arrow,  ///< The boundary at the end of the last segment
		                                seq_seg_vec   arg_fragments    ///< The (possibly empty) list of the boundaries associated with any gaps between this seq_seg_run's segments
		                                ) : start_arrow ( std::move( arg_start_arrow        ) ),
		                                    stop_arrow  ( std::move( arg_stop_arrow         ) ),
		                                    fragments   ( std::move( arg_fragments          ) ) {
			sanity_check();
		}

		/// \brief Return whether this seq_seg_run is discontiguous
		inline bool seq_seg_run::is_discontig() const {
			return ! fragments.empty();
		}

		/// \brief Return the number of segments in this seq_seg_run
		inline size_t seq_seg_run::get_num_segments() const {
			return fragments.size() + 1;
		}

		/// \brief Get the start boundary of the segment with the specified index
		inline const seq_arrow & seq_seg_run::get_start_arrow_of_segment(const size_t &arg_segment_index ///< The index of the segment whose start arrow should be returned
		                                                                 ) const {
			return ( arg_segment_index > 0                ) ? fragments[ arg_segment_index - 1 ].get_stop_arrow()
			                                                : start_arrow;
		}

		/// \brief Get the stop boundary of the segment with the specified index
		inline const seq_arrow & seq_seg_run::get_stop_arrow_of_segment(const size_t &arg_segment_index ///< The index of the segment whose stop arrow should be returned
		                                                                ) const {
			return ( arg_segment_index < fragments.size() ) ? fragments[ arg_segment_index     ].get_start_arrow()
			                                                : stop_arrow;
		}

		/// \brief Get the length of the specified seq_seg_run's segment corresponding to the specified index
		///
		/// \relates seq_seg_run
		inline size_t get_length_of_seq_seg(const seq_seg_run &arg_seq_seg_run, ///< The seq_seg_run to query
		                                    const size_t      &arg_seg_idx      ///< The index of the segment who length should be returned
		                                    ) {
			return static_cast<size_t>(
				arg_seq_seg_run.get_stop_arrow_of_segment ( arg_seg_idx )
				-
				arg_seq_seg_run.get_start_arrow_of_segment( arg_seg_idx )
			);
		}
		
		/// \brief Get the specified seq_seg_run's segment corresponding to the specified index
		///
		/// \relates seq_seg_run
		inline seq_seg get_seq_seg_of_seg_idx(const seq_seg_run &arg_seq_seg_run, ///< The seq_seg_run to query
		                                      const size_t      &arg_seg_idx      ///< The index of the segment to return
		                                      ) {
			return {
				arg_seq_seg_run.get_start_arrow_of_segment( arg_seg_idx ),
				arg_seq_seg_run.get_stop_arrow_of_segment ( arg_seg_idx )
			};
		}

		/// \brief Get a vector of the segments in this seq_seg_run
		///
		/// \relates seq_seg_run
		inline const seq_seg_vec get_seq_segs(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to query
		                                      ) {
			return common::transform_build<seq_seg_vec>(
				boost::irange( 0_z, arg_seq_seg_run.get_num_segments() ),
				[&] (const size_t &x) {
					return get_seq_seg_of_seg_idx( arg_seq_seg_run, x );
				}
			);
		}

		/// \brief Get the (possibly-repeated, non-sorted) segments from the specified seq_seg_runs
		///
		/// \relates seq_seg_run
		inline const seq_seg_vec get_seq_segs(const seq_seg_run_vec &arg_seq_seg_run_vec ///< The seq_seg_runs whose segments should be returned
		                                      ) {
			seq_seg_vec results;
			for (const seq_seg_run &the_seq_seg_run : arg_seq_seg_run_vec) {
				common::append( results, get_seq_segs( the_seq_seg_run ) );
			}
			return results;
		}

		/// \brief Get a vector of the specified seq_seg_runs' segments, sorted by their starts
		inline seq_seg_vec get_start_sorted_seq_segs(const seq_seg_run_vec &arg_seq_seg_run_vec ///< The vector of seq_seg_runs to query
		                                             ) {
			return start_sort_seq_segs_copy( get_seq_segs( arg_seq_seg_run_vec ) );
		}

		/// \brief Get the (first) start of this seq_seg_run
		inline const seq_arrow & seq_seg_run::get_start_arrow() const {
			return start_arrow;
		}

		/// \brief Get the (last) stop of this seq_seg_run
		inline const seq_arrow & seq_seg_run::get_stop_arrow()  const {
			return stop_arrow;
		}

		/// \brief Get the start residue index of the segment of specified index in the specified seq_seg_run
		///
		/// \relates seq_seg_run
		inline const residx_t & get_start_res_index_of_segment(const seq_seg_run &arg_seq_seg_run,  ///< The seq_seg_run to query
		                                                       const size_t      &arg_segment_index ///< The index of the segment to query
		                                                       ) {
			return arg_seq_seg_run.get_start_arrow_of_segment( arg_segment_index ).res_after();
		}

		/// \brief Get the stop residue index of the segment of specified index in the specified seq_seg_run
		///
		/// \relates seq_seg_run
		inline residx_t get_stop_res_index_of_segment(const seq_seg_run &arg_seq_seg_run,  ///< The seq_seg_run to query
		                                              const size_t      &arg_segment_index ///< The index of the segment to query
		                                              ) {
			return arg_seq_seg_run.get_stop_arrow_of_segment( arg_segment_index ).res_before();
		}

		/// \brief Get the start residue index of the specified seq_seg_run
		///
		/// \relates seq_seg_run
		inline const residx_t & get_start_res_index(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to query
		                                            ) {
			return arg_seq_seg_run.get_start_arrow().res_after();
		}

		/// \brief Get the stop residue index of the specified seq_seg_run
		///
		/// \relates seq_seg_run
		inline residx_t get_stop_res_index(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to query
		                                   ) {
			return arg_seq_seg_run.get_stop_arrow().res_before();
		}

		/// \brief Get the stop of the first segment in the specified seq_seg_run
		///
		/// \pre `arg_seq_seg_run.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates seq_seg_run
		inline seq_arrow get_stop_of_first_segment(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to query
		                                           ) {
			if ( ! arg_seq_seg_run.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_stop_of_first_segment of contiguous seq_seg_run"));
			}
			return arg_seq_seg_run.get_stop_arrow_of_segment( 0 );
		}

		/// \brief Get the start of the last segment in the specified seq_seg_run
		///
		/// \pre `arg_seq_seg_run.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates seq_seg_run
		inline seq_arrow get_start_of_last_segment(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to query
		                                           ) {
			if ( ! arg_seq_seg_run.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_start_of_last_segment of contiguous seq_seg_run"));
			}
			return arg_seq_seg_run.get_start_arrow_of_segment( arg_seq_seg_run.get_num_segments() - 1 );
		}

		/// \brief Get the total length of the specified seq_seg_run (ie the sum of its segments' lengths)
		///
		/// \relates seq_seg_run
		inline size_t get_total_length(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to query
		                               ) {
			return boost::accumulate(
				boost::irange( 0_z, arg_seq_seg_run.get_num_segments() )
					| boost::adaptors::transformed( [&] (const size_t &x) {
						return get_length_of_seq_seg( arg_seq_seg_run, x );
					} ),
				0_z
			);
		}

		/// \brief Make a continuous seq_seg_run from the residue indices
		///
		/// \relates seq_seg_run
		inline seq_seg_run make_seq_seg_run_from_res_indices(const residx_t &arg_start_res_idx, ///< The start residue index
		                                                     const residx_t &arg_stop_res_idx   ///< The stop residue index
		                                                     ) {
			return {
				arrow_before_res( arg_start_res_idx ),
				arrow_after_res ( arg_stop_res_idx  )
			};
		}

		/// \brief Make a seq_seg_run
		///
		/// \relates seq_seg_run
		inline seq_seg_run make_seq_seg_run_from_res_indices(const residx_residx_pair_vec &arg_residue_index_segments ///< The residue index start/stop pairs of the seq_seg_run's segments
		                                                     ) {
			return {
				common::transform_build<seq_seg_vec>(
					arg_residue_index_segments,
					seq_seg_of_res_idx_pair
				)
			};
		}

		/// \brief Return whether the either of the two specified seq_seg_runs overlaps, interleaves or straddles the other
		///
		/// \relates seq_seg_run
		inline bool any_interaction(const seq_seg_run &arg_seq_seg_run_a, ///< The first  seq_seg_run to query
		                            const seq_seg_run &arg_seq_seg_run_b  ///< The second seq_seg_run to query
		                            ) {
			return (
				arg_seq_seg_run_a.get_start_arrow() < arg_seq_seg_run_b.get_stop_arrow()
				&&
				arg_seq_seg_run_b.get_start_arrow() < arg_seq_seg_run_a.get_stop_arrow()
			);
		}
		
		/// \brief Return whether the two specified seq_seg_runs overlap with each other
		///
		/// This requires there to be a genuine overlap of segments, not just that one
		/// seq_seg_run interleaves or straddles the other
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \todo This can be made more efficient by iterating over the two segment
		///       lists simultaneously
		///
		/// \relates seq_seg_run
		inline bool are_overlapping(const seq_seg_run &arg_seq_seg_run_a, ///< The first  seq_seg_run to query
		                            const seq_seg_run &arg_seq_seg_run_b  ///< The second seq_seg_run to query
		                            ) {
			if ( ! any_interaction( arg_seq_seg_run_a, arg_seq_seg_run_b ) ) {
				return false;
			}
			for (const auto &seg_ctr_a : boost::irange( 0_z, arg_seq_seg_run_a.get_num_segments() ) ) {
				for (const auto &seg_ctr_b : boost::irange( 0_z, arg_seq_seg_run_b.get_num_segments() ) ) {
					const bool seg_overlap = are_overlapping(
						get_seq_seg_of_seg_idx( arg_seq_seg_run_a, seg_ctr_a ),
						get_seq_seg_of_seg_idx( arg_seq_seg_run_b, seg_ctr_b )
					);
					if ( seg_overlap ) {
						return true;
					}
				}
			}
			return false;
		}

	} // namespace seq
} // namespace cath

#endif
