/// \file
/// \brief The seq_seg class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_HIT_SEG_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_HIT_SEG_H

#include <boost/core/ignore_unused.hpp>
#include <boost/operators.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "common/cpp14/constexpr_min_max.hpp"
#include "common/json_style.hpp"
#include "common/rapidjson_addenda/rapidjson_writer.hpp"
#include "common/size_t_literal.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "seq/seq_arrow.hpp"

using namespace cath::common::literals;

namespace cath {
	namespace seq {

		/// \brief Represent a single segment of a sequence hit (ie domain)
		///
		/// This can also be used to represent a fragment
		///
		/// Like all cath-resolve-hits code, this assumes simple residue numbering
		/// and is hence unsuitable for use with raw PDB residue numbers
		class seq_seg final : private boost::equality_comparable<seq_seg> {
		private:
			/// \brief The start boundary of the segment
			seq_arrow start;

			/// \brief The stop boundary of the segment
			seq_arrow stop;

			static constexpr bool sanity_check(const seq_arrow &,
			                                   const seq_arrow &);

		public:
			constexpr seq_seg(const seq_arrow &,
			                  const seq_arrow &);

			constexpr seq_seg(const residx_t &,
			                  const residx_t &);

			constexpr const seq_arrow & get_start_arrow() const;
			constexpr const seq_arrow & get_stop_arrow() const;

			seq_seg & set_start_arrow(const seq_arrow &);
			seq_seg & set_stop_arrow(const seq_arrow &);

			static auto get_seq_seg_start_less() {
				return [] (const seq_seg &x, const seq_seg &y) {
					return ( x.get_start_arrow() < y.get_start_arrow() );
				};
			}
		};

		/// \brief Sanity check this seq_seg and throw an exception if a problem is detected
		///
		/// \todo Come GCC >= 5 (with relaxed constexpr), make this code nicer
		constexpr bool seq_seg::sanity_check(const seq_arrow &arg_start, ///< The start position
		                                     const seq_arrow &arg_stop   ///< The stop position
		                                     ) {
			return ( arg_start < arg_stop ) ? true
			                                : throw std::invalid_argument( "Cannot create seq_seg with start residue before the stop residue" );
		}

		/// \brief Ctor seq_seg from a start and stop
		///
		/// \todo GCC >= 5 (with relaxed constexpr), move the check out into a separate sanity_check() function
		inline constexpr seq_seg::seq_seg(const seq_arrow &arg_start, ///< The residue boundary of the segment's start
		                                  const seq_arrow &arg_stop   ///< The residue boundary of the segment's stop
		                                  ) : start{ arg_start },
		                                      stop {
		                                      	sanity_check( arg_start, arg_stop )
		                                      	? arg_stop
		                                      	: arg_stop
		                                      } {
		}

		/// \brief seq_seg ctor from standard start/stop indices
		inline constexpr seq_seg::seq_seg(const residx_t &arg_start_idx, ///< The residue index of the segment's start
		                                  const residx_t &arg_stop_idx   ///< The residue index of the segment's stop
		                                  ) : seq_seg{
		                                      	arrow_before_res( arg_start_idx ),
		                                      	arrow_after_res ( arg_stop_idx  )
		                                      } {
		}

		/// \brief Getter for the start boundary
		inline constexpr const seq_arrow & seq_seg::get_start_arrow() const {
			return start;
		}

		/// \brief Getter for the stop boundary
		inline constexpr const seq_arrow & seq_seg::get_stop_arrow() const {
			return stop;
		}

		/// \brief Setter for the start boundary
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr
		inline seq_seg & seq_seg::set_start_arrow(const seq_arrow &arg_start ///< The start boundary to set
		                                          ) {
			sanity_check( arg_start, stop );
			start = arg_start;
			return *this;
		}

		/// \brief Setter for the stop boundary
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr
		inline seq_seg & seq_seg::set_stop_arrow(const seq_arrow &arg_stop ///< The stop boundary to set
		                                         ) {
			sanity_check( start, arg_stop );
			stop = arg_stop;
			return *this;
		}

		/// \brief Get the start residue index of the specified segment
		///
		/// \relates seq_seg
		inline constexpr const residx_t & get_start_res_index(const seq_seg &arg_seq_seg ///< The seq_seg to query
		                                                      ) {
			return arg_seq_seg.get_start_arrow().res_after();
		}

		/// \brief Get the stop residue index of the specified segment
		///
		/// \relates seq_seg
		inline constexpr residx_t get_stop_res_index(const seq_seg &arg_seq_seg ///< The seq_seg to query
		                                             ) {
			return arg_seq_seg.get_stop_arrow().res_before();
		}

		/// \brief Get the length of the specified seq_seg
		///
		/// \relates seq_seg
		inline constexpr residx_t get_length(const seq_seg &arg_seq_seg ///< The seq_seg to query
		                                   ) {
			return (
				arg_seq_seg.get_stop_arrow ()
				-
				arg_seq_seg.get_start_arrow()
			);
		}

		/// \brief Get the total length of the specified seq_segs
		///
		/// \relates hit
		inline size_t get_total_length(const seq_seg_vec &arg_seq_segs ///< The hits to query
		                               ) {
			return boost::accumulate(
				arg_seq_segs | boost::adaptors::transformed( &cath::seq::get_length ),
				0_z
			);
		}

		/// \brief Build a seq_seg from a pair of start/stop residue indices
		///
		/// \relates seq_seg
		inline constexpr seq_seg seq_seg_of_res_idx_pair(const residx_residx_pair &arg_res_idx_pair ///< The segments start/stop residue indices
		                                                 ) {
			return {
				arg_res_idx_pair.first,
				arg_res_idx_pair.second
			};
		}

		/// \brief Return whether the two specified seq_segs are identical
		///
		/// \relates seq_seg
		inline constexpr bool operator==(const seq_seg &arg_seq_seg_a, ///< The first  seq_seg to compare
		                                 const seq_seg &arg_seq_seg_b  ///< The second seq_seg to compare
		                                 ) {
			return (
				arg_seq_seg_a.get_start_arrow() == arg_seq_seg_b.get_start_arrow()
				&&
				arg_seq_seg_a.get_stop_arrow()  == arg_seq_seg_b.get_stop_arrow()
			);
		}

		/// \brief Return whether the two specified seq_segs overlap
		///
		/// Note: don't call this `overlap` - that can cause problems with other `overlap` functions
		///
		/// \relates seq_seg
		inline constexpr bool are_overlapping(const seq_seg &arg_seq_seg_a, ///< The first  seq_seg to query
		                                      const seq_seg &arg_seq_seg_b  ///< The second seq_seg to query
		                                      ) {
			return (
				arg_seq_seg_a.get_start_arrow() <  arg_seq_seg_b.get_stop_arrow()
				&&
				arg_seq_seg_b.get_start_arrow() <  arg_seq_seg_a.get_stop_arrow()
			);
		}

		/// \brief Return the number of residues by which the two specified seq_segs overlap (or 0 if they don't overlap)
		///
		/// \todo Come relaxed, constexpr use this (but keep constexpr_max() / constexpr_min() if they haven't been retired yet):
		///     const auto max_start = std::max( arg_seq_seg_a.get_start_arrow(), arg_seq_seg_b.get_start_arrow() );
		///     const auto min_stop  = std::min( arg_seq_seg_a.get_stop_arrow(),  arg_seq_seg_b.get_stop_arrow()  );
		///     return ( std::max( min_stop, max_start ) - max_start );
		///
		/// \relates seq_seg
		inline constexpr residx_t overlap_by(const seq_seg &arg_seq_seg_a, ///< The first  seq_seg to query
		                                     const seq_seg &arg_seq_seg_b  ///< The second seq_seg to query
		                                     ) {
			return (
				common::constexpr_max(
					common::constexpr_min( arg_seq_seg_a.get_stop_arrow(),  arg_seq_seg_b.get_stop_arrow()  ),
					common::constexpr_max( arg_seq_seg_a.get_start_arrow(), arg_seq_seg_b.get_start_arrow() )
				)
				-
				common::constexpr_max( arg_seq_seg_a.get_start_arrow(), arg_seq_seg_b.get_start_arrow() )
			);
		}

		/// \brief Return the length of the shorter of the two specified seq_segs
		///
		/// \relates seq_seg
		inline constexpr residx_t shorter_length(const seq_seg &arg_seq_seg_a, ///< The first  seq_seg to query
		                                         const seq_seg &arg_seq_seg_b  ///< The second seq_seg to query
		                                         ) {
			return static_cast<residx_t>( common::constexpr_min(
				get_length( arg_seq_seg_a ),
				get_length( arg_seq_seg_b )
			) );
		}

		/// \brief Return the length of the longer of the two specified seq_segs
		///
		/// \relates seq_seg
		inline constexpr residx_t longer_length(const seq_seg &arg_seq_seg_a, ///< The first  seq_seg to query
		                                        const seq_seg &arg_seq_seg_b  ///< The second seq_seg to query
		                                        ) {
			return static_cast<residx_t>( common::constexpr_max(
				get_length( arg_seq_seg_a ),
				get_length( arg_seq_seg_b )
			) );
		}

		/// \brief Return the fraction overlap between the two specified seq_segs over the length of the shorter
		///
		/// \relates seq_seg
		inline constexpr double fraction_overlap_over_shorter(const seq_seg &arg_seq_seg_a, ///< The first  seq_seg to query
		                                                      const seq_seg &arg_seq_seg_b  ///< The second seq_seg to query
		                                                      ) {
			return
				static_cast<double>( overlap_by   ( arg_seq_seg_a, arg_seq_seg_b ) )
				/
				static_cast<double>( shorter_length( arg_seq_seg_a, arg_seq_seg_b ) );
		}

		/// \brief Return the fraction overlap between the two specified seq_segs over the length of the longer
		///
		/// \relates seq_seg
		inline constexpr double fraction_overlap_over_longer(const seq_seg &arg_seq_seg_a, ///< The first  seq_seg to query
		                                                     const seq_seg &arg_seq_seg_b  ///< The second seq_seg to query
		                                                     ) {
			return
				static_cast<double>( overlap_by   ( arg_seq_seg_a, arg_seq_seg_b ) )
				/
				static_cast<double>( longer_length( arg_seq_seg_a, arg_seq_seg_b ) );
		}

		/// \brief In-place sort the specified segments by their starts
		///
		/// \relates seq_seg
		inline void start_sort_seq_segs(seq_seg_vec &arg_seq_segs ///< The segments to sort
		                                ) {
			boost::range::sort(
				arg_seq_segs,
				seq_seg::get_seq_seg_start_less()
			);
		}

		/// \brief Return a copy of the specified segments, sorted by their starts
		///
		/// \relates seq_seg
		inline seq_seg_vec start_sort_seq_segs_copy(seq_seg_vec arg_seq_segs ///< The segments to sort
		                                            ) {
			start_sort_seq_segs( arg_seq_segs );
			return arg_seq_segs;
		}

		/// \brief Build a vector of seq_seg from the specified list of bounds representing start/stop residue index pairs
		///
		/// \pre `arg_bounds.size() % 2 != 0` else an invalid_argument_exception will be thrown
		///
		/// \pre Each start after the first must be on a residue later than the preceding stop residue
		///       else an invalid_argument_exception will be thrown
		///
		/// \pre Each seq_seq must be valid (ie stop must be on a residue not earlier than the corresponding start residue)
		///       else an invalid_argument will be thrown
		///
		/// \relates seq_seg
		inline seq_seg_vec segments_from_bounds(const residx_vec &arg_bounds ///< The start/stop residue index pairs from which to build the seq_segs
		                                        ) {
			// Check there's an even number of bounds
			if ( arg_bounds.size() % 2 != 0 ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("There should be an even number of bounds from which to make segments"));
			}

			// Prepare a vector of segments of the correct size
			seq_seg_vec segments;
			segments.reserve( arg_bounds.size() / 2 );

			// Loop over the pairs of bounds
			for (const size_t &bound_ctr : boost::irange( 0_z, arg_bounds.size(), 2_z ) ) {
				const auto start_arrow = arrow_before_res( arg_bounds[ bound_ctr     ] );
				const auto stop_arrow  = arrow_after_res ( arg_bounds[ bound_ctr + 1 ] );

				// Check the segment comes after any preceding segments, else throw
				if ( ! segments.empty() && segments.back().get_stop_arrow() > start_arrow ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Whilst building segments from bounds, found preceding stop after start"));
				}

				// Add the seg_seg (whose ctor will check start < stop)
				segments.emplace_back( start_arrow, stop_arrow );
			}

			return segments;
		}

		/// \brief Whether the segments in the first specified seq_seg_run never extend outside
		///        those in the second specified seq_seg_run
		///
		/// In other words: the first is within or equal to the second
		///
		/// \relates seq_seg
		inline bool first_is_not_outside_second(const seq_seg_vec &arg_seq_segs_a, ///< The first  seq_seg_vec to query
		                                        const seq_seg_vec &arg_seq_segs_b  ///< The second seq_seg_vec to query
		                                        ) {
			if ( arg_seq_segs_a.front().get_start_arrow() < arg_seq_segs_b.front().get_start_arrow() ) {
				return false;
			}
			if ( arg_seq_segs_a.back().get_stop_arrow  () > arg_seq_segs_b.back().get_stop_arrow  () ) {
				return false;
			}

			const size_t num_segments_lhs = arg_seq_segs_a.size();
			const size_t num_segments_rhs = arg_seq_segs_b.size();

			size_t rhs_ctr = 0;
			for (const auto &lhs_ctr : boost::irange( 0_z, num_segments_lhs ) ) {
				while ( rhs_ctr < num_segments_rhs && arg_seq_segs_a[ lhs_ctr ].get_stop_arrow() > arg_seq_segs_b[ rhs_ctr ].get_stop_arrow() ) {
					++rhs_ctr;
				}
				if ( rhs_ctr == num_segments_rhs ) {
					return false;
				}
				if ( arg_seq_segs_a[ lhs_ctr ].get_start_arrow() < arg_seq_segs_b[ rhs_ctr ].get_start_arrow() ) {
					return false;
				}
			}
			return true;
		}

		/// \brief Whether either of the specified seq_seg_run covers the other
		///
		/// \relates seq_seg
		inline bool one_covers_other(const seq_seg_vec &arg_seq_segs_a, ///< The first  seq_seg_vec to query
		                             const seq_seg_vec &arg_seq_segs_b  ///< The second seq_seg_vec to query
		                             ) {
			const residx_t length_a = static_cast<residx_t>( get_total_length( arg_seq_segs_a ) );
			const residx_t length_b = static_cast<residx_t>( get_total_length( arg_seq_segs_b ) );
			if ( length_a < length_b ) {
				return first_is_not_outside_second( arg_seq_segs_a, arg_seq_segs_b );
			}
			if ( length_a > length_b ) {
				return first_is_not_outside_second( arg_seq_segs_b, arg_seq_segs_a );
			}
			return ( arg_seq_segs_a == arg_seq_segs_b );
		}

		/// \brief Whether the segments in the first specified seq_seg_vec are shorter strictly
		///        within those  in the second specified seq_seg_vec
		///
		/// In other words: the first is within and not equal to the second
		///
		/// \relates seq_seg
		inline bool first_is_shorter_and_within_second(const seq_seg_vec &arg_seq_segs_a, ///< The first  calc_hit to query
		                                               const seq_seg_vec &arg_seq_segs_b  ///< The second calc_hit to query
		                                               ) {
			return (
				get_total_length( arg_seq_segs_a ) < get_total_length( arg_seq_segs_b )
				&&
				first_is_not_outside_second( arg_seq_segs_a, arg_seq_segs_b )
			);
		}

		/// \brief Calculate a hash number for the segments in the specified seq_seg_vec
		///
		/// \relates seq_seg
		inline size_t calc_hash(const seq_seg_vec &arg_seq_segs ///< The segments to hash
		                        ) {
			if ( arg_seq_segs.empty() ) {
				return 0;
			}
			const std::hash<resarw_t> hasher{};
			size_t result = hasher( arg_seq_segs.front().get_start_arrow().get_index() );
			const auto combine_fn = [&] (const resarw_t &x) {
				result ^= hasher( x ) + 0x9e3779b9 + ( result << 6 ) + ( result >> 2 );
			};
			for (const size_t &seg_ctr : boost::irange( 0_z, arg_seq_segs.size() ) ) {
				combine_fn( arg_seq_segs[ seg_ctr ].get_start_arrow().get_index() );
				combine_fn( arg_seq_segs[ seg_ctr ].get_stop_arrow() .get_index() );
			}
			combine_fn( arg_seq_segs.back().get_stop_arrow().get_index() );
			return result;
		}

		seq_seg_vec get_present_segments(const seq_seg_opt_vec &);
		std::string get_segments_string(const seq_seg_vec &);
		seq_seg_vec make_fragments_of_segments(seq_seg_vec);
		bool segments_are_start_sorted_and_non_overlapping(const seq_seg_vec &);
		seq_seg_vec make_fragments_of_start_sorted_segments(const seq_seg_vec &);
		bool midpoint_less(const seq_seg &,
		                   const seq_seg &);
		std::string to_simple_string(const seq_seg &);
		std::string to_simple_string(const seq_seg_opt &);
		std::string to_simple_seg_string(const seq_arrow &,
		                                 const seq_arrow &);
		std::string to_string(const seq_seg &);
		std::ostream & operator<<(std::ostream &,
		                          const seq_seg &);

		/// \brief Write the specified seq_seg to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson(common::rapidjson_writer<Style> &arg_writer, ///< The rapidjson_writer to which the seq_seg should be written
		                        const seq_seg                   &arg_seq_seg ///< The seq_seg to write
		                        ) {
			arg_writer.start_array();
			arg_writer.write_value( get_start_res_index( arg_seq_seg ) );
			arg_writer.write_value( get_stop_res_index ( arg_seq_seg ) );
			arg_writer.end_array();
		}

		/// \brief Write the specified seq_seg_vec to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson(common::rapidjson_writer<Style> &arg_writer,  ///< The rapidjson_writer to which the seq_seg_vec should be written
		                        const seq_seg_vec               &arg_seq_segs ///< The seq_seg_vec to write
		                        ) {
			arg_writer.start_array();
			for (const seq_seg &the_seq_seg : arg_seq_segs) {
				write_to_rapidjson( arg_writer, the_seq_seg );
			}
			arg_writer.end_array();
		}

	} // namespace seq
} // namespace cath

#endif
