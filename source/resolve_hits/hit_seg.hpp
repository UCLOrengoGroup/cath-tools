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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_HIT_SEG_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_HIT_SEG_H

#include <boost/core/ignore_unused.hpp>
#include <boost/operators.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/numeric.hpp>

#include "common/json_style.hpp"
#include "common/rapidjson_addenda/rapidjson_writer.hpp"
#include "common/size_t_literal.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "resolve_hits/res_arrow.hpp"

using namespace cath::common::literals;

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

			static constexpr bool sanity_check(const res_arrow &,
			                                   const res_arrow &);

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

		/// \brief Sanity check this hit_seg and throw an exception if a problem is detected
		///
		/// \todo Come GCC >= 5 (with relaxed constexpr), make this code nicer
		constexpr bool hit_seg::sanity_check(const res_arrow &arg_start, ///< The start position
		                                     const res_arrow &arg_stop   ///< The stop position
		                                     ) {
			return ( arg_start < arg_stop ) ? true
			                                : throw std::invalid_argument( "Cannot create hit_seg with start residue before the stop residue" );
		}

		/// \brief Ctor hit_seg from a start and stop
		///
		/// \todo GCC >= 5 (with relaxed constexpr), move the check out into a separate sanity_check() function
		inline constexpr hit_seg::hit_seg(const res_arrow &arg_start, ///< The residue boundary of the segment's start
		                                  const res_arrow &arg_stop   ///< The residue boundary of the segment's stop
		                                  ) : start{ arg_start },
		                                      stop {
		                                      	sanity_check( arg_start, arg_stop )
		                                      	? arg_stop
		                                      	: arg_stop
		                                      } {
		}

		/// \brief Getter for the start boundary
		inline constexpr const res_arrow & hit_seg::get_start_arrow() const {
			return start;
		}

		/// \brief Getter for the stop boundary
		inline constexpr const res_arrow & hit_seg::get_stop_arrow() const {
			return stop;
		}

		/// \brief Setter for the start boundary
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr
		inline hit_seg & hit_seg::set_start_arrow(const res_arrow &arg_start ///< The start boundary to set
		                                          ) {
			sanity_check( arg_start, stop );
			start = arg_start;
			return *this;
		}

		/// \brief Setter for the stop boundary
		///
		/// \todo GCC >= 5 (with relaxed constexpr), make this constexpr
		inline hit_seg & hit_seg::set_stop_arrow(const res_arrow &arg_stop ///< The stop boundary to set
		                                         ) {
			sanity_check( start, arg_stop );
			stop = arg_stop;
			return *this;
		}

		/// \brief Get the start residue index of the specified segment
		///
		/// \relates hit_seg
		inline constexpr const residx_t & get_start_res_index(const hit_seg &arg_hit_seg ///< The hit_seg to query
		                                                      ) {
			return arg_hit_seg.get_start_arrow().res_after();
		}

		/// \brief Get the stop residue index of the specified segment
		///
		/// \relates hit_seg
		inline constexpr residx_t get_stop_res_index(const hit_seg &arg_hit_seg ///< The hit_seg to query
		                                             ) {
			return arg_hit_seg.get_stop_arrow().res_before();
		}

		/// \brief Get the length of the specified hit_seg
		///
		/// \relates hit_seg
		inline constexpr residx_t get_length(const hit_seg &arg_hit_seg ///< The hit_seg to query
		                                   ) {
			return (
				arg_hit_seg.get_stop_arrow ()
				-
				arg_hit_seg.get_start_arrow()
			);
		}

		/// \brief Get the total length of the specified hit_segs
		///
		/// \relates hit
		inline size_t get_total_length(const hit_seg_vec &arg_hit_segs ///< The hits to query
		                               ) {
			return boost::accumulate(
				arg_hit_segs | boost::adaptors::transformed( &cath::rslv::get_length ),
				0_z
			);
		}

		/// \brief Build a hit_seg of the specified start/stop residue indices
		///
		/// \relates hit_seg
		inline constexpr hit_seg hit_seg_of_res_idcs(const residx_t &arg_start_res_idx, ///< The segment's start residue index
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
		inline constexpr hit_seg hit_seg_of_res_idx_pair(const residx_residx_pair &arg_res_idx_pair ///< The segments start/stop residue indices
		                                                 ) {
			return {
				arrow_before_res( arg_res_idx_pair.first  ),
				arrow_after_res ( arg_res_idx_pair.second )
			};
		}

		/// \brief Return whether the two specified hit_segs are identical
		///
		/// \relates hit_seg
		inline constexpr bool operator==(const hit_seg &arg_hit_seg_a, ///< The first  hit_seg to compare
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
		inline constexpr bool hit_segs_overlap(const hit_seg &arg_hit_seg_a, ///< The first  hit_seg to query
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
		bool midpoint_less(const hit_seg &,
		                   const hit_seg &);
		std::string to_simple_string(const hit_seg &);
		std::string to_simple_string(const hit_seg_opt &);
		std::string to_simple_seg_string(const res_arrow &,
		                                 const res_arrow &);
		std::string to_string(const hit_seg &);
		std::ostream & operator<<(std::ostream &,
		                          const hit_seg &);

		/// \brief Write the specified hit_seg to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson(common::rapidjson_writer<Style> &arg_writer, ///< The rapidjson_writer to which the hit_seg should be written
		                        const hit_seg                   &arg_hit_seg ///< The hit_seg to write
		                        ) {
			arg_writer.start_array();
			arg_writer.write_uint( get_start_res_index( arg_hit_seg ) );
			arg_writer.write_uint( get_stop_res_index ( arg_hit_seg ) );
			arg_writer.end_array();
		}

		/// \brief Write the specified hit_seg_vec to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson(common::rapidjson_writer<Style> &arg_writer,  ///< The rapidjson_writer to which the hit_seg_vec should be written
		                        const hit_seg_vec               &arg_hit_segs ///< The hit_seg_vec to write
		                        ) {
			arg_writer.start_array();
			for (const hit_seg &the_hit_seg : arg_hit_segs) {
				write_to_rapidjson( arg_writer, the_hit_seg );
			}
			arg_writer.end_array();
		}

	} // namespace rslv
} // namespace cath

#endif
