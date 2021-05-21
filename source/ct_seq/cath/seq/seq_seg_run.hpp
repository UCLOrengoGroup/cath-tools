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

#ifndef _CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_SEG_RUN_HPP
#define _CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_SEG_RUN_HPP

#include <array>
#include <functional>

#include "cath/common/algorithm/append.hpp"
#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/algorithm/variadic_and.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/hash/hash_value_combine.hpp"
#include "cath/common/type_traits.hpp"
#include "cath/seq/seq_arrow.hpp"
#include "cath/seq/seq_seg.hpp"

namespace cath {
	namespace seq {

		class seq_seg_run;
		inline seq_arrow get_stop_of_first_segment(const seq_seg_run &);
		inline seq_arrow get_start_of_last_segment(const seq_seg_run &);

		/// \brief Represent a series of non-overlapping, increasing segments
		///
		/// This stores the first start and last stop on the stack and any segments as
		/// gaps so that the first-start/last stop can be processed (without accessing
		/// external memory and hence) very quickly.
		///
		/// Many seq_seg_runs are single-segment, which can be handled completely locally
		///
		/// \todo Carefully review handling of consecutive segments that touch each other,
		///       particularly in the associated non-member, non-friend functions
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

			explicit seq_seg_run(const seq_seg_vec &);

			seq_seg_run(seq_arrow,
			            seq_arrow,
			            seq_seg_vec);

			[[nodiscard]] bool             is_discontig() const;
			[[nodiscard]] size_t           get_num_segments() const;
			[[nodiscard]] const seq_arrow &get_start_arrow_of_segment( const size_t & ) const;
			[[nodiscard]] const seq_arrow &get_stop_arrow_of_segment( const size_t & ) const;

			[[nodiscard]] const seq_arrow &get_start_arrow() const;
			[[nodiscard]] const seq_arrow &get_stop_arrow() const;

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

		/// \brief Calculate a hash number for the segments in the seq_seg_run
		inline size_t calc_hash(const seq_seg_run &prm_seq_seg_run ///< The segments to hash
		                        ) {
			const std::hash<resarw_t> hasher{};
			size_t result = hasher( prm_seq_seg_run.get_start_arrow().get_index() );
			const auto combine_fn = [&] (const resarw_t &x) {
				common::hash_value_combine( result, hasher( x ) );
			};
			for (const size_t &seg_ctr : common::indices( prm_seq_seg_run.get_num_segments() ) ) {
				combine_fn( prm_seq_seg_run.get_start_arrow_of_segment( seg_ctr ).get_index() );
				combine_fn( prm_seq_seg_run.get_stop_arrow_of_segment ( seg_ctr ).get_index() );
			}
			combine_fn( prm_seq_seg_run.get_stop_arrow().get_index() );
			return result;
		}

		std::string get_segments_string(const seq_seg_run &);
		std::string get_segments_suffix_string(const seq_seg_run_opt &);
		std::string to_string(const seq_seg_run &);
		std::string to_string(const seq_seg_run_opt &);
		std::ostream & operator<<(std::ostream &,
		                          const seq_seg_run &);
		bool operator==(const seq_seg_run &,
		                const seq_seg_run &);
		bool any_interaction(const seq_seg_run &,
		                     const seq_seg_run &);

		/// \brief Sanity check that the seq_seg_run is sensible and throw an exception if not
		inline void seq_seg_run::sanity_check() const {
			if ( stop_arrow <= start_arrow ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Start index must not be greater than or equal to the stop index"));
			}
			if ( ! fragments.empty() ) {
				if ( start_arrow >= fragments.front().get_start_arrow() ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create a seq_seg_run with fragments that don't start after the start"));
				}
				if ( stop_arrow  <= fragments.back ().get_stop_arrow () ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create a seq_seg_run with fragments that don't end before the end"));
				}
				const auto first_seq_seg_is_not_earlier = [] (const seq_seg &x, const seq_seg &y) {
					return x.get_stop_arrow() >= y.get_start_arrow();
				};
				if ( common::contains_adjacent_match( fragments, first_seq_seg_is_not_earlier ) ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create a seq_seg_run with fragments that aren't increasing"));
				}
			}
		}

		/// \brief Ctor for contiguous seq_seg_run
		inline seq_seg_run::seq_seg_run(seq_arrow   prm_start_arrow, ///< The start boundary of the continuous seq_seg_run
		                                seq_arrow   prm_stop_arrow   ///< The end boundary of the continuous seq_seg_run
		                                ) : start_arrow ( std::move( prm_start_arrow ) ),
		                                    stop_arrow  ( std::move( prm_stop_arrow  ) ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous seq_seg_run from segments
		inline seq_seg_run::seq_seg_run(const seq_seg_vec &prm_segments ///< The segments of the seq_seg_run
		                                ) : start_arrow ( prm_segments.front().get_start_arrow()      ),
		                                    stop_arrow  ( prm_segments.back ().get_stop_arrow ()      ),
		                                    fragments   ( make_fragments_of_segments( prm_segments )  ) {
			sanity_check();
		}

		/// \brief Ctor for a possibly discontinuous seq_seg_run from start, stop and fragments
		inline seq_seg_run::seq_seg_run(seq_arrow     prm_start_arrow, ///< The boundary at the start of the first segment
		                                seq_arrow     prm_stop_arrow,  ///< The boundary at the end of the last segment
		                                seq_seg_vec   prm_fragments    ///< The (possibly empty) list of the boundaries associated with any gaps between this seq_seg_run's segments
		                                ) : start_arrow ( std::move( prm_start_arrow        ) ),
		                                    stop_arrow  ( std::move( prm_stop_arrow         ) ),
		                                    fragments   ( std::move( prm_fragments          ) ) {
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
		inline const seq_arrow & seq_seg_run::get_start_arrow_of_segment(const size_t &prm_segment_index ///< The index of the segment whose start arrow should be returned
		                                                                 ) const {
			return ( prm_segment_index > 0                ) ? fragments[ prm_segment_index - 1 ].get_stop_arrow()
			                                                : start_arrow;
		}

		/// \brief Get the stop boundary of the segment with the specified index
		inline const seq_arrow & seq_seg_run::get_stop_arrow_of_segment(const size_t &prm_segment_index ///< The index of the segment whose stop arrow should be returned
		                                                                ) const {
			return ( prm_segment_index < fragments.size() ) ? fragments[ prm_segment_index     ].get_start_arrow()
			                                                : stop_arrow;
		}

		/// \brief Get the length of the specified seq_seg_run's segment corresponding to the specified index
		///
		/// \relates seq_seg_run
		inline residx_t get_length_of_seq_seg(const seq_seg_run &prm_seq_seg_run, ///< The seq_seg_run to query
		                                      const size_t      &prm_seg_idx      ///< The index of the segment who length should be returned
		                                      ) {
			return prm_seq_seg_run.get_stop_arrow_of_segment( prm_seg_idx )
			       - prm_seq_seg_run.get_start_arrow_of_segment( prm_seg_idx );
		}

		/// \brief Get the middle of the specified seq_seg_run's segment corresponding to the specified index
		///
		/// Eg,
		///   * the middle of a seg from 2 to 4 is 3
		///   * the middle of a seg from 2 to 5 is 3.5
		///
		/// \relates seq_seg_run
		inline double get_middle_of_seq_seg(const seq_seg_run &prm_seq_seg_run, ///< The seq_seg_run to query
		                                    const size_t      &prm_seg_idx      ///< The index of the segment who middle index should be returned
		                                    ) {
			// Use res_after() for the start_arrow so that this can handle a seq_seg that starts at 0
			// and then subtract 0.5 at the end to make up the difference
			return (
				static_cast<double>(
					prm_seq_seg_run.get_stop_arrow_of_segment ( prm_seg_idx ).res_after ()
					+
					prm_seq_seg_run.get_start_arrow_of_segment( prm_seg_idx ).res_after()
				)
				/ 2.0
			) - 0.5;
		}

		/// \brief Get the specified seq_seg_run's segment corresponding to the specified index
		///
		/// \relates seq_seg_run
		inline seq_seg get_seq_seg_of_seg_idx(const seq_seg_run &prm_seq_seg_run, ///< The seq_seg_run to query
		                                      const size_t      &prm_seg_idx      ///< The index of the segment to return
		                                      ) {
			return {
				prm_seq_seg_run.get_start_arrow_of_segment( prm_seg_idx ),
				prm_seq_seg_run.get_stop_arrow_of_segment ( prm_seg_idx )
			};
		}

		/// \brief Get a vector of the segments in this seq_seg_run
		///
		/// \relates seq_seg_run
		inline seq_seg_vec get_seq_segs( const seq_seg_run &prm_seq_seg_run ///< The seq_seg_run to query
		                                 ) {
			return common::transform_build<seq_seg_vec>(
				common::indices( prm_seq_seg_run.get_num_segments() ),
				[&] (const size_t &x) {
					return get_seq_seg_of_seg_idx( prm_seq_seg_run, x );
				}
			);
		}

		/// \brief Get the (possibly-repeated, non-sorted) segments from the specified seq_seg_runs
		///
		/// \relates seq_seg_run
		inline seq_seg_vec get_seq_segs( const seq_seg_run_vec &prm_seq_seg_run_vec ///< The seq_seg_runs whose segments should be returned
		                                 ) {
			seq_seg_vec results;
			for (const seq_seg_run &the_seq_seg_run : prm_seq_seg_run_vec) {
				common::append( results, get_seq_segs( the_seq_seg_run ) );
			}
			return results;
		}

		/// \brief Get a vector of the specified seq_seg_runs' segments, sorted by their starts
		inline seq_seg_vec get_start_sorted_seq_segs(const seq_seg_run_vec &prm_seq_seg_run_vec ///< The vector of seq_seg_runs to query
		                                             ) {
			return start_sort_seq_segs_copy( get_seq_segs( prm_seq_seg_run_vec ) );
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
		inline const residx_t & get_start_res_index_of_segment(const seq_seg_run &prm_seq_seg_run,  ///< The seq_seg_run to query
		                                                       const size_t      &prm_segment_index ///< The index of the segment to query
		                                                       ) {
			return prm_seq_seg_run.get_start_arrow_of_segment( prm_segment_index ).res_after();
		}

		/// \brief Get the stop residue index of the segment of specified index in the specified seq_seg_run
		///
		/// \relates seq_seg_run
		inline residx_t get_stop_res_index_of_segment(const seq_seg_run &prm_seq_seg_run,  ///< The seq_seg_run to query
		                                              const size_t      &prm_segment_index ///< The index of the segment to query
		                                              ) {
			return prm_seq_seg_run.get_stop_arrow_of_segment( prm_segment_index ).res_before();
		}

		/// \brief Get the start residue index of the specified seq_seg_run
		///
		/// \relates seq_seg_run
		inline const residx_t & get_start_res_index(const seq_seg_run &prm_seq_seg_run ///< The seq_seg_run to query
		                                            ) {
			return prm_seq_seg_run.get_start_arrow().res_after();
		}

		/// \brief Get the stop residue index of the specified seq_seg_run
		///
		/// \relates seq_seg_run
		inline residx_t get_stop_res_index(const seq_seg_run &prm_seq_seg_run ///< The seq_seg_run to query
		                                   ) {
			return prm_seq_seg_run.get_stop_arrow().res_before();
		}

		/// \brief Get the stop of the first segment in the specified seq_seg_run
		///
		/// \pre `prm_seq_seg_run.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates seq_seg_run
		inline seq_arrow get_stop_of_first_segment(const seq_seg_run &prm_seq_seg_run ///< The seq_seg_run to query
		                                           ) {
			if ( ! prm_seq_seg_run.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_stop_of_first_segment of contiguous seq_seg_run"));
			}
			return prm_seq_seg_run.get_stop_arrow_of_segment( 0 );
		}

		/// \brief Get the start of the last segment in the specified seq_seg_run
		///
		/// \pre `prm_seq_seg_run.is_discontig()` else an invalid_argument_exception will be thrown
		///
		/// \relates seq_seg_run
		inline seq_arrow get_start_of_last_segment(const seq_seg_run &prm_seq_seg_run ///< The seq_seg_run to query
		                                           ) {
			if ( ! prm_seq_seg_run.is_discontig() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_start_of_last_segment of contiguous seq_seg_run"));
			}
			return prm_seq_seg_run.get_start_arrow_of_segment( prm_seq_seg_run.get_num_segments() - 1 );
		}

		/// \brief Get the total length of the specified seq_seg_run (ie the sum of its segments' lengths)
		///
		/// \relates seq_seg_run
		inline residx_t get_total_length(const seq_seg_run &prm_seq_seg_run ///< The seq_seg_run to query
		                                 ) {
			return boost::accumulate(
				common::indices( prm_seq_seg_run.get_num_segments() )
					| boost::adaptors::transformed( [&] (const size_t &x) {
						return get_length_of_seq_seg( prm_seq_seg_run, x );
					} ),
				static_cast<residx_t>( 0 )
			);
		}

		/// \brief Get the middle index of the specified seq_seg_run
		///        (ie the length-weighted average of its segments' middles, as calculated by get_middle_of_seq_seg())
		///
		/// \relates seq_seg_run
		inline double middle_index(const seq_seg_run &prm_seq_seg_run ///< The seq_seg_run to query
		                           ) {
			return boost::accumulate(
				common::indices( prm_seq_seg_run.get_num_segments() )
					| boost::adaptors::transformed( [&] (const size_t &x) {
						return (
							get_length_of_seq_seg( prm_seq_seg_run, x )
							*
							get_middle_of_seq_seg( prm_seq_seg_run, x )
						);
					} ),
				0.0
			) / static_cast<double>( get_total_length( prm_seq_seg_run ) );
		}


		namespace detail {

			/// \brief Implementation function for returning a seq_seg_vec from an even list of starts/stops
			///
			/// \relates seq_seg_run
			template <size_t N, std::size_t... Index>
			seq_seg_vec make_seq_seg_run_from_res_indices_impl(const std::array<residx_t, N> &prm_array,
			                                                   std::index_sequence<Index...>
			                                                   ) {
				return seq_seg_vec{ {
					seq_seg{
						prm_array.at( 2_z * Index     ),
						prm_array.at( 2_z * Index + 1 )
					}...
				} };
			}
		} // namespace detail

		/// \brief Make a seq_seg_run from the even list of residue start/stop indices
		///
		/// \relates seq_seg_run
		template <typename... Ts>
		seq_seg_run make_seq_seg_run_from_res_indices(const Ts &...prm_res_idcs ///< The start/stop residue indices
		                                              ) {
			static_assert( sizeof...( Ts ) >= 2,
				"make_seq_seg_run_from_res_indices() requires at least two index parameters" );

			static_assert( ( sizeof...( Ts ) % 2 ) == 0,
				"make_seq_seg_run_from_res_indices() requires an even number of index parameters" );

			static_assert( ( ... && std::is_integral_v<common::remove_cvref_t<Ts>> ),
				"make_seq_seg_run_from_res_indices() requires all parameters to by of integral type" );

			return seq_seg_run{
				detail::make_seq_seg_run_from_res_indices_impl(
					std::array<residx_t, sizeof...( Ts )>{ { debug_numeric_cast<residx_t>( prm_res_idcs )... } },
					std::make_index_sequence< sizeof...( Ts ) / 2 >{}
				)
			};
		}


		/// \brief Make a continuous seq_seg_run from the residue indices
		///
		/// \relates seq_seg_run
		inline seq_seg_run make_seq_seg_run_from_res_indices(const residx_t &prm_start_res_idx, ///< The start residue index
		                                                     const residx_t &prm_stop_res_idx   ///< The stop residue index
		                                                     ) {
			return {
				arrow_before_res( prm_start_res_idx ),
				arrow_after_res ( prm_stop_res_idx  )
			};
		}

		/// \brief Make a seq_seg_run
		///
		/// \relates seq_seg_run
		inline seq_seg_run make_seq_seg_run_from_res_indices(const residx_residx_pair_vec &prm_residue_index_segments ///< The residue index start/stop pairs of the seq_seg_run's segments
		                                                     ) {
			return seq_seg_run{
				common::transform_build<seq_seg_vec>(
					prm_residue_index_segments,
					seq_seg_of_res_idx_pair
				)
			};
		}

		/// \brief Return whether the either of the two specified seq_seg_runs overlaps, interleaves or straddles the other
		///
		/// \relates seq_seg_run
		inline bool any_interaction(const seq_seg_run &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
		                            const seq_seg_run &prm_seq_seg_run_b  ///< The second seq_seg_run to query
		                            ) {
			return (
				prm_seq_seg_run_a.get_start_arrow() < prm_seq_seg_run_b.get_stop_arrow()
				&&
				prm_seq_seg_run_b.get_start_arrow() < prm_seq_seg_run_a.get_stop_arrow()
			);
		}

		namespace detail {

			/// \brief Get the index of the first segment in the specified seq_seg_run that isn't entirely before
			///        the specified arrow
			inline size_t index_of_first_seg_not_earlier_than_arrow(const seq_seg_run &prm_seq_seg_run, ///< The seq_seg_run to query
			                                                        const seq_arrow   &prm_arrow        ///< The arrow to compare to
			                                                        ) {
				const size_t num_segs       = prm_seq_seg_run.get_num_segments();
				const auto   seg_nums_range = common::indices( num_segs );
				const auto   result_itr     = lower_bound(
					seg_nums_range,
					prm_arrow,
					[&] (const size_t &seg_index, const seq_arrow &other_start_arrow) {
						return prm_seq_seg_run.get_stop_arrow_of_segment( seg_index ) <= other_start_arrow;
					}
				);
				return ( result_itr != ::std::cend( seg_nums_range ) )
					? *result_itr
					: num_segs;
			}

			/// \brief Apply the specified function to all overlapping pairs of seq_segs between the two specified seq_seg_runs
			///
			/// \todo In an ideal world, this would be a range over which functions could iterator, rather than a function
			template <typename Fn>
			inline void apply_to_overlaps_in_seq_seg_runs(const seq_seg_run  &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
			                                              const seq_seg_run  &prm_seq_seg_run_b, ///< The second seq_seg_run to query
			                                              Fn                &&prm_fn             ///< The function to apply
			                                              ) {
				const size_t num_segs_a = prm_seq_seg_run_a.get_num_segments();
				const size_t num_segs_b = prm_seq_seg_run_b.get_num_segments();

				size_t ctr_a = index_of_first_seg_not_earlier_than_arrow( prm_seq_seg_run_a, prm_seq_seg_run_b.get_start_arrow() );
				size_t ctr_b = index_of_first_seg_not_earlier_than_arrow( prm_seq_seg_run_b, prm_seq_seg_run_a.get_start_arrow() );

				while ( ctr_a != num_segs_a && ctr_b != num_segs_b) {
					const seq_seg seg_a = get_seq_seg_of_seg_idx( prm_seq_seg_run_a, ctr_a );
					const seq_seg seg_b = get_seq_seg_of_seg_idx( prm_seq_seg_run_b, ctr_b );
					if ( are_overlapping( seg_a, seg_b ) ) {
						::std::invoke( prm_fn, seg_a, seg_b );
					}
					if ( seg_a.get_stop_arrow() <= seg_b.get_stop_arrow() ) {
						++ctr_a;
					}
					if ( seg_b.get_stop_arrow() <= seg_a.get_stop_arrow() ) {
						++ctr_b;
					}
				}
			}
		} // namespace detail

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
		inline bool are_overlapping(const seq_seg_run &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
		                            const seq_seg_run &prm_seq_seg_run_b  ///< The second seq_seg_run to query
		                            ) {
			if ( ! any_interaction( prm_seq_seg_run_a, prm_seq_seg_run_b ) ) {
				return false;
			}
			const size_t num_segs_a = prm_seq_seg_run_a.get_num_segments();
			const size_t num_segs_b = prm_seq_seg_run_b.get_num_segments();

			// If there are non-trivial numbers of segments, use a linear strategy
			if ( num_segs_a + num_segs_b > 4 ) {
				for (size_t ctr_a = 0, ctr_b = 0 ; ctr_a != num_segs_a && ctr_b != num_segs_b ; ) {
					if ( are_overlapping( get_seq_seg_of_seg_idx( prm_seq_seg_run_a, ctr_a ),
					                      get_seq_seg_of_seg_idx( prm_seq_seg_run_b, ctr_b ) ) ) {
						return true;
					}

					const size_t orig_ctr_a = ctr_a;
					if ( prm_seq_seg_run_a.get_stop_arrow_of_segment( ctr_a ) <= prm_seq_seg_run_b.get_stop_arrow_of_segment(      ctr_b ) ) {
						++ctr_a;
					}
					if ( prm_seq_seg_run_b.get_stop_arrow_of_segment( ctr_b ) <= prm_seq_seg_run_a.get_stop_arrow_of_segment( orig_ctr_a ) ) {
						++ctr_b;
					}
				}
				return false;
			}
			// Otherwise, it turns out to be measurably faster to just do all-vs-all
			else {
				for (const auto &seg_ctr_a : common::indices( num_segs_a ) ) {
					for (const auto &seg_ctr_b : common::indices( num_segs_b ) ) {
						const bool seg_overlap = are_overlapping(
							get_seq_seg_of_seg_idx( prm_seq_seg_run_a, seg_ctr_a ),
							get_seq_seg_of_seg_idx( prm_seq_seg_run_b, seg_ctr_b )
						);
						if ( seg_overlap ) {
							return true;
						}
					}
				}
				return false;
			}
		}

		/// \brief Return the number of residues by which the two specified seq_seg_runs overlap (or 0 if they don't overlap)
		///
		/// \todo Consider seeing whether the lesson of are_overlapping() also applies here
		///       (ie use a dumb strategy for trivial numbers of segments)
		///
		/// \relates seq_seg
		inline residx_t overlap_by(const seq_seg_run &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
		                           const seq_seg_run &prm_seq_seg_run_b  ///< The second seq_seg_run to query
		                           ) {
			if ( ! any_interaction( prm_seq_seg_run_a, prm_seq_seg_run_b ) ) {
				return 0;
			}
			residx_t overlap_size = 0;
			detail::apply_to_overlaps_in_seq_seg_runs(
				prm_seq_seg_run_a,
				prm_seq_seg_run_b,
				[&] (const seq_seg &x, const seq_seg &y) {
					overlap_size += overlap_by( x, y );
				}
			);
			return overlap_size;
		}

		/// \brief Return the length of the shorter of the two specified seq_seg_runs
		///
		/// \relates seq_seg
		inline residx_t shorter_length(const seq_seg_run &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
		                               const seq_seg_run &prm_seq_seg_run_b  ///< The second seq_seg_run to query
		                               ) {
			return std::min(
				get_total_length( prm_seq_seg_run_a ),
				get_total_length( prm_seq_seg_run_b )
			);
		}

		/// \brief Return the length of the longer of the two specified seq_seg_runs
		///
		/// \relates seq_seg
		inline residx_t longer_length(const seq_seg_run &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
		                              const seq_seg_run &prm_seq_seg_run_b  ///< The second seq_seg_run to query
		                              ) {
			return std::max(
				get_total_length( prm_seq_seg_run_a ),
				get_total_length( prm_seq_seg_run_b )
			);
		}

		/// \brief Return the fraction overlap between the two specified seq_seg_runs over the length of the shorter
		///
		/// \relates seq_seg
		inline double fraction_overlap_over_shorter(const seq_seg_run &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
		                                            const seq_seg_run &prm_seq_seg_run_b  ///< The second seq_seg_run to query
		                                            ) {
			return
				static_cast<double>( overlap_by    ( prm_seq_seg_run_a, prm_seq_seg_run_b ) )
				/
				static_cast<double>( shorter_length( prm_seq_seg_run_a, prm_seq_seg_run_b ) );
		}

		/// \brief Return the fraction overlap between the two specified seq_seg_runs over the length of the longer
		///
		/// \relates seq_seg
		inline double fraction_overlap_over_longer(const seq_seg_run &prm_seq_seg_run_a, ///< The first  seq_seg_run to query
		                                           const seq_seg_run &prm_seq_seg_run_b  ///< The second seq_seg_run to query
		                                           ) {
			return
				static_cast<double>( overlap_by   ( prm_seq_seg_run_a, prm_seq_seg_run_b ) )
				/
				static_cast<double>( longer_length( prm_seq_seg_run_a, prm_seq_seg_run_b ) );
		}

		/// \brief Whether the segments in the first specified seq_seg_run never extend outside
		///        those in the second specified seq_seg_run
		///
		/// In other words: the first is within or equal to the second
		///
		/// \relates calc_hit
		inline bool first_is_not_outside_second(const seq_seg_run &prm_seq_seg_run_a, ///< The first  calc_hit to query
		                                        const seq_seg_run &prm_seq_seg_run_b  ///< The second calc_hit to query
		                                        ) {
			if ( prm_seq_seg_run_a.get_start_arrow() < prm_seq_seg_run_b.get_start_arrow() ) {
				return false;
			}
			if ( prm_seq_seg_run_a.get_stop_arrow () > prm_seq_seg_run_b.get_stop_arrow () ) {
				return false;
			}

			const size_t num_segments_lhs = prm_seq_seg_run_a.get_num_segments();
			const size_t num_segments_rhs = prm_seq_seg_run_b.get_num_segments();

			size_t rhs_ctr = 0;
			for (const auto &lhs_ctr : common::indices( num_segments_lhs ) ) {
				while ( rhs_ctr < num_segments_rhs && prm_seq_seg_run_a.get_stop_arrow_of_segment( lhs_ctr ) > prm_seq_seg_run_b.get_stop_arrow_of_segment( rhs_ctr ) ) {
					++rhs_ctr;
				}
				if ( rhs_ctr == num_segments_rhs ) {
					return false;
				}
				if ( prm_seq_seg_run_a.get_start_arrow_of_segment( lhs_ctr ) < prm_seq_seg_run_b.get_start_arrow_of_segment( rhs_ctr ) ) {
					return false;
				}
			}
			return true;
		}

		/// \brief Whether either of the specified seq_seg_run covers the other
		///
		/// \relates calc_hit
		inline bool one_covers_other(const seq_seg_run &prm_seq_seg_run_a, ///< The first  calc_hit to query
		                             const seq_seg_run &prm_seq_seg_run_b  ///< The second calc_hit to query
		                             ) {
			const residx_t length_a = get_total_length( prm_seq_seg_run_a );
			const residx_t length_b = get_total_length( prm_seq_seg_run_b );
			if ( length_a < length_b ) {
				return first_is_not_outside_second( prm_seq_seg_run_a, prm_seq_seg_run_b );
			}
			if ( length_a > length_b ) {
				return first_is_not_outside_second( prm_seq_seg_run_b, prm_seq_seg_run_a );
			}
			return ( prm_seq_seg_run_a == prm_seq_seg_run_b );
		}

		/// \brief Whether the segments in the first specified seq_seg_run are shorter strictly
		///        within those  in the second specified seq_seg_run
		///
		/// In other words: the first is within and not equal to the second
		///
		/// \relates calc_hit
		inline bool first_is_shorter_and_within_second(const seq_seg_run &prm_seq_seg_run_a, ///< The first  calc_hit to query
		                                               const seq_seg_run &prm_seq_seg_run_b  ///< The second calc_hit to query
		                                               ) {
			return (
				get_total_length( prm_seq_seg_run_a ) < get_total_length( prm_seq_seg_run_b )
				&&
				first_is_not_outside_second( prm_seq_seg_run_a, prm_seq_seg_run_b )
			);
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
		inline bool second_right_intersperses_first(const seq_seg_run &prm_seq_seg_run_a, ///< The first  calc_hit to query
		                                            const seq_seg_run &prm_seq_seg_run_b  ///< The second calc_hit to query
		                                            ) {
			return (
				prm_seq_seg_run_a.is_discontig()
				&&
				prm_seq_seg_run_b.is_discontig()
				&&
				prm_seq_seg_run_a.get_start_arrow() < prm_seq_seg_run_b.get_start_arrow()
				&&
				prm_seq_seg_run_a.get_stop_arrow () < prm_seq_seg_run_b.get_stop_arrow ()
				&&
				prm_seq_seg_run_b.get_start_arrow() < prm_seq_seg_run_a.get_stop_arrow ()
				&&
				! are_overlapping( prm_seq_seg_run_a, prm_seq_seg_run_b )
			);
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
		inline bool second_right_or_inside_intersperses_first(const seq_seg_run &prm_seq_seg_run_a, ///< The first  calc_hit to query
		                                                      const seq_seg_run &prm_seq_seg_run_b  ///< The second calc_hit to query
		                                                      ) {
			return (
				prm_seq_seg_run_a.is_discontig()
				&&
				prm_seq_seg_run_b.is_discontig()
				&&
				prm_seq_seg_run_a.get_start_arrow() < prm_seq_seg_run_b.get_start_arrow()
				&&
				prm_seq_seg_run_b.get_start_arrow() < prm_seq_seg_run_a.get_stop_arrow ()
				&&
				! are_overlapping( prm_seq_seg_run_a, prm_seq_seg_run_b )
			);
		}

	} // namespace seq
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_SEG_RUN_HPP
