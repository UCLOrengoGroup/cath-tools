/// \file
/// \brief The full_hit class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FULL_HIT_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FULL_HIT_H

#include <boost/operators.hpp>
#include <boost/optional.hpp>

#include "common/type_aliases.h"
#include "resolve_hits/file/alnd_rgn.h"
#include "resolve_hits/hit_output_format.h"
#include "resolve_hits/hit_score_type.h"
#include "resolve_hits/hit_seg.h"
#include "resolve_hits/resolve_hits_type_aliases.h"
#include "resolve_hits/score_functions.h"
#include "resolve_hits/trim/trim_spec.h"

#include <string>

namespace cath {
	namespace rslv {

		/// \brief Represent a single full_hit (ie one domain) with a score, label and one or more segments
		///
		/// This is the full hit with all information; whereas calc_hit is
		/// a cut down version to be used in the algorithm calculations
		///
		/// Like all cath-resolve_hits code, this assumes simple residue numbering
		/// and is hence unsuitable for use with raw PDB residue numbers.
		class full_hit final : private boost::equality_comparable<full_hit> {
		private:
			/// \brief The list of segments
			hit_seg_vec segments;

			/// \brief The index of the label for this full_hit (in some corresponding list labels)
			std::string label;

			/// \brief The score associated with this full_hit
			///
			/// This must be greater than 0.0
			double the_score;

			/// \brief The type of score that's being used here
			hit_score_type score_type;

			/// \brief Any aligned regions parsed from an hmmsearch output alignment or boost::none otherwise
			alnd_rgn_vec_opt alnd_rgns_opt;

			void sanity_check() const;

		public:
			full_hit(const hit_seg_vec &,
			         const std::string &,
			         const double &,
			         const hit_score_type & = hit_score_type::CRH_SCORE,
			         const alnd_rgn_vec_opt & = boost::none);

			full_hit(const hit_seg_vec &,
			         const std::string &&,
			         const double &,
			         const hit_score_type & = hit_score_type::CRH_SCORE,
			         const alnd_rgn_vec_opt && = boost::none);

			const hit_seg_vec & get_segments() const;
			const std::string & get_label() const;
			const double & get_score() const;
			const hit_score_type & get_score_type() const;
			const alnd_rgn_vec_opt & get_alnd_rgns_opt() const;
		};

		std::string get_score_string(const double &,
		                             const hit_score_type &,
		                             const size_t & = 4);

		std::string get_score_string(const full_hit &,
		                             const size_t & = 4);

		std::string get_segments_string(const full_hit &,
		                                const boost::optional<trim_spec> & = boost::none);
		std::string to_string(const full_hit &,
		                      const hit_output_format & = hit_output_format::CLASS,
		                      const std::string & = std::string{},
		                      const boost::optional<trim_spec> & = boost::none);

		/// \brief Sanity check that the full_hit is sensible and throw an exception if not
		inline void full_hit::sanity_check() const {
			if ( ! segments_are_start_sorted_and_non_overlapping( segments ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Hit's segments must be start-sorted and non-overlapping"));
			}
			if ( the_score <= 0 ) {
				if ( score_type != hit_score_type::FULL_EVALUE || the_score < 0 ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Hit's score must be strictly greater than 0 else the algorithm doesn't work (because there's no way to no how to trade scores off against empty space)"));
				}
			}
		}

		/// \brief Ctor
		inline full_hit::full_hit(const hit_seg_vec      &arg_segments,     ///< The segments of the full_hit
		                          const std::string      &arg_label,        ///< The label of the hits' match protein
		                          const double           &arg_score,        ///< The score associated with the full_hit
		                          const hit_score_type   &arg_score_type,   ///< The type of score stored in this hit (eg evalue / bitscore / crh-score)
		                          const alnd_rgn_vec_opt &arg_alnd_rgns_opt ///< Any aligned regions parsed from an hmmsearch output alignment or boost::none otherwise
		                          ) : segments     { arg_segments      },
		                              label        { arg_label         },
		                              the_score    { arg_score         },
		                              score_type   { arg_score_type    },
		                              alnd_rgns_opt{ arg_alnd_rgns_opt } {
			sanity_check();
		}

		/// \brief Ctor
		inline full_hit::full_hit(const hit_seg_vec       &arg_segments,     ///< The segments of the full_hit
		                          const std::string      &&arg_label,        ///< The label of the hits' match protein
		                          const double            &arg_score,        ///< The score associated with the full_hit
		                          const hit_score_type    &arg_score_type,   ///< The type of score stored in this hit (eg evalue / bitscore / crh-score)
		                          const alnd_rgn_vec_opt &&arg_alnd_rgns_opt ///< Any aligned regions parsed from an hmmsearch output alignment or boost::none otherwise
		                          ) : segments     { arg_segments           },
		                              label        { std::move( arg_label ) },
		                              the_score    { arg_score              },
		                              score_type   { arg_score_type         },
		                              alnd_rgns_opt{ arg_alnd_rgns_opt      } {
			sanity_check();
		}

		/// \brief Getter for the segments of the full_hit
		inline const hit_seg_vec & full_hit::get_segments() const {
			return segments;
		}
		
		/// \brief Getter for the label of the hits' match protein
		inline const std::string & full_hit::get_label() const {
			return label;
		}

		/// \brief Getter for the score associated with the full_hit
		inline const double & full_hit::get_score() const {
			return the_score;
		}

		/// \brief Getter for the type of score stored in this hit (eg evalue / bitscore / crh-score)
		inline const hit_score_type & full_hit::get_score_type() const {
			return score_type;
		}

		/// \brief Getter for aligned regions if they were parsed from a hmmsearch files
		inline const alnd_rgn_vec_opt & full_hit::get_alnd_rgns_opt() const {
			return alnd_rgns_opt;
		}

		/// \brief Return whether the two specified full_hits are identical
		///
		/// \relates full_hit
		inline bool operator==(const full_hit &arg_lhs, ///< The first  full_hit to compare
		                       const full_hit &arg_rhs  ///< The second full_hit to compare
		                       ) {
			return (
				( arg_lhs.get_segments()                  == arg_rhs.get_segments()                  )
				&&
				( arg_lhs.get_label()                     == arg_rhs.get_label()                     )
				&&
				( arg_lhs.get_score()                     == arg_rhs.get_score()                     )
				&&
				( arg_lhs.get_score_type()                == arg_rhs.get_score_type()                )
			);
		}

		/// \brief Return whether this full_hit is discontiguous
		inline bool is_discontig(const full_hit &arg_full_hit
		                         ) {
			return ( arg_full_hit.get_segments().size() > 1 );
		}

		/// \brief Get the specified full_hit's segment corresponding to the specified index
		///
		/// \relates full_hit
		inline const hit_seg & get_hit_seg_of_seg_idx(const full_hit &arg_full_hit, ///< The full_hit to query
		                                              const size_t   &arg_seg_idx   ///< The index of the segment to return
		                                              ) {
			return arg_full_hit.get_segments()[ arg_seg_idx ];
		}

		/// \brief Get the length of the specified full_hit's segment corresponding to the specified index
		///
		/// \relates full_hit
		inline size_t get_length_of_hit_seg(const full_hit &arg_full_hit, ///< The full_hit to query
		                                    const size_t   &arg_seg_idx   ///< The index of the segment who length should be returned
		                                    ) {
			return get_length( get_hit_seg_of_seg_idx( arg_full_hit, arg_seg_idx ) );
		}

		/// \brief Get the start residue index of the segment of specified index in the specified full_hit
		///
		/// \relates full_hit
		inline const residx_t & get_start_res_index_of_segment(const full_hit &arg_full_hit,     ///< The full_hit to query
		                                                       const size_t   &arg_segment_index ///< The index of the segment to query
		                                                       ) {
			return get_start_res_index( get_hit_seg_of_seg_idx( arg_full_hit, arg_segment_index ) );
		}

		/// \brief Get the stop residue index of the segment of specified index in the specified full_hit
		///
		/// \relates full_hit
		inline residx_t get_stop_res_index_of_segment(const full_hit &arg_full_hit,     ///< The full_hit to query
		                                              const size_t   &arg_segment_index ///< The index of the segment to query
		                                              ) {
			return get_stop_res_index( get_hit_seg_of_seg_idx( arg_full_hit, arg_segment_index ) );
		}

		/// \brief Get the start res_arrow of the specified full_hit
		///
		/// \relates full_hit
		inline const res_arrow & get_start_res_arrow(const full_hit &arg_full_hit ///< The full_hit to query
		                                             ) {
			return arg_full_hit.get_segments().front().get_start_arrow();
		}

		/// \brief Get the stop res_arrow of the specified full_hit
		///
		/// \relates full_hit
		inline const res_arrow & get_stop_res_arrow(const full_hit &arg_full_hit ///< The full_hit to query
		                                            ) {
			return arg_full_hit.get_segments().back().get_stop_arrow();
		}

		/// \brief Get the start residue index of the specified full_hit
		///
		/// \relates full_hit
		inline const residx_t & get_start_res_index(const full_hit &arg_full_hit ///< The full_hit to query
		                                            ) {
			return get_start_res_index( arg_full_hit.get_segments().front() );
		}

		/// \brief Get the stop residue index of the specified full_hit
		///
		/// \relates full_hit
		inline residx_t get_stop_res_index(const full_hit &arg_full_hit ///< The full_hit to query
		                                   ) {
			return get_stop_res_index( arg_full_hit.get_segments().back() );
		}

		/// \brief Get the stop of the first segment in the specified full_hit
		///
		/// \pre `is_discontig( arg_full_hit )` else an invalid_argument_exception will be thrown
		///
		/// \relates full_hit
		inline const res_arrow & get_stop_of_first_segment(const full_hit &arg_full_hit ///< The full_hit to query
		                                                   ) {
			if ( ! is_discontig( arg_full_hit ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_stop_of_first_segment of contiguous full_hit"));
			}
		
			return arg_full_hit.get_segments().front().get_stop_arrow();
		}

		/// \brief Get the start of the last segment in the specified full_hit
		///
		/// \pre `is_discontig( arg_full_hit )` else an invalid_argument_exception will be thrown
		///
		/// \relates full_hit
		inline res_arrow get_start_of_last_segment(const full_hit &arg_full_hit ///< The full_hit to query
		                                           ) {
			if ( ! is_discontig( arg_full_hit ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get_start_of_last_segment of contiguous full_hit"));
			}
			return arg_full_hit.get_segments().back().get_start_arrow();
		}

		/// \brief Get the total length of the specified full_hit (ie the sum of its segments' lengths)
		///
		/// \relates full_hit
		inline size_t get_total_length(const full_hit &arg_full_hit ///< The full_hit to query
		                               ) {
			return boost::accumulate(
				arg_full_hit.get_segments()
					| boost::adaptors::transformed( [&] (const hit_seg &x) {
						return get_length( x );
					} ),
				0_z
			);
		}

		/// \brief Get the crh-score associated with the specified hit, calculated according to the specified crh_score_spec
		inline resscr_t get_crh_score(const full_hit       &arg_full_hit,  ///< The full_hit to query
		                              const crh_score_spec &arg_score_spec ///< The crh_score_spec to use in any score calculations that are required
		                              ) {
			switch ( arg_full_hit.get_score_type() ) {
				case( hit_score_type::FULL_EVALUE ) : {
					return crh_score_of_evalue(
						arg_full_hit.get_score(),
						get_total_length( arg_full_hit ),
						arg_score_spec
					);
				}
				case( hit_score_type::BITSCORE    ) : {
					return crh_score_of_pseudo_bitscore(
						arg_full_hit.get_score(),
						get_total_length( arg_full_hit ),
						arg_score_spec
					);
				}
				case( hit_score_type::CRH_SCORE   ) : {
					return debug_numeric_cast<resscr_t>( arg_full_hit.get_score() );
				}
				default : {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of hit_score_type not recognised whilst getting crh_score from full_hit"));
					return 0; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
				}
			}
		}

	} // namespace rslv
} // namespace cath

#endif
