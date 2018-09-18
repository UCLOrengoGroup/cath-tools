/// \file
/// \brief The calc_hit_list class definitions

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

#include "calc_hit_list.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/equal.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/reverse.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/upper_bound.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/spirit/include/qi.hpp>

#include "common/algorithm/remove_itrs_from_range.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/make_string_ref.hpp"
#include "common/boost_addenda/range/adaptor/adjacented.hpp"
#include "common/boost_addenda/range/back.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/file/open_fstream.hpp"
#include "common/type_aliases.hpp"
#include "resolve_hits/detail/calc_hit_prune_builder.hpp"
#include "resolve_hits/first_hit_is_better.hpp"
#include "resolve_hits/full_hit_fns.hpp"
#include "resolve_hits/full_hit_list.hpp"
#include "resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"

#include <chrono>
#include <fstream>
#include <string>

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace cath::seq;

using boost::adaptors::filtered;
using boost::empty;
using boost::filesystem::path;
using boost::integer_range;
using boost::irange;
using boost::make_optional;
using boost::none;
using boost::numeric_cast;
using boost::range::lower_bound;
using boost::range::max_element;
using boost::range::reverse;
using boost::range::sort;
using boost::range::upper_bound;
using boost::spirit::qi::double_;
using boost::spirit::qi::omit;
using boost::spirit::qi::parse;
using boost::spirit::uint_;
using boost::string_ref;
using std::distance;
using std::find_if;
using std::ifstream;
using std::istream;
using std::ostream;
using std::pair;
using std::string;
using std::vector;

/// \brief Make a calc_hit_list from a full_hit_list
///
/// \relates calc_hit_list
calc_hit_vec cath::rslv::make_hit_list_from_full_hit_list(const full_hit_list    &prm_full_hit_list, ///< The full_hit_list to convert
                                                          const crh_score_spec   &prm_score_spec,    ///< The crh_score_spec to specify how the crh-scores are to be calculated from the full-hits
                                                          const crh_segment_spec &prm_segment_spec,  ///< The crh_segment_spec to specify how the segments are to be handled before being put into the hits for calculation
                                                          const crh_filter_spec  &prm_filter_spec    ///< The crh_filter_spec specifying how hits should be filtered
                                                          ) {
	const trim_spec &overlap_trim_spec  = prm_segment_spec.get_overlap_trim_spec();
	const residx_t  &min_seg_length     = prm_segment_spec.get_min_seg_length();
	const auto       seg_long_enough_fn = [&] (const seq_seg &x) { return ( get_length( x ) >= min_seg_length ); };
	bool hit_failed_seg_length = false;

	return transform_build<calc_hit_vec>(
		indices( prm_full_hit_list.size() )
			| filtered( [&] (const size_t &x) {
				const full_hit &full_hit_x = prm_full_hit_list[ x ];
				if ( empty( full_hit_x.get_segments() | filtered( seg_long_enough_fn ) ) ) {
					if ( ! hit_failed_seg_length ) {
						BOOST_LOG_TRIVIAL( warning )
							<< "At least one hit (with match ID "
							<< full_hit_x.get_label()
							<< ") in list has no segments that meet the min-seg-length "
							<< ::std::to_string( min_seg_length );
						hit_failed_seg_length = true;
					}
					return false;
				}
				return score_passes_filter( prm_filter_spec, full_hit_x.get_score(), full_hit_x.get_score_type() );
			} ),
		[&] (const size_t &x) {
			const full_hit &the_hit       = prm_full_hit_list[ x ];
			const auto      filtered_segs = the_hit.get_segments() | filtered( seg_long_enough_fn );
			return calc_hit{
				trim_seq_seg_copy( front( filtered_segs ), overlap_trim_spec ).get_start_arrow(),
				trim_seq_seg_copy( back ( filtered_segs ), overlap_trim_spec ).get_stop_arrow(),
				transform_build<seq_seg_vec>(
					filtered_segs | adjacented,
					[&] (const pair<seq_seg, seq_seg> &y) {
						return seq_seg{
							trim_seq_seg_copy( y.first,  overlap_trim_spec ).get_stop_arrow(),
							trim_seq_seg_copy( y.second, overlap_trim_spec ).get_start_arrow()
						};
					}
				),
				get_crh_score( the_hit, prm_score_spec ),
				debug_numeric_cast<hitidx_t>( x )
			};
		}
	);
}

/// \brief Make a vector of calc_hits from the specified full_hit_list that's both sorted
///        and pruned of any hits that are redundant (because there's at least one other
///        hit in the list that's strictly better than it)
///
///
/// \relates calc_hit_list
calc_hit_vec cath::rslv::make_sorted_pruned_calc_hit_vec(const full_hit_list       &prm_full_hit_list, ///< The full_hit_list to convert
                                                         const crh_score_spec      &prm_score_spec,    ///< The crh_score_spec to specify how the crh-scores are to be calculated from the full-hits
                                                         const crh_segment_spec    &prm_segment_spec,  ///< The crh_segment_spec to specify how the segments are to be handled before being put into the hits for calculation
                                                         const crh_filter_spec     &prm_filter_spec,   ///< The crh_filter_spec specifying how hits should be filtered
                                                         const seg_dupl_hit_policy &prm_policy         ///< Whether the strictly-worse hits should be pruned
                                                         ) {
	const trim_spec &overlap_trim_spec = prm_segment_spec.get_overlap_trim_spec();
	const residx_t  &min_seg_length    = prm_segment_spec.get_min_seg_length();

	bool failed_seg_length = false;
	detail::calc_hit_prune_builder the_builder{ prm_policy };
	the_builder.reserve( prm_full_hit_list.size() );

	seq_seg_vec trimmed_segs;

	for (const size_t &full_hit_ctr : indices( prm_full_hit_list.size() ) ) {
		const full_hit &the_full_hit = prm_full_hit_list[ full_hit_ctr ];

		if ( ! score_passes_filter( prm_filter_spec, the_full_hit.get_score(), the_full_hit.get_score_type() ) ) {
			continue;
		}

		const auto full_segs          = the_full_hit.get_segments();
		const auto seg_long_enough_fn = [&] (const seq_seg &x) { return ( get_length( x ) >= min_seg_length ); };
		const auto filtered_segs      = full_segs | filtered( seg_long_enough_fn );

		trimmed_segs.clear();
		for (const auto &the_seg : filtered_segs) {
			trimmed_segs.push_back( trim_seq_seg_copy( the_seg, overlap_trim_spec ) );
		}

		if ( empty( trimmed_segs ) ) {
			if ( ! failed_seg_length ) {
				BOOST_LOG_TRIVIAL( warning )
					<< "At least one hit (with match ID "
					<< the_full_hit.get_label()
					<< ") in list has no segments that meet the min-seg-length "
					<< ::std::to_string( min_seg_length );
				failed_seg_length = true;
			}
			continue;
		}

		seq_seg_vec fragments;
		if ( trimmed_segs.size() > 1 ) {
			fragments.reserve( trimmed_segs.size() - 1 );
		}
		for (const size_t seg_frag_ctr : boost::irange( 1_z, trimmed_segs.size() ) ) {
			fragments.emplace_back(
				trimmed_segs[ seg_frag_ctr - 1 ].get_stop_arrow(),
				trimmed_segs[ seg_frag_ctr     ].get_start_arrow()
			);
		}

		the_builder.add_hit(
			calc_hit{
				front( trimmed_segs ).get_start_arrow(),
				back ( trimmed_segs ).get_stop_arrow(),
				std::move( fragments ),
				get_crh_score( the_full_hit, prm_score_spec ),
				debug_numeric_cast<hitidx_t>( full_hit_ctr )
			},
			prm_full_hit_list
		);
	}

	calc_hit_vec result_hits{ std::move( the_builder.get_built_hits() ) };
	boost::range::sort(
		result_hits,
		calc_hit_list::get_less_than_fn( prm_full_hit_list )
	);

	return result_hits;
}

/// \brief Read a calc_hit_list from the specified file
///
/// \relates calc_hit_list
void cath::rslv::read_hit_list_from_file(read_and_process_mgr &prm_read_and_process_mgr, ///< The read_and_process_mgr to which each calc_hit should be sent when it's read
                                         const path           &prm_file,                 ///< The file from which to read the hits data
                                         const hit_score_type &prm_score_type            ///< The type of score
                                         ) {
	ifstream the_ifstream;
	open_ifstream( the_ifstream, prm_file );

	read_hit_list_from_istream( prm_read_and_process_mgr, the_ifstream, prm_score_type );

	the_ifstream.close();
}

/// \brief Read a calc_hit_list from the specified istream
///
/// \todo Alter this to use seq_seg_run_parser
///
/// \relates calc_hit_list
void cath::rslv::read_hit_list_from_istream(read_and_process_mgr &prm_read_and_process_mgr, ///< The read_and_process_mgr to which each calc_hit should be sent when it's read
                                            istream              &prm_istream,              ///< The istream from which to read the hits data
                                            const hit_score_type &prm_score_type            ///< The type of score
                                            ) {
	if ( prm_score_type != hit_score_type::FULL_EVALUE && prm_score_type != hit_score_type::CRH_SCORE ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}

	prm_read_and_process_mgr.process_all_outstanding();

	const auto is_space_char     = [ ] (const auto &x) { return ( ( x == ' ' ) || ( x == '\t' ) ); };
	const auto is_non_space_char = [ ] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); };
	const auto find_space        = [&] (const auto &b, const auto &e) {
		return find_if( b, e, is_space_char     );
	};
	const auto find_non_space    = [&] (const auto &b, const auto &e ) {
		return find_if( b, e, is_non_space_char );
	};

	double     score;
	residx_vec bounds;
	string     line;
	string     query_id;

	const auto bounds_pusher = [&] (const residx_t &x) { bounds.push_back( x ); };

	// Store the query IDs seen so far if the crh_filter_spec specifies a limit on the number of queries
	query_id_recorder seen_query_ids;

	while ( getline( prm_istream, line ) ) {
		const auto line_begin_itr        = common::cbegin ( line );
		const auto line_end_itr          = common::cend   ( line );

		const auto end_of_query_id_itr   = find_space     ( line_begin_itr,        line_end_itr        );
		const auto query_id_str_ref      = make_string_ref( line_begin_itr,        end_of_query_id_itr );

		// If this query ID should be skipped, then skip this entry.
		// The function also updates seen_query_ids if not skipping this query ID
		if ( should_skip_query_and_update( prm_read_and_process_mgr, query_id_str_ref, seen_query_ids ) ) {
			continue;
		}

		const auto begin_of_match_id_itr = find_non_space ( end_of_query_id_itr,   line_end_itr        );
		const auto end_of_match_id_itr   = find_space     ( begin_of_match_id_itr, line_end_itr        );
		const auto begin_of_score_itr    = find_non_space ( end_of_match_id_itr,   line_end_itr        );

		// If the line contains nothing but whitespace, skip it
		if ( end_of_query_id_itr == line_begin_itr && begin_of_match_id_itr == line_end_itr ) {
			continue;
		}

		bounds.clear();
		auto parse_itr = begin_of_score_itr;
		const bool ok = parse(
			parse_itr,
			line_end_itr,
			   double_
			>> omit[ +boost::spirit::qi::space ]
			>> uint_[ bounds_pusher ]
			>> omit[ '-' ] //
			>> uint_[ bounds_pusher ]
			>> *(
				','
				>> uint_[ bounds_pusher ]
				>> omit[ "-" ] //
				>> uint_[ bounds_pusher ]
			),
			score
		);

		if ( ! ok || parse_itr != line_end_itr ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception( "Error on attempt to parse line : " + line ));
		}
		if ( bounds.empty() ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception( "No bounds" ));
		}
		if ( bounds.size() % 2 != 0 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception( "Odd number of bounds" ));
		}

		prm_read_and_process_mgr.add_hit(
			query_id_str_ref,
			segments_from_bounds( bounds ),
			string{ begin_of_match_id_itr, end_of_match_id_itr },
			score,
			prm_score_type
		);
	}

	prm_read_and_process_mgr.process_all_outstanding();
}

/// \brief Generate a string describing the specified calc_hit_list
///
/// \relates calc_hit_list
string cath::rslv::to_string(const calc_hit_list &prm_calc_hit_list ///< The calc_hit_list to describe
                             ) {
	return "calc_hit_list["
		+ ::std::to_string( prm_calc_hit_list.size() )
		+ " hits]";
}

/// \brief Insert a description of the specified calc_hit_list into the specified ostream
///
/// \relates calc_hit_list
ostream & cath::rslv::operator<<(ostream             &prm_ostream,      ///< The ostream into which the description should be inserted
                                 const calc_hit_list &prm_calc_hit_list ///< The calc_hit_list to describe
                                 ) {
	prm_ostream << to_string( prm_calc_hit_list );
	return prm_ostream;
}

/// \brief Get the maximum stop residue of all the hits in the specified calc_hit_list
///        or none if the calc_hit_list is empty
///
/// \relates calc_hit_list
residx_opt cath::rslv::get_max_stop(const calc_hit_list &prm_calc_hit_list ///< The calc_hit_list to query
                                    ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if prm_full_hit_list's empty
	return prm_calc_hit_list.empty()
		? none
		: make_optional( get_stop_res_index( *max_element(
			prm_calc_hit_list,
			[] (const calc_hit &x, const calc_hit &y) {
				return get_stop_arrow( x ) < get_stop_arrow( y );
		} ) ) );
}


/// \brief Get the best (ie maximum) score of all the hits in the specified calc_hit_list
///        or none if the calc_hit_list is empty
///
/// \relates calc_hit_list
resscr_opt cath::rslv::get_best_score(const calc_hit_list &prm_calc_hit_list ///< The calc_hit_list to query
                                      ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if prm_full_hit_list's empty
	return prm_calc_hit_list.empty()
		? none
		: make_optional( max_element(
			prm_calc_hit_list,
			[] (const calc_hit &x, const calc_hit &y) {
				return x.get_score() < y.get_score();
			}
		)->get_score() );
}

/// \brief Find the first calc_hit in the specified calc_hit_list that stops at or after the specified residue boundary
///
/// Note: may return one-after-end iterator
///
/// \relates calc_hit_list
calc_hit_list::const_iterator cath::rslv::find_first_hit_stopping_at_or_after(const calc_hit_list &prm_calc_hit_list, ///< The calc_hit_list to query
                                                                              const seq_arrow     &prm_res_arrow      ///< The boundary that the calc_hit must stop at or after
                                                                              ) {
	return lower_bound(
		prm_calc_hit_list,
		prm_res_arrow,
		[] (const calc_hit &h, const seq_arrow &a) {
			return ( get_stop_arrow( h ) < a );
		}
	);
}

/// \brief Find the first calc_hit in the specified calc_hit_list that stops after the specified residue boundary
///
/// Note: may return one-after-end iterator
///
/// \relates calc_hit_list
calc_hit_list::const_iterator cath::rslv::find_first_hit_stopping_after(const calc_hit_list &prm_calc_hit_list, ///< The calc_hit_list to query
                                                                        const seq_arrow     &prm_res_arrow      ///< The boundary that the calc_hit must stop after
                                                                        ) {
	return upper_bound(
		prm_calc_hit_list,
		prm_res_arrow,
		[] (const seq_arrow &a, const calc_hit &h) {
			return ( a < get_stop_arrow( h ) );
		}
	);
}

/// \brief Get a range of indices corresponding to those hits in the specified calc_hit_list that stop in
///        the specified range
///
/// The range is inclusive at the end (ie the returned range may include hits that stop at prm_stop_arrow)
///
/// \relates calc_hit_list
integer_range<hitidx_t> cath::rslv::indices_of_hits_that_stop_in_range(const calc_hit_list &prm_calc_hit_list, ///< The calc_hit_list to query
                                                                       const seq_arrow     &prm_start_arrow,   ///< The start arrow that the hits corresponding to the returned indices must stop after
                                                                       const seq_arrow     &prm_stop_arrow     ///< The stop  arrow that the hits corresponding to the returned indices must stop before or at
                                                                       ) {
	return {
		numeric_cast<hitidx_t>( distance( cbegin( prm_calc_hit_list ), find_first_hit_stopping_after( prm_calc_hit_list, prm_start_arrow ) ) ),
		numeric_cast<hitidx_t>( distance( cbegin( prm_calc_hit_list ), find_first_hit_stopping_after( prm_calc_hit_list, prm_stop_arrow  ) ) )
	};
}

/// \brief Generate a list of iterators to the calc_hits in the specified list that are redundant
///        (because there's at least one other hit in the list that's strictly better than it)
///
/// \relates calc_hit_list
calc_hit_vec_citr_vec cath::rslv::identify_redundant_hits(const calc_hit_vec  &prm_calc_hits, ///< The calc_hit_list to query
                                                          const full_hit_list &prm_full_hits  ///< The corresponding full_hit_list that is used for names that can be used to pick a "better" of two otherwise identical hits
                                                          ) {
	const auto rend_itr = common::rend( prm_calc_hits );
	using chl_citr      = calc_hit_list::const_iterator;
	using chl_citr_vec  = vector<chl_citr>;

	chl_citr_vec active_hit_itrs;
	chl_citr_vec to_be_removed_itrs;
	to_be_removed_itrs.reserve( prm_calc_hits.size() );

	for (auto the_ritr = common::rbegin( prm_calc_hits ); the_ritr != rend_itr; ++the_ritr) {
		const auto &this_hit = *the_ritr;

		const auto active_hitrs_end = common::cend( active_hit_itrs );
		auto active_hit_write_itr   = std::begin  ( active_hit_itrs );

		for (const auto &active_hit_itr : active_hit_itrs ) {
			const auto &active_hit = *active_hit_itr;

			if ( get_start_arrow( active_hit ) < get_stop_arrow( this_hit ) ) {
				const auto result = first_hit_is_better( this_hit, active_hit, prm_full_hits );
				if ( is_true ( result ) ) {
					to_be_removed_itrs.push_back( active_hit_itr );
				}
				else {
					if ( &*active_hit_write_itr != &active_hit_itr ) {
						*active_hit_write_itr = active_hit_itr;
					}
					++active_hit_write_itr;

					if ( is_false ( result ) ) {
						to_be_removed_itrs.push_back( next( the_ritr ).base() );
						break;
					}
				}
			}
		}

		// If any actives have been removed then erase the spares at the end
		if ( active_hit_write_itr != active_hitrs_end ) {
			active_hit_itrs.erase( active_hit_write_itr, active_hitrs_end );
		}

		// If this element hasn't been marked for removal, then put it into the actives for future comparisons
		if ( to_be_removed_itrs.empty() || to_be_removed_itrs.back() != next( the_ritr ).base() ) {
			active_hit_itrs.push_back( next( the_ritr ).base() );
		}
	}

	reverse( to_be_removed_itrs );
	sort( to_be_removed_itrs );

	return to_be_removed_itrs;
}

/// \brief Remove all redundant hits of the specified presorted calc_hit_vec
///
/// \relates calc_hit_list
void cath::rslv::remove_redundant_hits(calc_hit_vec        &prm_calc_hits, ///< The calc_hit_list to query, which must be presorted under calc_hit_list::get_less_than_fn
                                       const full_hit_list &prm_full_hits  ///< The full_hits to use for calc_hit comparisons
                                       ) {
	prm_calc_hits.erase(
		remove_itrs_from_range(
			prm_calc_hits,
			identify_redundant_hits(
				prm_calc_hits,
				prm_full_hits
			)
		),
		common::cend( prm_calc_hits )
	);
}
