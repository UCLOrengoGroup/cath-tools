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
#include <boost/range/algorithm/upper_bound.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/spirit/include/qi.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/make_string_ref.hpp"
#include "common/boost_addenda/range/adaptor/adjacented.hpp"
#include "common/boost_addenda/range/back.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/file/open_fstream.hpp"
#include "common/type_aliases.hpp"
#include "resolve_hits/full_hit_list.hpp"
#include "resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"

#include <chrono>
#include <fstream>
#include <string>

using namespace cath::common;
using namespace cath::rslv;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::is_any_of;
using boost::algorithm::is_space;
using boost::contains;
using boost::empty;
using boost::filesystem::path;
using boost::integer_range;
using boost::irange;
using boost::make_optional;
using boost::none;
using boost::numeric_cast;
using boost::range::equal;
using boost::range::lower_bound;
using boost::range::max_element;
using boost::range::upper_bound;
using boost::spirit::qi::double_;
using boost::spirit::qi::omit;
using boost::spirit::qi::parse;
using boost::spirit::uint_;
using boost::string_ref;
using boost::sub_range;
using boost::token_compress_on;
using std::distance;
using std::find_if;
using std::ifstream;
using std::istream;
using std::make_pair;
using std::ostream;
using std::pair;
using std::string;
using std::vector;

/// \brief Make a calc_hit_list from a full_hit_list
///
/// \relates calc_hit_list
calc_hit_vec cath::rslv::make_hit_list_from_full_hit_list(const full_hit_list    &arg_full_hit_list, ///< The full_hit_list to convert
                                                          const crh_score_spec   &arg_score_spec,    ///< The crh_score_spec to specify how the crh-scores are to be calculated from the full-hits
                                                          const crh_segment_spec &arg_segment_spec,  ///< The crh_segment_spec to specify how the segments are to be handled before being put into the hits for calculation
                                                          const crh_filter_spec  &arg_filter_spec    ///< The crh_filter_spec specifying how hits should be filtered
                                                          ) {
	const trim_spec &overlap_trim_spec  = arg_segment_spec.get_overlap_trim_spec();
	const residx_t  &min_seg_length     = arg_segment_spec.get_min_seg_length();
	const auto       seg_long_enough_fn = [&] (const hit_seg &x) { return ( get_length( x ) >= min_seg_length ); };
	bool hit_failed_seg_length = false;

	return transform_build<calc_hit_vec>(
		irange( 0_z, arg_full_hit_list.size() )
			| filtered( [&] (const size_t &x) {
				const full_hit &full_hit_x = arg_full_hit_list[ x ];
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
				return score_passes_filter( arg_filter_spec, full_hit_x.get_score(), full_hit_x.get_score_type() );
			} ),
		[&] (const size_t &x) {
			const full_hit &the_hit       = arg_full_hit_list[ x ];
			const auto      filtered_segs = the_hit.get_segments() | filtered( seg_long_enough_fn );
			return calc_hit{
				trim_hit_seg_copy( front( filtered_segs ), overlap_trim_spec ).get_start_arrow(),
				trim_hit_seg_copy( back ( filtered_segs ), overlap_trim_spec ).get_stop_arrow(),
				transform_build<hit_seg_vec>(
					filtered_segs | adjacented,
					[&] (const pair<hit_seg, hit_seg> &x) {
						return hit_seg{
							trim_hit_seg_copy( x.first,  overlap_trim_spec ).get_stop_arrow(),
							trim_hit_seg_copy( x.second, overlap_trim_spec ).get_start_arrow()
						};
					}
				),
				get_crh_score( the_hit, arg_score_spec ),
				debug_numeric_cast<hitidx_t>( x )
			};
		}
	);
}

/// \brief Read a calc_hit_list from the specified file
///
/// \relates calc_hit_list
void cath::rslv::read_hit_list_from_file(read_and_process_mgr &arg_read_and_process_mgr, ///< The read_and_process_mgr to which each calc_hit should be sent when it's read
                                         const path           &arg_file,                 ///< The file from which to read the hits data
                                         const hit_score_type &arg_score_type            ///< The type of score
                                         ) {
	ifstream the_ifstream;
	open_ifstream( the_ifstream, arg_file );

	read_hit_list_from_istream( arg_read_and_process_mgr, the_ifstream, arg_score_type );

	the_ifstream.close();
}

/// \brief Read a calc_hit_list from the specified istream
///
/// \relates calc_hit_list
void cath::rslv::read_hit_list_from_istream(read_and_process_mgr &arg_read_and_process_mgr, ///< The read_and_process_mgr to which each calc_hit should be sent when it's read
                                            istream              &arg_istream,              ///< The istream from which to read the hits data
                                            const hit_score_type &arg_score_type            ///< The type of score
                                            ) {
	if ( arg_score_type != hit_score_type::FULL_EVALUE && arg_score_type != hit_score_type::CRH_SCORE ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
	}

	arg_read_and_process_mgr.process_all_outstanding();

	const auto is_space_char     = [ ] (const auto &x) { return ( ( x == ' ' ) || ( x == '\t' ) ); };
	const auto is_non_space_char = [ ] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); };
	const auto find_space        = [&] (const auto &b, const auto &e) {
		return find_if( b, e, is_space_char     );
	};
	const auto find_non_space    = [&] (const auto &b, const auto &e ) {
		return find_if( b, e, is_non_space_char );
	};

	double           score;
	hit_seg_vec      segments;
	vector<residx_t> bounds;
	string           line;
	string           query_id;

	const auto bounds_pusher = [&] (const residx_t &x) { bounds.push_back( x ); };

	while ( getline( arg_istream, line ) ) {
		const auto line_begin_itr        = common::cbegin ( line );
		const auto line_end_itr          = common::cend   ( line );

		const auto end_of_query_id_itr   = find_space     ( line_begin_itr,        line_end_itr        );
		const auto query_id_str_ref      = make_string_ref( line_begin_itr,        end_of_query_id_itr );

		// If there are filter query IDs but this isn't amongst them, then skip this entry
		if ( should_skip_query_id( arg_read_and_process_mgr, query_id_str_ref ) ) {
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

		// *** SHOULD CHECK THAT BOUNDS IS STRICTLY ASCENDING ****
		segments.clear();
		for (const size_t &bound_ctr : irange( 0_z, bounds.size(), 2_z ) ) {
			segments.emplace_back(
				arrow_before_res( bounds[ bound_ctr     ] ),
				arrow_after_res ( bounds[ bound_ctr + 1 ] )
			);
		}

		arg_read_and_process_mgr.add_hit(
			query_id_str_ref,
			segments,
			string{ begin_of_match_id_itr, end_of_match_id_itr },
			score,
			arg_score_type
		);
	}

	arg_read_and_process_mgr.process_all_outstanding();
}

/// \brief Generate a string describing the specified calc_hit_list
///
/// \relates calc_hit_list
string cath::rslv::to_string(const calc_hit_list &arg_calc_hit_list ///< The calc_hit_list to describe
                             ) {
	return "calc_hit_list["
		+ ::std::to_string( arg_calc_hit_list.size() )
		+ "hits]";
}

/// \brief Insert a description of the specified calc_hit_list into the specified ostream
///
/// \relates calc_hit_list
ostream & cath::rslv::operator<<(ostream             &arg_ostream,      ///< The ostream into which the description should be inserted
                                 const calc_hit_list &arg_calc_hit_list ///< The calc_hit_list to describe
                                 ) {
	arg_ostream << to_string( arg_calc_hit_list );
	return arg_ostream;
}

/// \brief Get the maximum stop residue of all the hits in the specified calc_hit_list
///        or none if the calc_hit_list is empty
///
/// \relates calc_hit_list
residx_opt cath::rslv::get_max_stop(const calc_hit_list &arg_calc_hit_list ///< The calc_hit_list to query
                                    ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if arg_full_hit_list's empty
	return arg_calc_hit_list.empty()
		? none
		: make_optional( get_stop_res_index( *max_element(
			arg_calc_hit_list,
			[] (const calc_hit &x, const calc_hit &y) {
				return x.get_stop_arrow() < y.get_stop_arrow();
		} ) ) );
}


/// \brief Get the best (ie maximum) score of all the hits in the specified calc_hit_list
///        or none if the calc_hit_list is empty
///
/// \relates calc_hit_list
resscr_opt cath::rslv::get_best_score(const calc_hit_list &arg_calc_hit_list ///< The calc_hit_list to query
                                      ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if arg_full_hit_list's empty
	return arg_calc_hit_list.empty()
		? none
		: make_optional( max_element(
			arg_calc_hit_list,
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
calc_hit_list::const_iterator cath::rslv::find_first_hit_stopping_at_or_after(const calc_hit_list &arg_calc_hit_list, ///< The calc_hit_list to query
                                                                              const res_arrow     &arg_res_arrow      ///< The boundary that the calc_hit must stop at or after
                                                                              ) {
	return lower_bound(
		arg_calc_hit_list,
		arg_res_arrow,
		[] (const calc_hit &h, const res_arrow &a) {
			return ( h.get_stop_arrow() < a );
		}
	);
}

/// \brief Find the first calc_hit in the specified calc_hit_list that stops after the specified residue boundary
///
/// Note: may return one-after-end iterator
///
/// \relates calc_hit_list
calc_hit_list::const_iterator cath::rslv::find_first_hit_stopping_after(const calc_hit_list &arg_calc_hit_list, ///< The calc_hit_list to query
                                                                        const res_arrow     &arg_res_arrow      ///< The boundary that the calc_hit must stop after
                                                                        ) {
	return upper_bound(
		arg_calc_hit_list,
		arg_res_arrow,
		[] (const res_arrow &a, const calc_hit &h) {
			return ( a < h.get_stop_arrow() );
		}
	);
}

/// \brief Get a range of indices corresponding to those hits in the specified calc_hit_list that stop in
///        the specified range
///
/// The range is inclusive at the end (ie the returned range may include hits that stop at arg_stop_arrow)
///
/// \relates calc_hit_list
integer_range<hitidx_t> cath::rslv::indices_of_hits_that_stop_in_range(const calc_hit_list &arg_calc_hit_list, ///< The calc_hit_list to query
                                                                       const res_arrow     &arg_start_arrow,   ///< The start arrow that the hits corresponding to the returned indices must stop after
                                                                       const res_arrow     &arg_stop_arrow     ///< The stop  arrow that the hits corresponding to the returned indices must stop before or at
                                                                       ) {
	return {
		numeric_cast<hitidx_t>( distance( cbegin( arg_calc_hit_list ), find_first_hit_stopping_after( arg_calc_hit_list, arg_start_arrow ) ) ),
		numeric_cast<hitidx_t>( distance( cbegin( arg_calc_hit_list ), find_first_hit_stopping_after( arg_calc_hit_list, arg_stop_arrow  ) ) )
	};
}