/// \file
/// \brief The hit_list class definitions

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

#include "hit_list.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/equal.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/upper_bound.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/spirit/include/qi.hpp>

#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/file/open_fstream.h"
#include "common/type_aliases.h"
#include "resolve_hits/read_and_resolve_mgr.h"

#include <chrono>
#include <chrono>
#include <fstream>
#include <string>

using namespace cath::common;
using namespace cath::rslv;

using boost::algorithm::is_any_of;
using boost::algorithm::is_space;
using boost::contains;
using boost::filesystem::path;
using boost::integer_range;
using boost::irange;
using boost::numeric_cast;
using boost::range::equal;
using boost::range::lower_bound;
using boost::range::max_element;
using boost::range::upper_bound;
using boost::spirit::qi::double_;
using boost::spirit::qi::omit;
using boost::spirit::qi::parse;
using boost::spirit::uint_;
using boost::sub_range;
using boost::token_compress_on;
using std::distance;
using std::find_if;
using std::ifstream;
using std::istream;
using std::make_pair;
using std::ostream;
using std::string;
using std::vector;


/// \brief Read a hit_list from the specified file
///
/// \relates hit_list
void cath::rslv::read_hit_list_from_file(read_and_resolve_mgr &arg_read_and_resolve_mgr, ///< The read_and_resolve_mgr to which each hit should be sent when it's read
                                         const path           &arg_file                  ///< The file from which to read the hits data
                                         ) {

	ifstream the_ifstream;
	open_ifstream( the_ifstream, arg_file );

	read_hit_list_from_istream( arg_read_and_resolve_mgr, the_ifstream );

	the_ifstream.close();
}

/// \brief Read a hit_list from the specified istream
///
/// \relates hit_list
void cath::rslv::read_hit_list_from_istream(read_and_resolve_mgr &arg_read_and_resolve_mgr, ///< The read_and_resolve_mgr to which each hit should be sent when it's read
                                            istream              &arg_istream               ///< The istream from which to read the hits data
                                            ) {
	if ( arg_read_and_resolve_mgr.is_active() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot read_hit_list_from_file() with a read_and_resolve_mgr that's already active"));
	}

	const auto is_space_char     = [ ] (const auto &x) { return ( ( x == ' ' ) || ( x == '\t' ) ); };
	const auto is_non_space_char = [ ] (const auto &x) { return ( ( x != ' ' ) && ( x != '\t' ) ); };
	const auto find_space        = [&] (const auto &b, const auto &e) {
		return find_if( b, e, is_space_char     );
	};
	const auto find_non_space    = [&] (const auto &b, const auto &e ) {
		return find_if( b, e, is_non_space_char );
	};

	resscr_t         score;
	hit_seg_vec      fragments;
	vector<residx_t> bounds;
	string           line;
	string           query_id;

	const auto bounds_pusher = [&] (const residx_t &x) { bounds.push_back( x ); };

	while ( getline( arg_istream, line ) ) {
		const auto line_begin_itr        = common::cbegin( line );
		const auto line_end_itr          = common::cend  ( line );

		const auto end_of_query_id_itr   = find_space    ( line_begin_itr,        line_end_itr );
		const auto begin_of_match_id_itr = find_non_space( end_of_query_id_itr,   line_end_itr );
		const auto end_of_match_id_itr   = find_space    ( begin_of_match_id_itr, line_end_itr );
		const auto begin_of_score_itr    = find_non_space( end_of_match_id_itr,   line_end_itr );

		if ( ! arg_read_and_resolve_mgr.is_active() || ! boost::range::equal( arg_read_and_resolve_mgr.get_query_id(), sub_range<const string>{ line_begin_itr, end_of_query_id_itr } ) ) {
			complete_if_active( arg_read_and_resolve_mgr );
			arg_read_and_resolve_mgr.set_query_id( string{ line_begin_itr, end_of_query_id_itr } );
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

		fragments.clear();
		for (const size_t &bound_ctr : irange( 1_z, bounds.size() - 1, 2_z ) ) {
			fragments.emplace_back(
				arrow_after_res ( bounds[ bound_ctr     ] ),
				arrow_before_res( bounds[ bound_ctr + 1 ]  )
			);
		}

		arg_read_and_resolve_mgr.add_hit(
			arrow_before_res( bounds.front() ),
			arrow_after_res ( bounds.back () ),
			fragments,
			score,
			string{ begin_of_match_id_itr, end_of_match_id_itr }
		);
	}

	complete_if_active( arg_read_and_resolve_mgr );
	arg_read_and_resolve_mgr.final_wait();
}

/// \brief Generate a string describing the specified hit_list
///
/// \relates hit_list
string cath::rslv::to_string(const hit_list &arg_hit_list ///< The hit_list to describe
                             ) {
	return "hit_list["
		+ ::std::to_string( arg_hit_list.size() )
		+ "hits]";
}

/// \brief Insert a description of the specified hit_list into the specified ostream
///
/// \relates hit_list
ostream & cath::rslv::operator<<(ostream        &arg_ostream, ///< The ostream into which the description should be inserted
                                 const hit_list &arg_hit_list ///< The hit_list to describe
                                 ) {
	arg_ostream << to_string( arg_hit_list );
	return arg_ostream;
}

/// \brief Get the maximum stop residue of all the hits in the specified hit_list
///
/// \relates hit_list
residx_t cath::rslv::get_max_stop(const hit_list &arg_hit_list ///< The hit_list to query
                                  ) {
	return get_stop_res_index( *max_element(
		arg_hit_list,
		[] (const hit &x, const hit &y) {
			return x.get_stop_arrow() < y.get_stop_arrow();
	} ) );
}

/// \brief Find the first hit in the specified hit_list that stops at or after the specified residue boundary
///
/// Note: may return one-after-end iterator
///
/// \relates hit_list
hit_list::const_iterator cath::rslv::find_first_hit_stopping_at_or_after(const hit_list  &arg_hit_list, ///< The hit_list to query
                                                                         const res_arrow &arg_res_arrow ///< The boundary that the hit must stop at or after
                                                                         ) {
	return lower_bound(
		arg_hit_list,
		arg_res_arrow,
		[] (const hit &h, const res_arrow &a) {
			return ( h.get_stop_arrow() < a );
		}
	);
}

/// \brief Find the first hit in the specified hit_list that stops after the specified residue boundary
///
/// Note: may return one-after-end iterator
///
/// \relates hit_list
hit_list::const_iterator cath::rslv::find_first_hit_stopping_after(const hit_list  &arg_hit_list, ///< The hit_list to query
                                                                   const res_arrow &arg_res_arrow ///< The boundary that the hit must stop after
                                                                   ) {
	return upper_bound(
		arg_hit_list,
		arg_res_arrow,
		[] (const res_arrow &a, const hit &h) {
			return ( a < h.get_stop_arrow() );
		}
	);
}

/// \brief Get a range of indices corresponding to those hits in the specified hit_list that stop in
///        the specified range
///
/// The range is inclusive at the end (ie the returned range may include hits that stop at arg_stop_arrow)
///
/// \relates hit_list
integer_range<hitidx_t> cath::rslv::indices_of_hits_that_stop_in_range(const hit_list  &arg_hit_list,    ///< The hit_list to query
                                                                       const res_arrow &arg_start_arrow, ///< The start arrow that the hits corresponding to the returned indices must stop after
                                                                       const res_arrow &arg_stop_arrow   ///< The stop  arrow that the hits corresponding to the returned indices must stop before or at
                                                                       ) {
	return {
		numeric_cast<hitidx_t>( distance( cbegin( arg_hit_list ), find_first_hit_stopping_after( arg_hit_list, arg_start_arrow ) ) ),
		numeric_cast<hitidx_t>( distance( cbegin( arg_hit_list ), find_first_hit_stopping_after( arg_hit_list, arg_stop_arrow  ) ) )
	};
}