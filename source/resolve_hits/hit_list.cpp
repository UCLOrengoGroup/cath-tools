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
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/sub_range.hpp>

#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/file/open_fstream.h"
#include "common/type_aliases.h"

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
using boost::range::max_element;
using boost::sub_range;
using boost::token_compress_on;
using std::distance;
using std::ifstream;
using std::make_pair;
using std::ostream;
using std::string;
using std::vector;

/// \brief TODOCUMENT
///
/// \relates hit_list
hit_list cath::rslv::read_hit_list_from_file(const path &arg_file ///< TODOCUMENT
                                             ) {
	// cerr << "Starting\n" << std::flush;
	// const auto parse_start_time = std::chrono::high_resolution_clock::now();

	// const path the_file{ "/cath-tools/other_stuff/jon_dom_finder_gubbins/uniref100_sub_4135.faa.jon" };
	// const path the_file{ "/dev/shm/uniref100_sub_4135.faa.jon" };

	ifstream the_ifstream;
	open_ifstream( the_ifstream, arg_file );
	
//	size_t num_multi_segs_skipped = 0;

	vector<hit> the_hits;

// 	string line_string;
// 	while ( getline( the_ifstream, line_string ) ) {
// 		const auto line_parts = split_build<str_vec>( line_string, is_space(), token_compress_on );
// 		if ( line_parts.size() != 3 ) {
// 			BOOST_THROW_EXCEPTION(runtime_error_exception("Line doesn't contain three parts"));
// 		}
// 		const auto &label   = line_parts[ 0 ];
// 		const auto  score   = stof( line_parts[ 1 ] );
// 		const auto &seg_str = line_parts[ 2 ];
// //		if ( contains( seg_str, "," ) ) {
// //			++num_multi_segs_skipped;
// //			continue;
// //		}

// 		const auto seg_parts = split_build<str_vec>( seg_str, is_any_of( ",-" ) );
// 		if ( seg_parts.size() % 2 != 0 ) {
// 			BOOST_THROW_EXCEPTION(runtime_error_exception(
// 				"Segments str "
// 				+ seg_str
// 				+ " doesn't contain two parts"
// 			));
// 		}
// 		// if ( seg_parts.size() == 2 ) {
// //			the_hits.emplace_back(
// //				numeric_cast<residx_t>( stol( seg_parts[ 0 ] ) ),
// //				numeric_cast<residx_t>( stol( seg_parts[ 1 ] ) ),
// //				score,
// //				label
// //			);
// //		}
// //		else {
// 			the_hits.emplace_back(
// 				transform_build<hit_seg_vec>(
// 					irange( 0_z, seg_parts.size(), 2 ),
// 					[&] (const size_t &x) {
// 						return hit_seg{
// 							arrow_before_res( numeric_cast<residx_t>( stol( seg_parts[ x     ] ) ) ),
// 							arrow_after_res ( numeric_cast<residx_t>( stol( seg_parts[ x + 1 ] ) ) )
// 						};
// 					}
// 				),
// 				score,
// 				label
// 			);
// 		// }
// 	}

	string      label;
	resscr_t    score;
	string      parts_string;
	hit_seg_vec segments;
	while ( ! the_ifstream.eof() ) {

		the_ifstream >> label;
		the_ifstream >> score;
		the_ifstream >> parts_string;

		segments.clear();
		residx_t start;
		residx_t stop;
		auto str_ptr = &parts_string.front();
		int num_read = 0;
		while ( sscanf( str_ptr, "%u-%u%n", &start, &stop, &num_read ) != EOF ) {
			str_ptr += num_read;
			segments.emplace_back(
				arrow_before_res( start ),
				arrow_after_res ( stop  )
			);
			if ( *str_ptr != ',' ) {
				break;
			}
			++str_ptr;
		}

		the_hits.emplace_back(
			segments,
			score,
			label
		);
	}


	// BOOST_LOG_TRIVIAL( warning ) << "Parsed  " << the_hits.size()        << " single-segment hits";
//	BOOST_LOG_TRIVIAL( warning ) << "Skipped " << num_multi_segs_skipped << " multi-segment hits";
	// BOOST_LOG_TRIVIAL( warning )
	// 	<< "Took Parsing "
	// 	<< the_hits.size()
	// 	<< " hits took "
	// 	<< durn_to_seconds_string( std::chrono::high_resolution_clock::now() - parse_start_time );

	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return hit_list{ move( the_hits ) };
}

/// \brief TODOCUMENT
///
/// \relates hit_list
string cath::rslv::to_string(const hit_list &arg_hit_list ///< TODOCUMENT
                             ) {
	return "hit_list["
		+ ::std::to_string( arg_hit_list.size() )
		+ "hits]";
}

/// \brief TODOCUMENT
///
/// \relates hit_list
ostream & cath::rslv::operator<<(ostream   &arg_ostream,      ///< TODOCUMENT
                                 const hit_list &arg_hit_list ///< TODOCUMENT
                                 ) {
	arg_ostream << to_string( arg_hit_list );
	return arg_ostream;
}

/// \brief TODOCUMENT
///
/// \relates hit_list
residx_t cath::rslv::get_max_stop(const hit_list &arg_hit_list ///< TODOCUMENT
                                  ) {
	return get_stop_res_index( *max_element(
		arg_hit_list,
		[] (const hit &x, const hit &y) {
			return x.get_stop_arrow() < y.get_stop_arrow();
	} ) );
}

/// \brief TODOCUMENT
///
/// Note: may return one-after-end iterator
///
/// \relates hit_list
hit_list::const_iterator cath::rslv::find_first_hit_stopping_at_or_after(const hit_list  &arg_hit_list, ///< TODOCUMENT
                                                                         const res_arrow &arg_res_arrow ///< TODOCUMENT
                                                                         ) {
	return lower_bound(
		arg_hit_list,
		arg_res_arrow,
		[] (const hit &h, const res_arrow &a) {
			return ( h.get_stop_arrow() < a );
		}
	);
}

/// \brief TODOCUMENT
///
/// Note: may return one-after-end iterator
///
/// \relates hit_list
hit_list::const_iterator cath::rslv::find_first_hit_stopping_after(const hit_list  &arg_hit_list, ///< TODOCUMENT
                                                                   const res_arrow &arg_res_arrow ///< TODOCUMENT
                                                                   ) {
	return upper_bound(
		arg_hit_list,
		arg_res_arrow,
		[] (const res_arrow &a, const hit &h) {
			return ( a < h.get_stop_arrow() );
		}
	);
}

/// \brief TODOCUMENT
///
/// The range is inclusive (returned range may include hits that
/// stop at arg_start_arrow or arg_stop_arrow)
///
/// \relates hit_list
integer_range<hitidx_t> cath::rslv::indices_of_hits_that_stop_in_range(const hit_list  &arg_hit_list,    ///< TODOCUMENT
                                                                       const res_arrow &arg_start_arrow, ///< TODOCUMENT
                                                                       const res_arrow &arg_stop_arrow   ///< TODOCUMENT
                                                                       ) {
	return {
		numeric_cast<hitidx_t>( distance( cbegin( arg_hit_list ), find_first_hit_stopping_after( arg_hit_list, arg_start_arrow ) ) ),
		numeric_cast<hitidx_t>( distance( cbegin( arg_hit_list ), find_first_hit_stopping_after( arg_hit_list, arg_stop_arrow  ) ) )
	};
}