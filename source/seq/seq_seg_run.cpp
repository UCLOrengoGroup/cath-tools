/// \file
/// \brief The seq_seg_run class definitions

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

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/boost_addenda/range/indices.hpp"
#include "seq/seq_seg_run.hpp"

#include <string>

using namespace cath::common;
using namespace cath::seq;
using namespace std::literals::string_literals;

using boost::adaptors::transformed;
using boost::algorithm::all_of;
using boost::algorithm::join;
using std::ostream;
using std::string;

/// \brief Generate a string describing the segments of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates seq_seg_run
string cath::seq::get_segments_string(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run whose segments should be described
                                      ) {
	return join(
		indices( arg_seq_seg_run.get_num_segments() )
			| boost::adaptors::transformed( [&] (const size_t &x) {
				return ::std::to_string( get_start_res_index_of_segment( arg_seq_seg_run, x ) )
				     + "-"
				     + ::std::to_string( get_stop_res_index_of_segment ( arg_seq_seg_run, x ) );
			} ),
		","
	);
}

/// \brief Generate a string describing the segments of the specified string or an empty string if none
///
/// If not none, the result is the same as returned from get_segments_string() prefixed by '/'
///
/// \relates seq_seg_run
string cath::seq::get_segments_suffix_string(const seq_seg_run_opt &arg_seq_seg_run ///< The seq_seg_run whose segments should be described
                                             ) {
	return arg_seq_seg_run ? ( "/" + get_segments_string( *arg_seq_seg_run ) ) : "";
}

/// \brief Generate a string describing the specified seq_seg_run
///
/// \relates seq_seg_run
string cath::seq::to_string(const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to describe
                            ) {
	return "seq_seg_run[" + get_segments_string( arg_seq_seg_run ) + "]";
}

/// \brief Generate a string describing the specified seq_seg_run
///
/// \relates seq_seg_run
string cath::seq::to_string(const seq_seg_run_opt &arg_seq_seg_run_opt ///< The seq_seg_run_opt to describe
                            ) {
	return "seq_seg_run["
		+ (
			arg_seq_seg_run_opt
				? get_segments_string( *arg_seq_seg_run_opt )
				: "WCD"
		)
		+ "]";
}

/// \brief Insert a description of the specified seq_seg_run into the specified ostream
///
/// \relates seq_seg_run
ostream & cath::seq::operator<<(ostream           &arg_os,         ///< The ostream into which the description should be inserted
                                const seq_seg_run &arg_seq_seg_run ///< The seq_seg_run to describe
                                ) {
	arg_os << to_string( arg_seq_seg_run );
	return arg_os;
}

/// \brief Return whether the two specified seq_seg_runs are identical
///
/// \relates seq_seg_run
bool cath::seq::operator==(const seq_seg_run &arg_lhs, ///< The first  seq_seg_run to compare
                           const seq_seg_run &arg_rhs  ///< The second seq_seg_run to compare
                           ) {
	return (
		( arg_lhs.get_num_segments() == arg_rhs.get_num_segments() )
		&&
		all_of(
			indices( arg_lhs.get_num_segments() ),
			[&] (const size_t &seg_ctr) {
				return (
					arg_lhs.get_start_arrow_of_segment( seg_ctr ) == arg_rhs.get_start_arrow_of_segment( seg_ctr )
					&&
					arg_lhs.get_stop_arrow_of_segment ( seg_ctr ) == arg_rhs.get_stop_arrow_of_segment ( seg_ctr )
				);
			}
		)
	);
}
