/// \file
/// \brief The full_hit functions definitions

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

#include "full_hit_fns.hpp"

#include "resolve_hits/file/alnd_rgn.hpp"
#include "resolve_hits/full_hit_fns.hpp"
#include "resolve_hits/full_hit_list_fns.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::make_optional;
using boost::none;
using std::string;

/// \brief Generate a string describing the segments of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates full_hit
hit_seg_vec cath::rslv::get_segments(const full_hit      &arg_full_hit,     ///< The full_hit whose segments should be described
                                     const trim_spec_opt &arg_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                     ) {
	return ::cath::rslv::get_segments( arg_full_hit.get_segments(), arg_trim_spec_opt );
}

/// \brief Generate a string describing the segments of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates full_hit
string cath::rslv::get_segments_string(const full_hit      &arg_full_hit,     ///< The full_hit whose segments should be described
                                       const trim_spec_opt &arg_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                       ) {
	return get_segments_string( arg_full_hit.get_segments(), arg_trim_spec_opt );
}

/// \brief Generate a string describing the specified full_hit
///
/// \relates full_hit
string cath::rslv::to_string(const full_hit             &arg_full_hit,         ///< The full_hit to describe
                             const hit_output_format    &arg_format,           ///< The format in which to generate the output
                             const string               &arg_prefix,           ///< A prefix string, typically used to put the query_id at the front. (Any non-empty string will have a space appended.)
                             const crh_segment_spec_opt &arg_segment_spec_opt, ///< An optional crh_segment_spec which can be used for including each full_hit's trimmed boundaries and resolved boundaries
                             const full_hit_list_opt    &arg_full_hits,        ///< An optional full_hit_list (from which the specified full_hit is drawn), which can be used for including the full_hit's resolved boundaries
                             const hit_boundary_output  &arg_boundary_output   ///< Whether to trim the boundaries before outputting them
                             ) {
	if ( arg_format != hit_output_format::JON && ! arg_prefix.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot specify prefix for any full_hit format other than JON"));
	}

	const auto trim_spec_opt            = get_trim_spec_opt( arg_segment_spec_opt );
	const auto trim_spec_if_trim_output =
		trim_spec_opt
		? make_optional( arg_boundary_output == hit_boundary_output::TRIMMED, *trim_spec_opt )
		: none;

	switch ( arg_format ) {
		case( hit_output_format::JON ) : {
			return arg_prefix
				+ ( arg_prefix.empty() ? ""s : " "s )
				+ arg_full_hit.get_label()
				+ " "
				+ get_score_string( arg_full_hit, 6 )
				+ " "
				+ get_segments_string( arg_full_hit, trim_spec_if_trim_output )
				+ (
					( arg_full_hits && arg_segment_spec_opt )
					? " " + get_all_resolved_segments_string( arg_full_hit, *arg_full_hits, *arg_segment_spec_opt )
					: ""s
				)
				+ (
					arg_full_hit.get_alnd_rgns_opt()
					? " " + to_string( *arg_full_hit.get_alnd_rgns_opt() )
					: ""s
				);
		}
		case( hit_output_format::CLASS ) : {
			return "full_hit["
				+ get_segments_string( arg_full_hit, trim_spec_if_trim_output )
				+ "; score: "
				+ get_score_string( arg_full_hit, 6 )
				+ "; label: \""
				+ arg_full_hit.get_label()
				+ "\"]";
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of full_hit_output_format not recognised in to_string() for full_hit"));
}


