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

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "cath/common/optional/make_optional_if.hpp"
#include "cath/resolve_hits/file/alnd_rgn.hpp"
#include "cath/resolve_hits/full_hit_fns.hpp"
#include "cath/resolve_hits/full_hit_list_fns.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::rslv;
using namespace ::cath::seq;

using ::boost::format;
using ::std::literals::string_literals::operator""s;
using ::std::make_optional;
using ::std::nullopt;
using ::std::string;

/// \brief Generate a string describing the segments of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates full_hit
seq_seg_vec cath::rslv::get_segments(const full_hit      &prm_full_hit,     ///< The full_hit whose segments should be described
                                     const trim_spec_opt &prm_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                     ) {
	return ::cath::rslv::get_segments( prm_full_hit.get_segments(), prm_trim_spec_opt );
}

/// \brief Generate a string describing the segments of the specified string
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates full_hit
string cath::rslv::get_segments_string(const full_hit      &prm_full_hit,     ///< The full_hit whose segments should be described
                                       const trim_spec_opt &prm_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                       ) {
	return get_segments_string( prm_full_hit.get_segments(), prm_trim_spec_opt );
}

/// \brief Generate a list of headers corresponding to the fields that will be generated for data as specified
///
/// \relates full_hit
str_vec cath::rslv::get_field_headers(const full_hit &prm_full_hit,                      ///< An example hit
                                      const bool     &has_prefix,                        ///< The data has a non-empty prefix (ie query-id)
                                      const bool     &prm_has_full_hits_and_segment_spec ///< The data has full_hits and segment_spec available
                                      ) {
	str_vec headers;
	headers.reserve(
		  ( has_prefix                         ? 1_z : 0_z )
		+ 3_z
		+ ( prm_has_full_hits_and_segment_spec ? 1_z : 0_z )
		+ prm_full_hit.get_extras_store().size()
	);
	if ( has_prefix ) {
		headers.push_back( full_hit::get_prefix_name()           );
	}
	headers.push_back( full_hit::get_label_name()            );
	headers.push_back( full_hit::get_score_name()            );
	headers.push_back( full_hit::get_segments_name()         );
	if ( prm_has_full_hits_and_segment_spec ) {
		headers.push_back( full_hit::get_resolved_name()         );
	}
	if ( get_first< hit_extra_cat::ALND_RGNS >( prm_full_hit.get_extras_store() ) ) {
		headers.push_back( to_string( hit_extra_cat::ALND_RGNS ) );
	}
	if ( get_first< hit_extra_cat::COND_EVAL >( prm_full_hit.get_extras_store() ) ) {
		headers.push_back( to_string( hit_extra_cat::COND_EVAL ) );
	}
	if ( get_first< hit_extra_cat::INDP_EVAL >( prm_full_hit.get_extras_store() ) ) {
		headers.push_back( to_string( hit_extra_cat::INDP_EVAL ) );
	}
	return headers;
}

/// \brief Generate a string describing the specified full_hit
///
/// \relates full_hit
string cath::rslv::to_string(const full_hit             &prm_full_hit,         ///< The full_hit to describe
                             const hit_output_format    &prm_format,           ///< The format in which to generate the output
                             const string               &prm_prefix,           ///< A prefix string, typically used to put the query_id at the front. (Any non-empty string will have a space appended.)
                             const crh_segment_spec_opt &prm_segment_spec_opt, ///< An optional crh_segment_spec which can be used for including each full_hit's trimmed boundaries and resolved boundaries
                             const full_hit_list_opt    &prm_full_hits,        ///< An optional full_hit_list (from which the specified full_hit is drawn), which can be used for including the full_hit's resolved boundaries
                             const hit_boundary_output  &prm_boundary_output   ///< Whether to trim the boundaries before outputting them
                             ) {
	if ( prm_format != hit_output_format::JON && ! prm_prefix.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot specify prefix for any full_hit format other than JON"));
	}

	const auto trim_spec_opt            = get_trim_spec_opt( prm_segment_spec_opt );
	const auto trim_spec_if_trim_output = if_then_optional(
		trim_spec_opt.has_value() && ( prm_boundary_output == hit_boundary_output::TRIMMED ),
		*trim_spec_opt
	);

	switch ( prm_format ) {
		case( hit_output_format::JON ) : {
			const str_opt  alnd_rgns_str_opt = get_first< hit_extra_cat::ALND_RGNS >( prm_full_hit.get_extras_store() );
			const doub_opt cond_eval_val_opt = get_first< hit_extra_cat::COND_EVAL >( prm_full_hit.get_extras_store() );
			const doub_opt indp_eval_val_opt = get_first< hit_extra_cat::INDP_EVAL >( prm_full_hit.get_extras_store() );
			return prm_prefix
				+ ( prm_prefix.empty() ? ""s : " "s )
				+ prm_full_hit.get_label()
				+ " "
				+ get_score_string( prm_full_hit, 6 )
				+ " "
				+ get_segments_string( prm_full_hit, trim_spec_if_trim_output )
				+ (
					( prm_full_hits && prm_segment_spec_opt )
					? " " + get_all_resolved_segments_string( prm_full_hit, *prm_full_hits, *prm_segment_spec_opt )
					: ""s
				)
				+ ( alnd_rgns_str_opt ? ( " " +                      *alnd_rgns_str_opt         ) : ""s )
				+ ( cond_eval_val_opt ? ( " " + ( format( "%.5g" ) % *cond_eval_val_opt ).str() ) : ""s )
				+ ( indp_eval_val_opt ? ( " " + ( format( "%.5g" ) % *indp_eval_val_opt ).str() ) : ""s );
		}
		case( hit_output_format::CLASS ) : {
			return "full_hit["
				+ get_segments_string( prm_full_hit, trim_spec_if_trim_output )
				+ "; score: "
				+ get_score_string( prm_full_hit, 6 )
				+ "; label: \""
				+ prm_full_hit.get_label()
				+ "\"]";
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of full_hit_output_format not recognised in to_string() for full_hit"));
}


