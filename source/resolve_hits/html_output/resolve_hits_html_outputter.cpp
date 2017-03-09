/// \file
/// \brief The resolve_hits_html_outputter class definitions

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

#include "resolve_hits_html_outputter.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/sort_uniq_build.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/debug_numeric_cast.hpp"
#include "display_colour/display_colour_gradient.hpp"
#include "resolve_hits/calc_hit_list.hpp"
#include "resolve_hits/full_hit_list.hpp"
#include "resolve_hits/hit_resolver.hpp"
#include "resolve_hits/html_output/html_hit.hpp"
#include "resolve_hits/html_output/html_segment.hpp"
#include "resolve_hits/options/options_block/crh_html_options_block.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"
#include "resolve_hits/scored_hit_arch.hpp"
#include "resolve_hits/trim/trim_spec.hpp"

#include <string>

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::any_of;
using boost::algorithm::join;
using boost::algorithm::replace_all;
using boost::algorithm::to_lower_copy;
using boost::algorithm::to_upper_copy;
using boost::format;
using boost::irange;
using boost::make_optional;
using boost::none;
using std::get;
using std::make_pair;
using std::string;
using std::tie;
using std::tuple;
using std::vector;

string upper_first_lower_rest(const string &arg_string
                              ) {
	return arg_string.empty()
		? arg_string
		: (
			to_upper_copy( string{ arg_string.front() } )
			+ to_lower_copy( string{
				std::next( common::cbegin( arg_string ) ),
				           common::cend  ( arg_string )
			} )
		);
}

/// \brief Perform a dumb HTML escaping on the specified string
void dumb_html_escape(string &arg_string ///< The string to escape
                      ) {
	replace_all( arg_string, R"(&)", "&amp;"  );
	replace_all( arg_string, R"(")", "&quot;" );
	replace_all( arg_string, R"(')", "&apos;" );
	replace_all( arg_string, R"(<)", "&lt;"   );
	replace_all( arg_string, R"(>)", "&gt;"   );
}

/// \brief Perform a dumb HTML escaping on a copy of the specified string
string dumb_html_escape_copy(string arg_string ///< The source string to copy
                             ) {
	dumb_html_escape( arg_string );
	return arg_string;
}

/// \brief Get the row CSS class of the specified hit_row_context
///
/// \relates hit_row_context
string cath::rslv::row_css_class_of_hit_row_context(const hit_row_context &arg_row_context ///< The hit_row_context for which to return the appropriate row CSS class
                                                    ) {
	switch ( arg_row_context ) {
		case ( hit_row_context::HIGHLIGHT   ) : { return "crh-row-highlight"   ; }
		case ( hit_row_context::NORMAL      ) : { return "crh-row-normal"      ; }
		case ( hit_row_context::RESULT      ) : { return "crh-row-result"      ; }
		case ( hit_row_context::RESULT_FULL ) : { return "crh-row-result_full" ; }
		default : {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of hit_row_context not recognised whilst getting row_css_class_of_hit_row_context()"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}

}

/// \brief Get the first cell CSS class of the specified hit_row_context
///
/// \relates hit_row_context
string cath::rslv::first_cell_css_class_of_hit_row_context(const hit_row_context &arg_row_context ///< The hit_row_context for which to return the appropriate first cell CSS class
                                                           ) {
	return ( arg_row_context == hit_row_context::RESULT || arg_row_context == hit_row_context::HIGHLIGHT )
		? "crh-cell-first-highlight"s
		: "crh-cell-first-norm"s;
}

/// \brief Generate the HTML string for the total-score row of the table
string resolve_hits_html_outputter::total_score_row(const resscr_t &arg_total_score ///< The total score
                                                    ) {
	return R"(<tr class="crh-row-result">
	<td class="crh-cell crh-cell-data crh-cell-first-norm"></td>
	<td class="crh-cell crh-cell-data"></td>
	<td class="crh-cell crh-cell-data"><strong>= )"
	+ ( format("%.4g") % arg_total_score ).str()
	+ R"(</strong></td>
	<td class="crh-cell crh-cell-data"></td>
	<td class="crh-cell crh-cell-data"></td>
	<td class="crh-cell crh-cell-data"></td>
</tr>)";
}

/// \brief Generate the HTML string for the markers row (ie the residue numbers)
string resolve_hits_html_outputter::markers_row(const size_t  &arg_sequence_length,   ///< The length of the full sequence on which this full_hit appears
                                                const str_opt &arg_score_header_lbl   ///< The string for the score header or none if headers shouldn't be used
                                                ) {
	const double length_mult  = 100.0 / debug_numeric_cast<double>( arg_sequence_length );
	return R"(<tr class="crh-row-colhead">
	<td class="crh-cell crh-cell-first-norm">)" + ( arg_score_header_lbl ? "ID"s : ""s ) + R"(</td>
	<td class="crh-cell crh-cell-colhead-figure">
		<div class="crh-figure-div-empty">
)"
	+ join(
		irange( 0_z, arg_sequence_length + 1, step_for_length( arg_sequence_length ) )
			| transformed( [&] (const size_t &x) {
				const auto left_str = ::std::to_string( length_mult * debug_numeric_cast<double>( x ) );
				// 10000000 8 digits -> -20px;
				//        0 1 digit  ->  -2px;
				// ( 4 - 18 x  ) / 7
				const int num_digits            = debug_numeric_cast<int>( ::std::to_string( x ).length() );
				const string margin_left_px_str = ::std::to_string( ( 4 - ( 18 * num_digits ) ) / 7 );

				return R"(		<span class="crh-figure-marker-num" style="left: )"
					+ left_str
					+ R"(%; margin-left: )"
					+ margin_left_px_str
					+ R"(px;">)"
					+ ::std::to_string( x )
					+ R"(</span>
		<span class="crh-figure-marker-tick" style="left: )"
					+ left_str
					+ R"(%;"></span>)" + "\n";
			} ),
		""
	)
	+ R"(		</div>
	</td>
	<td class="crh-cell">)" + ( arg_score_header_lbl ? "Calc Score"s : ""s )   + R"(</td>
	<td class="crh-cell">)" + ( arg_score_header_lbl.value_or(         ""s ) ) + R"(</td>
	<td class="crh-cell">)" + ( arg_score_header_lbl ? "Regions"s    : ""s )   + R"(</td>
	<td class="crh-cell">)" + ( arg_score_header_lbl ? "Length"s     : ""s )   + R"(</td>
</tr>)";
}

/// \brief Merge the original segment boundaries with an optional set of resolved
///        boundaries, replacing originals with resolveds if they exist
hit_seg_opt_vec merge_opt_resolved_boundaries(const hit_seg_vec               &arg_segs,              ///< The original boundaries
                                              const seg_boundary_pair_vec_opt &arg_result_boundaries, ///< The optional resolved boundaries, where that has been required
                                              const crh_segment_spec          &arg_segment_spec       ///< The crh_segment_spec defining how the segments will be handled (eg trimmed) by the algorithm
                                              ) {
	return arg_result_boundaries
		? merge_boundaries( arg_segs, *arg_result_boundaries, arg_segment_spec )
		: transform_build<hit_seg_opt_vec>(
			arg_segs,
			[] (const hit_seg &x) { return make_optional( x ); }
		);
}

/// \brief Generate a HTML to describe the specified full_hit
string resolve_hits_html_outputter::hit_html(const html_hit         &arg_html_hit,        ///< The hit to render in HTML
                                             const crh_segment_spec &arg_segment_spec,    ///< The crh_segment_spec defining how the segments will be handled (eg trimmed) by the algorithm
                                             const size_t           &arg_sequence_length, ///< The length of the full sequence on which this full_hit appears
                                             const hit_row_context  &arg_row_context      ///< The context of the hits
                                             ) {
	const full_hit                  &the_full_hit          = arg_html_hit.hit_ref.get();
	const seg_boundary_pair_vec_opt &the_result_boundaries = arg_html_hit.result_boundaries;

	const auto resolved_boundaries = merge_opt_resolved_boundaries( the_full_hit.get_segments(), the_result_boundaries, arg_segment_spec );
	const auto boundaries_strs     = transform_build<str_vec>(
		resolved_boundaries
			| filtered( [] (const hit_seg_opt &x) { return static_cast<bool>( x ); } ),
		[] (const hit_seg_opt &x) { return to_simple_string( x ); }
	);

	return join(
		irange( 0_z, the_full_hit.get_segments().size() )
			| transformed( [&] (const size_t &x_idx) {
				const hit_seg &x = the_full_hit.get_segments()[ x_idx ];

				// Grab the result of applying the crh_segment_spec to the segment
				// (which may be boost::none if the segment is shorter than min-seg-length)
				const hit_seg_opt trimmed_seg_opt = apply_spec_to_seg_copy( x, arg_segment_spec );

				// Prepare batch index and hit index strings
				const string      batch_idx_str   = ::std::to_string( arg_html_hit.batch_idx + 1 );
				const string      hit_idx_str     = ::std::to_string( arg_html_hit.hit_idx   + 1 );

				return "\t\t"
					+ join(
						html_segment{
							x.get_start_arrow(),
							( trimmed_seg_opt ? make_optional( trimmed_seg_opt->get_start_arrow() ) : none ),
							( trimmed_seg_opt ? make_optional( trimmed_seg_opt->get_stop_arrow () ) : none ),
							x.get_stop_arrow (),
							( the_result_boundaries ? ( *the_result_boundaries )[ x_idx ].first  : none ),
							( the_result_boundaries ? ( *the_result_boundaries )[ x_idx ].second : none ),
							arg_html_hit.colour,
							// data fields, should allow something like:
							//
							//     $(".crh-hit-ill-core").each(function(n) {
							//         let d = $(n).data();
							//         d.crh-hit-id;  # batch3-hit11
							//         d.crh-hit-match-id; # 1cukA01
							//     })
							{
								make_pair( "crh-hit-id"s,         "batch" + batch_idx_str + "-hit" + hit_idx_str ),
								make_pair( "crh-hit-match-id"s,   the_full_hit.get_label() ),
								make_pair( "crh-hit-boundaries"s, join( boundaries_strs, ", " ) ),
								make_pair( "crh-seg-num"s,        ::std::to_string( x_idx + 1 ) ),
								make_pair( "crh-seg-boundaries"s, to_simple_string( resolved_boundaries[ x_idx ] ) ),
							},
							arg_sequence_length
						}.get_all_span_html_strs( arg_row_context != hit_row_context::RESULT_FULL ),
						"\n\t\t"
					)
					+ "\n";
			} ),
		""
	);
}

/// \brief Generate an HTML fragment to describe the specified full_hit with the specified crh_segment_spec applied
string resolve_hits_html_outputter::hits_row_html(const html_hit_vec     &arg_full_hits_data,    ///< The detail of the hit(s) to render in the row (maybe more than one if the row is combining all result hits)
                                                  const crh_segment_spec &arg_segment_spec,      ///< The crh_segment_spec defining how the segments will be handled (eg trimmed) by the algorithm
                                                  const crh_score_spec   &arg_score_spec,        ///< The crh_score_spec to use to calculate the crh-score
                                                  const size_t           &arg_sequence_length,   ///< The length of the full sequence on which this full_hit appears
                                                  const hit_row_context  &arg_row_context        ///< The context of the hits
                                                  ) {
	if ( arg_full_hits_data.empty() ) {
		return "";
	}

	const double    length_mult             = 100.0 / debug_numeric_cast<double>( arg_sequence_length );
	const bool      isnt_full_result        = arg_row_context != hit_row_context::RESULT_FULL;
	const full_hit &first_hit               = arg_full_hits_data.front().hit_ref.get();
	const auto     &first_result_boundaries = arg_full_hits_data.front().result_boundaries;
	const str_vec   boundaries_strs         = (
		isnt_full_result
			? transform_build<str_vec>(
				merge_opt_resolved_boundaries(
					first_hit.get_segments(),
					first_result_boundaries,
					arg_segment_spec
				) | filtered( [] (const hit_seg_opt &x) { return static_cast<bool>( x ); } ),
				[] (const hit_seg_opt &x) { return to_simple_string( x ); }
			)
			: str_vec{}
	);

	if ( isnt_full_result && arg_full_hits_data.size() != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot generate an HTML row of hits data for multiple hits when not generating a full result"));
	}

	// For strictly-worse rows, can set: background-color: #ddd; color: #999;
	return R"(<tr class=")" + row_css_class_of_hit_row_context( arg_row_context ) + R"(">
	<td class="crh-cell crh-cell-data )" + first_cell_css_class_of_hit_row_context( arg_row_context ) + R"(">
		)" + ( isnt_full_result ? dumb_html_escape_copy( first_hit.get_label() ) : "&nbsp;"s ) + R"(
	</td>
	<td class="crh-cell crh-cell-data">
		<div class="crh-figure-div-line">
)"
		+ join(
			irange( 0_z, arg_sequence_length + 1, step_for_length( arg_sequence_length ) )
				| transformed( [&] (const size_t &x) {
					return R"(		<span class="crh-figure-tick" style="left: )"
						+ ::std::to_string( length_mult * debug_numeric_cast<double>( x ) )
						+ R"(%;"></span>)" + "\n";
				} ),
			""
		)
		+ join(
			arg_full_hits_data
				| transformed( [&] (const html_hit &x) {
					return hit_html(
						x,
						arg_segment_spec,
						arg_sequence_length,
						arg_row_context
					);
				} ),
			""
		)
		+ R"(		</div>
	</td>
	<td class="crh-cell crh-cell-data">)"
		+ ( isnt_full_result && first_result_boundaries ? "+ "s : ""s )
		+ ( isnt_full_result ? ( format("%.3g") % get_crh_score( first_hit, arg_score_spec ) ).str() : ""s )
		+ R"(</td>
	<td class="crh-cell crh-cell-data">)"
		+ ( isnt_full_result ? get_score_string( first_hit, 4 ) : ""s )
		+ R"(</td>
	<td class="crh-cell crh-cell-data">)"
		+ (
			R"(
		<div class="scan-result-regions">
			)"
			+ join(
				boundaries_strs
					| transformed( [] (const string &x) {
						return R"(<span class="crh-chopping-region-text">)" + x + "</span>";
					} ),
				R"(,
			)"
			)
			+ R"(
		</div>
	)"
		)
		+ R"(</td>
	<td class="crh-cell crh-cell-data">)"
	+ (
		R"(
		)"
			+ ( isnt_full_result ? ::std::to_string( get_total_length( first_hit ) ) : ""s )
			+ R"(
	)"
	)
	+ R"(</td>
</tr>)";
}

/// \brief Generate the HTML prefix string
string resolve_hits_html_outputter::html_prefix() {
	return R"(<!DOCTYPE html>
<html lang="en">

<head>
	<meta charset="utf-8">
	<title>cath-resolve-hits</title>
</head>

<body class="crh-body">

<style>

)"
	+ css_string()
	+ R"(
</style>
)";
}

/* .crh-hit-pill-core:hover {
	border-width      : 2px;
	margin-left       : -1px;
	margin-top        : -4px;
} */

/* .crh-hit-pill-ends:hover {
	border-width      : 2px;
	margin-left       : -1px;
	margin-top        : -3px;
} */

/// \brief Generate the HTML suffix string
string resolve_hits_html_outputter::html_key() {
	return R"(

<div class="crh-key-div">

<h3 class="crh-key-header">Key</h3>
<ul>
	<li><strong>ends   </strong> : The pale ends of the segments represent the trimmed regions that cath-resolve-hits may let overlap. Segments shorter than the min-seg-length are shown all pale.</li>
	<li><strong>colours</strong> : The hits' colours represent their calculation score: red is 0; green is the best score in the input hits (which may nevertheless be a &quot;bad&quot; hit).</li>
	<li><strong>grey</strong>    : Grey hits represent hits with scores that have not met the --worst-permissible-[...] filters.</li>
</ul>
</div>
)";
}

/// \brief Generate the HTML suffix string
string resolve_hits_html_outputter::html_suffix() {
	return R"(

</body>

</html>
)";
}

/// \brief Generate the CSS string
string resolve_hits_html_outputter::css_string() {
	return R"(/* --- Start of simple reset --- */
body.crh-body {
	line-height       : 1.15;
}
body.crh-body table {
	border-collapse   : separate;
	border-spacing    : 4px;
}
/* --- End of simple reset --- */


.crh-body {
	background-color  : #e6e6e6;
	font-family       : Helvetica, Arial, sans-serif;
	font-size         : 90%;
}

.crh-results-wrapper {
	background-color  : white;
	border-radius     : 20px;
	display           : inline-block;
	margin            : 20px;
	padding           : 20px;
	padding-left      : 15px;
}

.crh-query-header {
	margin            : 0px;
}

.crh-table {
	text-align        : center;
}

.crh-table-subheading-first, .crh-table-subheading-later {
	font-weight       : bolder;
	padding-bottom    : 12px;
	text-align        : center;
}

.crh-table-subheading-later {
	padding-top       : 15px;
}

.crh-table-subheading-span {
	background-color  : royalblue;
	border-color      : darkslateblue;
	border-radius     : 7px;
	border-style      : solid;
	border-width      : 1px;
	box-shadow        : inset 1px 1px 1px rgba(255,255,255,0.2);
	color             : white;
	font-size         : large;
	padding           : 7px 12px 7px 12px;
	text-shadow       : 2px 2px 1px rgba(0,0,0,0.2);
}

.crh-table-subheading-span-symb {
	font-size         : 115%;
}

.crh-table-subheading-uparrow {
	color             : #e6e6e6;
	font-size         : 200%;
}

.crh-table-subheading-padding {
	display           : inline-block;
	width             : 75px;
}

.crh-row-normal {}

.crh-row-result {
	font-size         : 75%;
}

.crh-row-highlight, .crh-row-colhead {
	font-weight       : bold;
}

crh-row-colhead {}

.crh-row-soln-break {}

.crh-cell-soln-break {}

.crh-row-soln-break-uparrow {
	color             : #e6e6e6;
	font-size         : 133%;
}

.crh-row-soln-break-padding {
	display           : inline-block;
	width             : 17.5%;
}

.crh-cell {
	padding           : 0 3px 0 3px;
}

.crh-cell-data {
	font-size         : 85%;
}


.crh-cell-first-norm, .crh-cell-first-highlight {
	border-left-style : solid;
	border-left-width : 7px;
}

.crh-cell-first-norm {
	border-color      : transparent;
}

.crh-cell-first-highlight {
	border-color      : royalblue;
}



.crh-cell-colhead-figure {
	font-weight       : normal;
}


.crh-figure-marker-num, .crh-figure-marker-tick, .crh-figure-tick {
	position          : absolute;
	width             : 2px;
}


.crh-figure-marker-tick, .crh-figure-tick, .crh-figure-div-line {
	background-color  : #e6e6e6;
}


.crh-figure-marker-num {
	color             : #a2a2a2;
	font-size         : 65%;
	height            : 11px;
	margin-top        : -2px;
}

.crh-figure-marker-tick {
	height            : 8px;
	margin-left       : -1px;
	margin-top        : 12px;
}

.crh-figure-tick {
	height            : 11px;
	margin-left       : -1px;
	margin-top        : -4px;
}


.crh-figure-div-empty, .crh-figure-div-line {
	margin-left       : 10px;
	position          : relative;
	width             : 550px;
}

.crh-figure-div-empty {
	background-color  : transparent;
	height            : 20px;
}

.crh-figure-div-line {
	height            : 3px;
}

.crh-hit-pill-core, .crh-hit-pill-ends, .crh-hit-boundary {
	box-shadow        : 1px 1px 2px rgba(0, 0, 0, 0.2);
}

.crh-hit-pill-core, .crh-hit-pill-ends, .crh-hit-pill-tail {
	border-radius     : 5px;
	border-style      : solid;
	border-width      : 1px;
	position          : absolute;
}

.crh-hit-pill-ends, .crh-hit-pill-tail {
	height            : 5px;
	margin-top        : -2px;
}

.crh-hit-pill-tail {
	border-color      : #888;
	border-style      : dotted;
	z-index           : 2;
}

.crh-hit-boundary {
	height            : 11px;
	margin-top        : -4px;
	position          : absolute;
	width             : 2px;
	z-index           : 6;
}

.crh-hit-pill-ends {
	z-index           : 4;
}

.crh-hit-pill-core {
	height            : 7px;
	margin-top        : -3px;
	z-index           : 8;
}

.crh-chopping-region-text { }

.crh-key-div {
	font-size         : 75%;
}

.crh-key-header {
	padding-top       : 10px;
}

.crh-advert-div {
	font-size         : 70%;
	opacity           : 0.6;
	position          : relative;
}

.crh-advert-span {
	position          : absolute;
	right             : -6px;
	top               : -6px;
}

.crh-exclusion-note {
	color             : #777;
	font-size         : 80%;
	padding           : 15px 0 0 40px;
}
)";
}


/// \brief Generate the bumber of residues to use between residue-numbering-markers for a sequence of the specified length
///
/// Want step sizes to be simple: one of 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000...
///
/// Idea:
///  * get naive answer by dividing the length by the aim_min_num_steps (eg 1253 / 4 = 313.25)
///  * get the order of magnitude of that guess (eg for 313.25, the OoM is 100)
///  * choose one of a few simple multipliers (1, 2, 5) for that order of magnitude that multiplies up to
///    less than or equal to the original naive answer
///  * return the order of magnitude multiplied by that multiplier
size_t resolve_hits_html_outputter::step_for_length(const size_t &arg_sequence_length ///< The length of the sequence for which the markers' step-size should be calculated
                                                    ) {
	constexpr size_t aim_min_num_steps = 4;
	if ( arg_sequence_length <= 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to calculate step length for sequence of 0 length"));
	}
	if ( arg_sequence_length < aim_min_num_steps ) {
		return 1;
	}

	const double naive_guess    = debug_numeric_cast<double>( arg_sequence_length ) / debug_numeric_cast<double>( aim_min_num_steps );
	const size_t order_of_mgntd = debug_numeric_cast<size_t>( pow( 10, floor( log10( naive_guess ) ) ) );
	const double naive_mantissa = naive_guess / debug_numeric_cast<double>( order_of_mgntd );
	const size_t mult           = (
		( naive_mantissa < 2.0 ) ? 1_z :
		( naive_mantissa < 5.0 ) ? 2_z :
		                           5_z
	);

	return static_cast<size_t>( order_of_mgntd ) * mult;
}

/// \brief Generate HTML to describe the specified full_hit_list with the specified specs applied 
string resolve_hits_html_outputter::output_html(const string            &arg_query_id,         ///< The query ID
                                                full_hit_list          &&arg_full_hit_list,    ///< The full_hit_list to describe
                                                const crh_score_spec    &arg_score_spec,       ///< The crh_score_spec to use to calculate the crh-score
                                                const crh_segment_spec  &arg_segment_spec,     ///< The crh_segment_spec defining how the segments will be handled (eg trimmed) by the algorithm
                                                const crh_html_spec     &arg_html_spec,        ///< The specification for how to render the HTML
                                                const bool              &arg_output_head_tail, ///< Whether to include the head and tail (ie prefix and suffix) in the output
                                                const crh_filter_spec   &arg_filter_spec,      ///< The crh_filter_spec defining which input hits will be skipped by the algorithm
                                                const size_t            &arg_batch_index       ///< The index of the batch of hits being output (used to allow hits' HTML to have unique data attributes)
                                                ) {
	return output_html(
		arg_query_id,
		calc_hit_list{ std::move( arg_full_hit_list ), arg_score_spec, arg_segment_spec, arg_filter_spec },
		arg_score_spec,
		arg_segment_spec,
		arg_html_spec,
		arg_output_head_tail,
		arg_filter_spec,
		arg_batch_index
	);
}

/// \brief Generate HTML to describe the specified full_hit_list with the specified trim_spec applied
string resolve_hits_html_outputter::output_html(const string           &arg_query_id,         ///< The query ID
                                                const full_hit_list    &arg_full_hit_list,    ///< The full_hit_list to describe
                                                const crh_score_spec   &arg_score_spec,       ///< The crh_score_spec to use to calculate the crh-score
                                                const crh_segment_spec &arg_segment_spec,     ///< The crh_segment_spec defining how the segments will be handled (eg trimmed) by the algorithm
                                                const crh_html_spec    &arg_html_spec,        ///< The specification for how to render the HTML
                                                const bool             &arg_output_head_tail, ///< Whether to include the head and tail (ie prefix and suffix) in the output
                                                const crh_filter_spec  &arg_filter_spec,      ///< The crh_filter_spec defining which input hits will be skipped by the algorithm
                                                const size_t           &arg_batch_index       ///< The index of the batch of hits being output (used to allow hits' HTML to have unique data attributes)
                                                ) {
	return output_html(
		arg_query_id,
		calc_hit_list{ arg_full_hit_list, arg_score_spec, arg_segment_spec, arg_filter_spec },
		arg_score_spec,
		arg_segment_spec,
		arg_html_spec,
		arg_output_head_tail,
		arg_filter_spec,
		arg_batch_index
	);
}

/// \brief Generate HTML to describe the specified full_hit_list with the specified trim_spec applied
string resolve_hits_html_outputter::output_html(const string           &arg_query_id,         ///< The query ID
                                                const calc_hit_list    &arg_calc_hit_list,    ///< The calc_hit_list to describe
                                                const crh_score_spec   &arg_score_spec,       ///< The crh_score_spec to use to calculate the crh-score
                                                const crh_segment_spec &arg_segment_spec,     ///< The crh_segment_spec defining how the segments will be handled (eg trimmed) by the algorithm
                                                const crh_html_spec    &arg_html_spec,        ///< The specification for how to render the HTML
                                                const bool             &arg_output_head_tail, ///< Whether to include the head and tail (ie prefix and suffix) in the output
                                                const crh_filter_spec  &arg_filter_spec,      ///< The crh_filter_spec defining which input hits will be skipped by the algorithm
                                                const size_t           &arg_batch_index       ///< The index of the batch of hits being output (used to allow hits' HTML to have unique data attributes)
                                                ) {
	const auto  filtered_grey     = display_colour{ 0.666, 0.666, 0.666 };
	const auto &the_full_hit_list = arg_calc_hit_list.get_full_hits();
	const auto  best_result       = resolve_hits( arg_calc_hit_list );
	const auto  chosen_full_hits  = full_hit_list{ transform_build<full_hit_vec>(
		best_result.get_arch(),
		[&] (const calc_hit &x) {
			return the_full_hit_list[ x.get_label_idx() ];
		}
	) };
	const auto sorted_indices = sort_build<size_vec>(
		irange( 0_z, the_full_hit_list.size() ),
		[&] (const size_t &x, const size_t &y) {
			const auto neg_x_score = -get_crh_score( the_full_hit_list[ x ], arg_score_spec );
			const auto neg_y_score = -get_crh_score( the_full_hit_list[ y ], arg_score_spec );
			return
				tie( neg_x_score, x )
				<
				tie( neg_y_score, y );
		}
	);

	const auto gradient = display_colour_gradient{
		{
			display_colour::LIGHT_RED,
			display_colour::YELLOW,
			display_colour::LIGHT_GREEN,
		},
		255
	};
	const auto   best_crh_score = get_best_crh_score( the_full_hit_list, arg_score_spec );
	const auto   max_stop       = get_max_stop      ( the_full_hit_list );
	const size_t seq_length     = max_stop.value_or( 0_z );
	const auto   orig_score_str = the_full_hit_list.empty() ? "Score"s
	                                                        : upper_first_lower_rest( to_string( front( the_full_hit_list ).get_score_type() ) );

	// Variable to keep track of exclusions
	size_set non_soln_hit_indices;
	size_set excluded_non_soln_hit_indices;

	const string main_return_string = ( arg_output_head_tail ? html_prefix() : string{} ) + R"(
<div class="crh-results-wrapper">

<div class="crh-advert-div">
	<span class="crh-advert-span">
		Generated by <a href="http://cath-tools.readthedocs.io/en/latest/tools/cath-resolve-hits/">cath-resolve-hits</a>,
		one of the <a href="https://github.com/UCLOrengoGroup/cath-tools">cath-tools</a>
	</span>
</div>

<h3 class="crh-query-header">)" + dumb_html_escape_copy( arg_query_id ) + R"(</h3>
<table class="crh-table">

<tr>
	<td colspan="6" class="crh-table-subheading-first">
		<span class="crh-table-subheading-span">
			<strong class="crh-table-subheading-span-symb">&hearts;&nbsp;</strong>
			Resolved hits
		</span>
	</td>
</tr>

)"
	+ markers_row( seq_length, none )
	+ hits_row_html(
		transform_build<html_hit_vec>(
			best_result.get_arch(),
			[&] (const calc_hit &x) {
				const auto &the_index    = x.get_label_idx();
				const auto &the_full_hit = the_full_hit_list[ the_index ];
				return html_hit{
					the_full_hit,
					arg_batch_index,
					the_index,
					score_passes_filter( arg_filter_spec, the_full_hit.get_score(), the_full_hit.get_score_type() )
						? get_colour_of_fraction(
							gradient,
							get_crh_score( the_full_hit, arg_score_spec ) / *best_crh_score
						)
						: filtered_grey,
					make_optional( resolved_boundaries(
						the_full_hit,
						chosen_full_hits,
						arg_segment_spec
					) )
				};
			}
		),
		arg_segment_spec,
		arg_score_spec,
		seq_length,
		hit_row_context::RESULT_FULL
	)
	+ R"(<tr class="crh-row-soln-break">
	<td />
	<td class="crh-cell-soln-break">
		<span class="crh-row-soln-break-uparrow">&#11014;</span>
		<span class="crh-row-soln-break-padding">&nbsp;</span>
		<span class="crh-row-soln-break-uparrow">&#11014;</span>
		<span class="crh-row-soln-break-padding">&nbsp;</span>
		<span class="crh-row-soln-break-uparrow">&#11014;</span>
		<span class="crh-row-soln-break-padding">&nbsp;</span>
		<span class="crh-row-soln-break-uparrow">&#11014;</span>
	</td>
</tr>
)"
	+ join(
		best_result.get_arch()
			| transformed( [&] (const calc_hit &x) {
				const auto &the_index    = x.get_label_idx();
				const auto &the_full_hit = the_full_hit_list[ the_index ];
				return hits_row_html(
					{ html_hit{
						the_full_hit,
						arg_batch_index,
						the_index,
						score_passes_filter( arg_filter_spec, the_full_hit.get_score(), the_full_hit.get_score_type() )
							? get_colour_of_fraction(
								gradient,
								get_crh_score( the_full_hit, arg_score_spec ) / *best_crh_score
							)
							: filtered_grey,
						make_optional( resolved_boundaries(
							the_full_hit,
							chosen_full_hits,
							arg_segment_spec
						) )
					} },
					arg_segment_spec,
					arg_score_spec,
					seq_length,
					hit_row_context::RESULT
				);
			} ),
		"\n"
	)
	+ total_score_row( best_result.get_score() )
	+ R"(
<tr>
	<td colspan="6" class="crh-table-subheading-later">
		<span class="crh-table-subheading-uparrow">&#11014;</span>
		<span class="crh-table-subheading-padding">&nbsp;</span>
		<span class="crh-table-subheading-uparrow">&#11014;</span>
		<span class="crh-table-subheading-padding">&nbsp;</span>
		<span class="crh-table-subheading-span">
			<strong class="crh-table-subheading-span-symb">&#9921;&nbsp;</strong>
			Input hits
		</span>
		<span class="crh-table-subheading-padding">&nbsp;</span>
		<span class="crh-table-subheading-uparrow">&#11014;</span>
		<span class="crh-table-subheading-padding">&nbsp;</span>
		<span class="crh-table-subheading-uparrow">&#11014;</span>
	</td>
</tr>

)"
	+ markers_row( seq_length, make_optional( orig_score_str ) )
	+ "\n\n"
	+ join(
		sorted_indices
			| transformed( [&] (const size_t &x) {
				const auto            &hit_x     = the_full_hit_list[ x ];
				const bool             in_result = any_of( best_result.get_arch(), [&] (const calc_hit &y) { return y.get_label_idx() == x; } );
				const hit_row_context  context   = in_result ? hit_row_context::HIGHLIGHT
				                                             : hit_row_context::NORMAL;
				const bool             rejected  = ! score_passes_filter( arg_filter_spec, hit_x.get_score(), hit_x.get_score_type() );

				if ( rejected && arg_html_spec.get_exclude_rejected_hits() ) {
					return ""s;
				}
				if ( ! in_result ) {
					if ( ! contains( non_soln_hit_indices, x ) ) {
						if ( non_soln_hit_indices.size() >= arg_html_spec.get_max_num_non_soln_hits() ) {
							excluded_non_soln_hit_indices.insert( x );
							return ""s;
						}
						non_soln_hit_indices.insert( x );
					}
				}
				return hits_row_html(
					{ html_hit{
						hit_x,
						arg_batch_index,
						x,
						rejected
							? filtered_grey
							: get_colour_of_fraction(
								gradient,
								get_crh_score( hit_x, arg_score_spec ) / *best_crh_score
							),
						none
					} },
					arg_segment_spec,
					arg_score_spec,
					seq_length,
					context
				);
			} )
			| filtered( [] (const string &x) { return ! x.empty() ; } ),
		"\n"
	)
	+ R"(
</table>
)";
	// Have to fully calculate main_return_string and then concatenate extra stuff here
	// because otherwise C++'s order of evaluation rules don't guarantee that the code above
	// has run and potentially populated excluded_non_soln_hit_indices before it's used below.
	// I have seen this issue exhibited in code compiled by GCC (but not by Clang).
	return
		main_return_string
		+ (
			excluded_non_soln_hit_indices.empty()
				? ""s
				: (
					R"(<div class="crh-exclusion-note">...hiding another )"
					+ std::to_string( excluded_non_soln_hit_indices.size() )
					+ R"( non-solution results (current limit is )"
					+ std::to_string( arg_html_spec.get_max_num_non_soln_hits() )
					+ R"(; use <code>--)"
					+ crh_html_options_block::PO_MAX_NUM_NON_SOLN_HITS
					+ R"(</code> to change)</div>)"
				)
		)
		+ R"(
</div>
)"
		+ ( arg_output_head_tail ? html_suffix() : string{} );
}
