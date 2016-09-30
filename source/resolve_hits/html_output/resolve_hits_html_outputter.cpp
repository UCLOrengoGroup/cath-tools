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

#include "resolve_hits_html_outputter.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include "common/debug_numeric_cast.h"
#include "display_colour/display_colour_gradient.h"
#include "resolve_hits/hit_list.h"
#include "resolve_hits/trim/trim_spec.h"

#include <string>

using namespace cath;
using namespace cath::rslv;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::format;
using boost::irange;
using std::string;

/// \brief Generate HTML to describe the specified hit_list with the specified trim_spec applied
string resolve_hits_html_outputter::output_html(const hit_list  &arg_hit_list, ///< The hit_list to describe
                                                const trim_spec &arg_trim_spec ///< The trim_spec that should be shown on the hits
                                                ) {
	const auto gradient = display_colour_gradient{
		{
			display_colour::LIGHT_RED,
			display_colour::YELLOW,
			display_colour::LIGHT_GREEN,
		},
		255
	};
	const auto best_score = get_best_score( arg_hit_list );
	const auto max_stop   = get_max_stop  ( arg_hit_list );
	return R"(<!DOCTYPE html>
<html lang="en">

<head>
	<meta charset="utf-8">
	<title>cath-resolve-hits</title>
</head>

<body>

<div id="domain-scan-results" class="domain">

<table class="table table-condensed">

<thead></thead>
<tbody><tr></tr></tbody>
<th>Match</th>
<th>Length</th>
<th>Regions</th>
<th>Architecture</th>
<th>Score</th>
<tbody>

)"
	// + ::std::to_string( arg_hit_list.size() )
	// + " \n "
	+ join(
		irange( 0_z, arg_hit_list.size() )
			| transformed( [&] (const size_t &x) {
				return output_html_fragment(
					arg_hit_list[ x ],
					arg_hit_list[ x ].get_label( arg_hit_list.get_labels() ),
					get_colour_of_fraction(
						gradient,
						arg_hit_list[ x ].get_score() / *best_score
					),
					arg_trim_spec,
					*max_stop
				);
			} ),
		"\n"
	)
	+ R"(

</table>

</div>

</body>
</html>
)";
}

// /// \brief Generate an HTML fragment to describe the specified hit with the specified trim_spec applied
string resolve_hits_html_outputter::output_html_fragment(const hit            &arg_hit,            ///< The hit to desribe
                                                         const string         &arg_label,          ///< The match ID / label associated with the hit
                                                         const display_colour &arg_colour,         ///< The colour in which the hit should be rendered
                                                         const trim_spec      &arg_trim_spec,      ///< ///< The trim_spec that should be shown on the hit
                                                         const size_t         &arg_sequence_length ///< The length of the full sequence on which this hit appears
                                                         ) {
	const double length_mult = 100.0 / debug_numeric_cast<double>( arg_sequence_length );

	// For strictly-worse rows, can set: background-color: #ddd; color: #999;
	return R"(

<tr>
	<td>
		)"
		+ arg_label
		+ R"(
	</td>
	<td>
		)"
		+ ::std::to_string( get_total_length( arg_hit ) )
		+ R"(
	</td>
	<td>
		<div class="scan-result-regions">
			)"
		+ join(
			get_hit_segs( arg_hit )
				| transformed( [] (const hit_seg &x) {
					return R"(<span class="scan-result-hsp-region" style="margin-right: 5px;">)"
						+ ::std::to_string( get_start_res_index( x ) )
						+ "-"
						+ ::std::to_string( get_stop_res_index ( x ) )
						+ "</span>";
				} ),
			""
		)
		+ R"(
		</div>
	</td>
	<td>
		<div class="scan-result-figure" style="width: 500px; height: 3px; position: relative; background-color: #ddd;">
)"
		+ join(
			get_hit_segs( arg_hit )
				| transformed( [&] (const hit_seg &x) {
					const auto  trimmed_start_stop = trim_copy_start_stop(
						arg_trim_spec,
						get_start_res_index( x ),
						get_stop_res_index ( x )
					);
					const auto &trimmed_start      = trimmed_start_stop.first;
					const auto &trimmed_stop       = trimmed_start_stop.second;
					const auto  trimmed_length     = trimmed_stop + 1 - trimmed_start;

// 					return R"(		<span class="scan-result-hsp" style="height: 7px; margin-top: -3px; )"
// 						"border-style: solid; border-color: #"
// 						+ hex_string_of_colour( darken_by_fraction ( arg_colour, 0.6 ) )
// 						+ "; border-width: 1px; position: absolute; left: "
// 						+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_start_res_index( x ) -  1 ) )
// 						+ R"(%; width: )"
// 						+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_length         ( x )      ) )
// 						+ R"(%; background-color: #)"
// 						+ hex_string_of_colour( arg_colour )
// 						+ R"(;"></span>
// 		<span class="scan-result-hsp" style="height: 7px; margin-top: -2px; border-style: solid; border-color: #)"
// 						+ hex_string_of_colour( darken_by_fraction ( arg_colour, 0.15 ) )
// 						+ "; border-width: 0px; border-right-width: 1px; position: absolute; left: "
// 						+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_start_res_index( x ) -  1 ) )
// 						+ R"(%; width: )"
// 						+ ::std::to_string( length_mult * debug_numeric_cast<double>( 15 ) )
// 						+ R"(%; background-color: #)"
// 						+ hex_string_of_colour( lighten_by_fraction( arg_colour, 0.85 ) )
// 						+ R"(;\"></span>
// 		<span class="scan-result-hsp" style="height: 7px; margin-top: -2px; border-style: solid; border-color: #)"
// 						+ hex_string_of_colour( darken_by_fraction ( arg_colour, 0.15 ) )
// 						+ "; border-width: 0px; border-left-width: 1px; position: absolute; left: "
// 						+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_stop_res_index( x ) - 15 ) )
// 						+ R"(%; width: )"
// 						+ ::std::to_string( length_mult * debug_numeric_cast<double>( 15 ) )
// 						+ R"(%; background-color: #)"
// 						+ hex_string_of_colour( lighten_by_fraction( arg_colour, 0.85 ) )
// 						+ R"(;\"></span>
// )";

					return R"(		<span class="scan-result-hsp" style="height: 5px; margin-top: -2px; )"
						"border-style: solid; border-color: #"
						+ hex_string_of_colour( darken_by_fraction( arg_colour, 0.65 ) )
						+ "; border-width: 1px; position: absolute; left: "
						+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_start_res_index( x ) - 1 ) )
						+ R"(%; width: )"
						+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_length         ( x )     ) )
						+ R"(%; background-color: #)"
						+ hex_string_of_colour( lighten_by_fraction( arg_colour, 0.75 ) )
						+ R"(;"></span>
		<span class="scan-result-hsp" style="height: 7px; margin-top: -3px; border-style: solid; border-color: #)"
						+ hex_string_of_colour( darken_by_fraction ( arg_colour, 0.65 ) )
						+ "; border-width: 1px; position: absolute; left: "
						+ ::std::to_string( length_mult * debug_numeric_cast<double>( trimmed_start - 1 ) )
						+ R"(%; width: )"
						+ ::std::to_string( length_mult * debug_numeric_cast<double>( trimmed_length    ) )
						+ R"(%; background-color: #)"
						+ hex_string_of_colour( arg_colour )
						+ ";\"></span>\n";


		// 			return R"(		<span class="scan-result-hsp" style="height: 7px; margin-top: -3px; )"
		// 				"border-style: solid; border-color: #"
		// 				+ hex_string_of_colour( darken_by_fraction( arg_colour, 0.65 ) )
		// 				+ "; border-width: 1px; position: absolute; left: "
		// 				+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_start_res_index( x ) - 1 ) )
		// 				+ R"(%; width: )"
		// 				+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_length         ( x )     ) )
		// 				+ R"(%; background-color: #)"
		// 				+ hex_string_of_colour( lighten_by_fraction( arg_colour, 0.85 ) )
		// 				+ R"(;"></span>
		// <span class="scan-result-hsp" style="height: 7px; margin-top: -2px; border-style: solid; border-color: #)"
		// 				+ hex_string_of_colour( darken_by_fraction ( arg_colour, 0.15 ) )
		// 				+ "; border-width: 0px; border-left-width: 1px; border-right-width: 1px; position: absolute; left: "
		// 				+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_start_res_index( x ) + 14 ) )
		// 				+ R"(%; width: )"
		// 				+ ::std::to_string( length_mult * debug_numeric_cast<double>( get_length         ( x ) - 30 ) )
		// 				+ R"(%; background-color: #)"
		// 				+ hex_string_of_colour( arg_colour )
		// 				+ ";\"></span>\n";


						// <span class="scan-result-hsp" style="height: 7px; margin-top: -2px; position: absolute; width:  20px; left:  10px; background-color: #5f5;"></span>
				} ),
			""
		)
			// <span class="scan-result-hsp" style="height: 7px; margin-top: -3px; border-style: solid; border-color: black; border-width: 1px; position: absolute; width: 135.3319px; left: 1.07066px; background-color: #efe;"></span>
			// <span class="scan-result-hsp" style="height: 7px; margin-top: -3px; border-style: solid; border-color: black; border-width: 1px; position: absolute; width: 87.7944px; left: 89.9358px; background-color: #fee;"></span>
			// <span class="scan-result-hsp" style="height: 7px; margin-top: -2px; position: absolute; width:  20px; left:  10px; background-color: #5f5;"></span>
			// <span class="scan-result-hsp" style="height: 7px; margin-top: -2px; position: absolute; width:  70px; left: 100px; background-color: #f55;"></span>
		+ R"(		</div>
	</td>
	<td>)"
		// + ::std::to_string( arg_hit.get_score() )
		+ ( format("%.3g") % arg_hit.get_score() ).str()
		+ R"(</td>
</tr>

)";
}

