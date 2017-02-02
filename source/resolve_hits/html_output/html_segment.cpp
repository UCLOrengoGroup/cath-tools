/// \file
/// \brief The html_segment class definitions

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

#include "html_segment.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/debug_numeric_cast.hpp"
#include "exception/invalid_argument_exception.hpp"

#include <string>

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::algorithm::starts_with;
using boost::none;
using std::string;

/// \brief Get an HTML span string to represent some aspect of a segment
string html_segment::get_html_string(const res_arrow          &arg_start,           ///< The start of the segment to render
                                     const res_arrow_opt      &arg_stop,            ///< The stop of the segment to render (or none for a boundary)
                                     const string             &arg_css_class,       ///< The CSS classes with which the HTML span should be marked
                                     const str_str_pair_vec   &arg_data_key_values, ///< A set of key/value pairs to be inserted as data attributes in the span (keys are prefixed with "data-" if not already)
                                     const display_colour_opt &arg_border_colour,   ///< The colour with which the border should be rendered
                                     const display_colour     &arg_fill_colour,     ///< The colour with which to fill the pill
                                     const size_t             &arg_full_seq_length, ///< The length of the full sequence on which this hit appears
                                     const pill_rounding      &arg_pill_rounding    ///< Which side of the pill (or neither/both) should be rounded
                                     ) {
	const double length_mult = 100.0 / debug_numeric_cast<double>( arg_full_seq_length );
	return R"(<span class=")"
		+ arg_css_class
		+ R"(" )"
		+ join(
			arg_data_key_values
				| transformed( [] (const str_str_pair &x) {
					return
						( starts_with( x.first, "data-" ) ? ""s : "data-" )
						+ x.first  + R"(=")"
						+ x.second + R"(")";
				} ),
			" "
		)
		+ R"( style="background-color: #)"
		+ hex_string_of_colour( arg_fill_colour )
		+ ";" 
		+ (
			arg_border_colour
			? ( " border-color: #" + hex_string_of_colour( *arg_border_colour ) + ";" )
			: ""s
		)
		+ " left: "
		+ ::std::to_string( length_mult * debug_numeric_cast<double>( arg_start.res_after() - 1 ) )
		+ R"(%;)"
		+ (
			arg_stop
				? " width: "
					+ ::std::to_string( length_mult * debug_numeric_cast<double>( *arg_stop - arg_start ) )
					+ R"(%;)"
				: ""s
		)
		+ (
			( arg_pill_rounding == pill_rounding::NEITHER || arg_pill_rounding == pill_rounding::RIGHT_ONLY )
			? " border-top-left-radius: 0; border-bottom-left-radius: 0;"s
			: ""s
		)
		+ (
			( arg_pill_rounding == pill_rounding::NEITHER || arg_pill_rounding == pill_rounding::LEFT_ONLY )
			? " border-top-right-radius: 0; border-bottom-right-radius: 0;"s
			: ""s
		)
		+ R"("></span>)";
}

/// \brief Get a resolved boundary HTML string
string html_segment::get_resolve_boundary_html_string(const res_arrow      &arg_point,          ///< The location of the arrow
                                                      const display_colour &arg_colour,         ///< The colour in which this boundary should be rendered
                                                      const size_t         &arg_full_seq_length ///< The length of the full sequence on which this hit appears
                                                      ) {
	return get_html_string(
		arg_point,
		none,
		"crh-hit-boundary",
		{},
		display_colour::BLACK,
		darken_by_fraction( arg_colour, 0.60 ),
		arg_full_seq_length
	);
}

/// \brief Get a grey back pill HTML string
///
/// This appears underneath the hit and shows the full extent of the hit
/// in cases of parts missing from results with overlapping hits
string html_segment::get_grey_back_html_string() const {
	return get_html_string(
		start,
		stop,
		"crh-hit-pill-tail",
		data_key_values,
		none,
		display_colour::WHITE,
		full_seq_length
	);
}

/// \brief Get a lightened back pill HTML string
///
/// This is used to show the region of the full segment, before trimming
string html_segment::get_lightened_back_html_string() const {
	return get_html_string(
		resolved_start.value_or( start ),
		resolved_stop .value_or( stop  ),
		"crh-hit-pill-ends",
		data_key_values,
		darken_by_fraction ( colour, 0.60 ),
		lighten_by_fraction( colour, 0.75 ),
		full_seq_length,
		resolved_start
			? ( resolved_stop ? pill_rounding::NEITHER   : pill_rounding::RIGHT_ONLY )
			: ( resolved_stop ? pill_rounding::LEFT_ONLY : pill_rounding::BOTH       )
	);
}

/// \brief Get a strong front pill HTML string for use in a full result
string html_segment::get_full_result_html_string() const {
	return get_html_string(
		resolved_start.value_or( start ),
		resolved_stop .value_or( stop  ),
		"crh-hit-pill-core",
		data_key_values,
		darken_by_fraction ( colour, 0.60 ),
		colour,
		full_seq_length,
		resolved_start
			? ( resolved_stop ? pill_rounding::NEITHER   : pill_rounding::RIGHT_ONLY )
			: ( resolved_stop ? pill_rounding::LEFT_ONLY : pill_rounding::BOTH       )
	);
}

/// \brief Get a strong front pill HTML string
///
/// This is used to show the region of the trimmed segment
string html_segment::get_strong_front_html_string() const {
	if ( ! trimmed_start || ! trimmed_stop ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot render strong front of segment in HTML if there is no trimmed start/stop"));
	}
	return get_html_string(
		*trimmed_start,
		*trimmed_stop,
		"crh-hit-pill-core",
		data_key_values,
		darken_by_fraction ( colour, 0.60 ),
		colour,
		full_seq_length
	);
}

/// \brief Get the list of HTML span strings to represent this segment
str_vec html_segment::get_all_span_html_strs(const bool &arg_do_layers ///< Whether or not render multiple layers or a single layer as in a full result
                                             ) const {
	if ( ! arg_do_layers ) {
		return { get_full_result_html_string() };
	}
	if ( ! trimmed_start || ! trimmed_stop ) {
		if ( resolved_start || resolved_stop ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("A segment for HTML rendering shouldn't have result-resolved-boundaries if it hasn't got a trimmed core"));
		}
		return { get_lightened_back_html_string() };
	}
	str_vec span_strs;
	span_strs.push_back( get_grey_back_html_string     () );
	if ( resolved_start ) {
		span_strs.push_back( get_resolve_boundary_html_string( *resolved_start, colour, full_seq_length ) );
	}
	if ( resolved_stop ) {
		span_strs.push_back( get_resolve_boundary_html_string( *resolved_stop,  colour, full_seq_length ) );
	}
	span_strs.push_back( get_lightened_back_html_string() );
	span_strs.push_back( get_strong_front_html_string  () );
	return span_strs;
}
