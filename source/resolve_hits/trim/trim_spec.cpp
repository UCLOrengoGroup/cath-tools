/// \file
/// \brief The trim_spec class definitions

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

#include "trim_spec.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/program_options/validator.hpp"

#include <string>

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::seq;

using boost::adaptors::transformed;
using boost::algorithm::is_any_of;
using boost::algorithm::join;
using boost::any;
using boost::lexical_cast;
using boost::optional;
using boost::adaptors::filtered;
using std::istream;
using std::ostream;
using std::stoul;
using std::string;

/// \brief Generate a short string describing the specified trim_spec
///        that can, for example, be used to describe the default trim_spec
///        in the program options description
///
/// \relates trim_spec
string cath::rslv::to_options_string(const trim_spec &arg_trim_spec ///< The trim_spec to describe
                                     ) {
	return
		  ::std::to_string( arg_trim_spec.get_full_length()    )
		+ "/"
		+ ::std::to_string( arg_trim_spec.get_total_trimming() );
}

/// \brief Generate a string describing the specified trim_spec
///
/// \relates trim_spec
string cath::rslv::to_string(const trim_spec &arg_trim_spec ///< The trim_spec to describe
                             ) {
	return "trim_spec[full_length: "
		+ ::std::to_string( arg_trim_spec.get_full_length()    )
		+ "; total_trimming: "
		+ ::std::to_string( arg_trim_spec.get_total_trimming() )
		+ "]";
}

/// \brief Insert a description of the specified trim_spec into the specified ostream
///
/// \relates trim_spec
ostream & cath::rslv::operator<<(ostream         &arg_os,       ///< The ostream into which the description should be inserted
                                 const trim_spec &arg_trim_spec ///< The trim_spec to describe
                                 ) {
	arg_os << to_string( arg_trim_spec );
	return arg_os;
}

/// \brief Attempt to parse a trim_spec from the specified string (eg "30/10")
///
/// \relates trim_spec
trim_spec cath::rslv::parse_trim_spec(const string &arg_trim_spec_str ///< The string from which the trim_spec should be parsed
                                      ) {
	const auto parts = split_build<str_vec>( arg_trim_spec_str, is_any_of( "/" ) );
	if ( parts.size() != 2 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Cannot find two parts separated by a '/'' in trim specification"));
	}
	try {
		const trim_spec the_result{
			lexical_cast<residx_t>( parts.front() ),
			lexical_cast<residx_t>( parts.back()  )
		};
		return the_result;
	}
	catch (...) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Problem making trim_spec"));
	}
}

/// \brief Extract into the specified trim_spec from the specified stream
///
/// \relates trim_spec
istream & cath::rslv::operator>>(istream   &arg_is,       ///< The stream from which the value should be extracted
                                 trim_spec &arg_trim_spec ///< The trim_spec to populate from the specified stream
                                 ) {
	string input_string;
	arg_is >> input_string;
	arg_trim_spec = parse_trim_spec( input_string );
	return arg_is;
}

/// \brief Extract into the specified trim_spec from the specified stream
///
/// \relates seq_seg
///
/// \alsorelates trim_spec
string cath::rslv::to_possibly_trimmed_simple_string(const seq_seg       &arg_seq_seg,      ///< The seq_seg to describe
                                                     const trim_spec_opt &arg_trim_spec_opt ///< The optional specification describing the possible trimming
                                                     ) {
	return to_simple_string( trim_seq_seg_copy( arg_seq_seg, arg_trim_spec_opt ) );
}

/// \brief Generate a string describing the segments of the specified segments
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates full_hit
seq_seg_vec cath::rslv::get_segments(const seq_seg_vec   &arg_segs,         ///< The segments to be described
                                     const trim_spec_opt &arg_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                     ) {
	return transform_build<seq_seg_vec>(
		arg_segs,
		[&] (const seq_seg &x) {
			return trim_seq_seg_copy( x, arg_trim_spec_opt );
		}
	);
}

/// \brief Generate a string describing the segments of the specified segments
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments.
///
/// \relates full_hit
string cath::rslv::get_segments_string(const seq_seg_vec   &arg_segs,         ///< The segments to be described
                                       const trim_spec_opt &arg_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                       ) {
	return get_segments_string( get_segments( arg_segs, arg_trim_spec_opt ) );
}

/// \brief Generate a string describing the segments of the specified optional segments
///
/// This is the numbers of the start/stop residues, separated by a '-' between the start and stop
/// and by a ',' between segments. Values of none are excluded from the output.
///
/// \relates full_hit
string cath::rslv::get_segments_string(const seq_seg_opt_vec &arg_segs,         ///< The segments to be described
                                       const trim_spec_opt   &arg_trim_spec_opt ///< An optional trim_spec which may be used to specify trimming for the segments in the string
                                       ) {
	return get_segments_string(
		get_present_segments( arg_segs),
		arg_trim_spec_opt
	);
}

/// \brief Provide Boost program_options validation for trim_spec
///
/// \relates trim_spec
void cath::rslv::validate(any           &arg_value,         ///< The value to populate
                          const str_vec &arg_value_strings, ///< The string values to validate
                          trim_spec *, int) {
	arg_value = lex_castable_validator<trim_spec>::perform_validate( arg_value, arg_value_strings );
}
