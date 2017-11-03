/// \file
/// \brief The hits_input_format_tag definitions

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

#include "hits_input_format_tag.hpp"

#include <boost/algorithm/string/case_conv.hpp>

#include "common/exception/out_of_range_exception.hpp"
#include "common/program_options/validator.hpp"

#include <map>

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using boost::algorithm::to_lower;
using boost::any;
using std::istream;
using std::map;
using std::ostream;
using std::string;

/// \brief Getter for a map from name to hits_input_format_tag
///
/// \relates hits_input_format_tag
map<string, hits_input_format_tag> hits_input_format_tag_by_name::get() {
	map<string, hits_input_format_tag> result;
	for (const auto &x : all_hits_input_format_tags) {
		result.emplace( to_string( x ), x );
	}
	return result;
}


/// \brief Get a list of all the hits_input_format_tag names
str_vec all_hits_input_format_tag_names::get() {
	str_vec result;
	for (const auto &x : all_hits_input_format_tags) {
		result.emplace_back( to_string( x ) );
	}
	return result;
}

/// \brief Generate a string describing the specified hits_input_format_tag (ie its name)
///
/// \relates hits_input_format_tag
string cath::rslv::to_string(const hits_input_format_tag &arg_format_tag ///< The hits_input_format_tag to describe in a string
                             ) {
	switch ( arg_format_tag ) {
		case ( hits_input_format_tag::HMMER_DOMTBLOUT  ) : { return "hmmer_domtblout"  ; }
		case ( hits_input_format_tag::HMMSCAN_OUT      ) : { return "hmmscan_out"      ; }
		case ( hits_input_format_tag::HMMSEARCH_OUT    ) : { return "hmmsearch_out"    ; }
		case ( hits_input_format_tag::RAW_WITH_SCORES  ) : { return "raw_with_scores"  ; }
		case ( hits_input_format_tag::RAW_WITH_EVALUES ) : { return "raw_with_evalues" ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("hits_input_format_tag value not recognised"));
}

/// \brief Insert a description of the specified hits_input_format_tag into the specified stream
///
/// \relates hits_input_format_tag
ostream & cath::rslv::operator<<(ostream                     &arg_os,        ///< The ostream into which the hits_input_format_tag's description should be inserted
                                 const hits_input_format_tag &arg_format_tag ///< The hits_input_format_tag to describe in the ostream
                                 ) {
	arg_os << to_string( arg_format_tag );
	return arg_os;
}

/// \brief Extract into the specified hits_input_format_tag from the specified stream
///
/// \relates hits_input_format_tag
istream & cath::rslv::operator>>(istream               &arg_is,        ///< The stream from which the hits_input_format_tag should be extracted
                                 hits_input_format_tag &arg_format_tag ///< The hits_input_format_tag to populate from the specified stream
                                 ) {
	string input_string;
	arg_is >> input_string;
	to_lower( input_string );

	const auto all_layout_tags_by_name = hits_input_format_tag_by_name::get();
	arg_format_tag = all_layout_tags_by_name.at( input_string );
	return arg_is;
}

/// \brief Generate a string containing a description of the specified hits_input_format_tag
///
/// \relates hits_input_format_tag
string cath::rslv::description_of_input_format(const hits_input_format_tag &arg_format_tag ///< The hits_input_format_tag to describe
                                               ) {
	switch ( arg_format_tag ) {
		case ( hits_input_format_tag::HMMER_DOMTBLOUT  ) : { return "HMMER domtblout format (must assume all hits are continuous)"             ; }
		case ( hits_input_format_tag::HMMSCAN_OUT      ) : { return "HMMER hmmscan output format (can be used to deduce discontinuous hits)"   ; }
		case ( hits_input_format_tag::HMMSEARCH_OUT    ) : { return "HMMER hmmsearch output format (can be used to deduce discontinuous hits)" ; }
		case ( hits_input_format_tag::RAW_WITH_SCORES  ) : { return "\"raw\" format with scores"                                               ; }
		case ( hits_input_format_tag::RAW_WITH_EVALUES ) : { return "\"raw\" format with evalues"                                              ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("hits_input_format_tag value not recognised"));
}

/// \brief Provide Boost program_options validation for hits_input_format_tag
void cath::rslv::validate(any           &arg_value,         ///< The value to populate
                          const str_vec &arg_value_strings, ///< The string values to validate
                          hits_input_format_tag *, int) {
	arg_value = lex_castable_validator<hits_input_format_tag>::perform_validate( arg_value, arg_value_strings );
}
