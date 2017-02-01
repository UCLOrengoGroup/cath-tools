/// \file
/// \brief The sillitoe_chopping_format class definitions

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

#include "sillitoe_chopping_format.hpp"

#include "chopping/domain/domain.hpp"
#include "common/boost_addenda/make_string_ref.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/debug_numeric_cast.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/not_implemented_exception.hpp" // ***** TEMPORARY *****

#include <cctype>
#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::chop;
using namespace cath::common;

using boost::string_ref;
using std::find;
using std::isdigit;
using std::next;
using std::prev;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<chopping_format> sillitoe_chopping_format::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
bool sillitoe_chopping_format::do_represents_fragments() const {
	return false;
}

/// \brief TODOCUMENT
domain sillitoe_chopping_format::do_parse_domain(const string &arg_domain_chopping_string ///< TODOCUMENT
                                                 ) const {
	std::cerr << "domain_chopping_string is " << arg_domain_chopping_string << "\n";

	BOOST_THROW_EXCEPTION(not_implemented_exception("sillitoe_chopping_format::do_parse_domain()"));

	return domain( region_vec() );
}

/// \brief Parse a segment from the specified segment string
///
/// Example valid inputs: "7-232:K", "1B-99C:S"
region sillitoe_chopping_format::parse_segment(const string_ref &arg_segment_string ///< The string from which to parse the segment
                                               ) const {
	constexpr char   CHAIN_DELIM_COLON      = ':';
	constexpr char   RESIDUE_NAME_DELIM     = '-';
	constexpr size_t CHAIN_DELIM_NEG_OFFSET = 2;

	const auto length = arg_segment_string.length();
	if ( length < 2 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse sillitoe-chopping-format segment from a string of fewer than two characters"));
	}
	if ( arg_segment_string[ length - CHAIN_DELIM_NEG_OFFSET ] != CHAIN_DELIM_COLON ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse sillitoe-chopping-format segment from a string that doesn't have a colon in the penultimate character"));
	}

	const chain_label the_chain_label{ arg_segment_string.back() };

	if ( length == 2 ) {
		/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
		return region{ the_chain_label };
	}

	const auto begin_itr          = common::cbegin( arg_segment_string );
	const auto end_itr            = common::cend  ( arg_segment_string );
	const auto begin_plus_one_itr = next( begin_itr );
	const auto res_end_itr        = next( end_itr, 0 - static_cast<int>( CHAIN_DELIM_NEG_OFFSET ) );
	const auto dash_itr           = find( begin_plus_one_itr, res_end_itr, RESIDUE_NAME_DELIM );

	if ( dash_itr == res_end_itr ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse sillitoe-chopping-format segment from a string with no dash separating residue names"));
	}

	const auto dash_plus_one_itr = next( dash_itr );

	return {
		the_chain_label,
		parse_residue( make_string_ref( begin_itr,         dash_itr    ) ),
		parse_residue( make_string_ref( dash_plus_one_itr, res_end_itr ) )
	};
}

/// \brief Parse a residue from the specified residue string
///
/// Example valid inputs: "232", "99C"
residue_name sillitoe_chopping_format::parse_residue(const string_ref &arg_string_ref ///< The string from which to parse the residue
                                                     ) const {
	if ( arg_string_ref.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse sillitoe-chopping-format residue from an empty string"));
	}
	const auto begin_itr   = common::cbegin( arg_string_ref );
	const auto end_itr     = common::cend  ( arg_string_ref );
	if ( isdigit( arg_string_ref.back() ) ) {
		/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
		return residue_name{ stoi( string{ begin_itr, end_itr } ) };
	}
	else {
		if ( arg_string_ref.length() <= 1 ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse sillitoe-chopping-format residue from a string containing a single, non-numeric character"));
		}
		return {
			stoi( string{ begin_itr, prev( end_itr ) } ),
			arg_string_ref.back()
		};
	}
}
