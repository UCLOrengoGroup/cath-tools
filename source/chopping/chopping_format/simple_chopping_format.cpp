/// \file
/// \brief The simple_chopping_format class definitions

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

#include "simple_chopping_format.hpp"

#include "chopping/domain/domain.hpp"
#include "common/boost_addenda/make_string_ref.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/not_implemented_exception.hpp" // ***** TEMPORARY *****

#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::chop;
using namespace cath::common;

using boost::string_ref;
using std::distance;
using std::find;
using std::next;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<chopping_format> simple_chopping_format::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
bool simple_chopping_format::do_represents_fragments() const {
	return false;
}

/// \brief TODOCUMENT
domain simple_chopping_format::do_parse_domain(const string &arg_domain_chopping_string ///< TODOCUMENT
                                               ) const {
	std::cerr << "domain_chopping_string is " << arg_domain_chopping_string << "\n";

	BOOST_THROW_EXCEPTION(not_implemented_exception("simple_chopping_format::do_parse_domain()"));

	return domain( region_vec() );
}

/// \brief Parse a segment from the specified segment string
region simple_chopping_format::parse_segment(const string_ref &arg_segment_string ///< The string from which to parse the segment
                                             ) const {
	constexpr char   CHAIN_OPEN_SQ_BR                = '[';
	constexpr char   CHAIN_CLOSE_SQ_BR               = ']';
	constexpr char   RESIDUE_NAME_DELIM              = '-';
	constexpr size_t MIN_VALID_CHARS                 = 6;
	constexpr size_t CHAIN_NEG_OFFSET                = 2;
	constexpr size_t CHAIN_OPEN_SQ_BR_END_NEG_OFFSET = 3;

	const auto length = arg_segment_string.length();
	if ( length < MIN_VALID_CHARS ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argh"));
	}
	if ( arg_segment_string.back() != CHAIN_CLOSE_SQ_BR || arg_segment_string[ length - CHAIN_OPEN_SQ_BR_END_NEG_OFFSET ] != CHAIN_OPEN_SQ_BR ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argh two"));
	}

	const auto begin_itr          = common::cbegin( arg_segment_string );
	const auto end_itr            = common::cend  ( arg_segment_string );
	const auto begin_plus_one_itr = next( begin_itr );
	const auto res_end_itr        = next( end_itr, 0 - static_cast<int>( CHAIN_OPEN_SQ_BR_END_NEG_OFFSET ) );
	const auto dash_itr           = find( begin_plus_one_itr, res_end_itr, RESIDUE_NAME_DELIM );

	if ( dash_itr == res_end_itr ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argh again"));
	}

	const auto dash_plus_one_itr = next( dash_itr );

	return {
		chain_label{ arg_segment_string[ length - CHAIN_NEG_OFFSET ] },
		parse_residue( make_string_ref( begin_itr,         dash_itr    ) ),
		parse_residue( make_string_ref( dash_plus_one_itr, res_end_itr ) )
	};
}

/// \brief Parse a residue from the specified residue string
residue_name simple_chopping_format::parse_residue(const string_ref &arg_string_ref ///< The string from which to parse the residue
                                                   ) const {
	constexpr char   INS_CODE_OPEN_BR             = '(';
	constexpr char   INS_CODE_CLOSE_BR            = ')';
	constexpr size_t MIN_INS_CODE_CHARS           = 4;
	constexpr size_t INS_CODE_OFFSET              = 2;
	constexpr size_t CHAIN_OPEN_BR_END_NEG_OFFSET = 3;

	const auto length = arg_string_ref.length();
	if ( length < MIN_INS_CODE_CHARS || arg_string_ref.back() != INS_CODE_CLOSE_BR ) {
		return residue_name{ stoi( arg_string_ref.to_string() ) };
	}

	if ( arg_string_ref[ length - CHAIN_OPEN_BR_END_NEG_OFFSET ] != INS_CODE_OPEN_BR ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argh yet again"));
	}
	const auto begin_itr   = common::cbegin( arg_string_ref );
	const auto end_itr     = common::cend  ( arg_string_ref );
	const auto res_end_itr = next( end_itr, 0 - static_cast<int>( CHAIN_OPEN_BR_END_NEG_OFFSET ) );

	return {
		stoi( string{ begin_itr, res_end_itr } ),
		arg_string_ref[ length - INS_CODE_OFFSET ]
	};
}