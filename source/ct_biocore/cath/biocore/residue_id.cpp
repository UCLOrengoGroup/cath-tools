/// \file
/// \brief The residue_id class definitions

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

#include <algorithm>
#include <iostream>
#include <iterator>
#include <optional>
#include <string>

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/throw_exception.hpp>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/biocore/residue_id.hpp"
#include "cath/biocore/residue_name.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/optional/make_optional_if.hpp"

using namespace ::cath;
using namespace ::cath::common;

using ::boost::algorithm::any_of;
using ::std::istream;
using ::std::nullopt;
using ::std::ostream;
using ::std::string;

/// \brief Generate a string describing the specified residue_id
///
/// \relates residue_id
string cath::to_string(const residue_id &prm_residue_id ///< The residue_id to describe
                       ) {
	return prm_residue_id.get_chain_label().to_string()
		+ ":"
		+ to_string( prm_residue_id.get_residue_name() );
}

/// \brief Insert a description of the specified residue_id into the specified ostream
///
/// \relates residue_id
ostream & cath::operator<<(ostream          &prm_os,        ///< The ostream into which the description should be inserted
                           const residue_id &prm_residue_id ///< The residue_id to describe
                           ) {
	prm_os << to_string( prm_residue_id );
	return prm_os;
}

/// \brief Extract into the specified residue_id from the specified stream
///
/// Expects the format to be a chain code character, followed by a colon, followed by a residue_name string
/// eg ("A:324A")
///
/// \relates residue_id
istream & cath::operator>>(istream    &prm_istream,   ///< The stream from which the residue_id should be extracted
                           residue_id &prm_residue_id ///< The residue_id to populate from the specified stream
                           ) {
	string input_string;
	prm_istream >> input_string;
	prm_residue_id = make_residue_id( input_string );
	return prm_istream;
}

/// \brief Generate a residue_id from the specified string
///
/// Expects the format to be a chain code character, followed by a colon, followed by a residue_name string
/// eg ("A:324A")
///
/// \relates residue_id
residue_id cath::make_residue_id(const string &prm_input_string ///< The string from which the residue_id should be parsed
                                 ) {
	if ( prm_input_string.empty() ) {
		return {};
	}
	if ( prm_input_string.length() < 3 || prm_input_string.at( 1 ) != ':' ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to parse residue ID from string"));
	}
	return {
		chain_label{ prm_input_string.front() },
		make_residue_name( prm_input_string.substr( 2 ) )
	};
}

/// \brief Return the specified residue_ids, grouped by chain_label
///
/// \relates residue_id
chain_label_residue_id_vec_map cath::get_residue_id_by_chain_label(const residue_id_vec &prm_residue_ids ///< The vector of residue_ids to query
                                                                   ) {
	chain_label_residue_id_vec_map results;
	for (const residue_id &res_id : prm_residue_ids) {
		results[ res_id.get_chain_label() ].emplace_back( res_id );
	}
	return results;
}

/// \brief Return a chain label that is used consistently in all of the specified residue_ids
///        or nullopt otherwise
///
/// Returns nullopt if `prm_residue_ids.empty()`
///
/// \relates residue_id
chain_label_opt cath::consistent_chain_label(const residue_id_vec &prm_residue_ids ///< The vector of residue_ids to query
                                             ) {
	if ( prm_residue_ids.empty() ) {
		return nullopt;
	}
	const auto &front_chain_label = prm_residue_ids.front().get_chain_label();
	return make_optional_if(
		all_of(
			next( cbegin( prm_residue_ids ) ),
			cend( prm_residue_ids ),
			[&] (const residue_id &x) { return x.get_chain_label() == front_chain_label; }
		),
		front_chain_label
	);
}

/// \brief Return whether the specified residue IDs have a consistent chain label
///
/// Returns true if `prm_residue_ids.empty()`
///
/// \relates residue_id
bool cath::have_consistent_chain_labels(const residue_id_vec &prm_residue_ids ///< The vector of residue_ids to query
                                        ) {
	return prm_residue_ids.empty() || consistent_chain_label( prm_residue_ids );
}

/// \brief Return whether any of the specified residue_ids has a strictly negative residue number
///
/// \relates residue_id
bool cath::has_any_strictly_negative_residue_numbers(const residue_id_vec &prm_residue_ids ///< The residue_ids to query
                                                     ) {
	return any_of(
		prm_residue_ids,
		[] (const residue_id &x) { return has_strictly_negative_residue_number( x ); }
	);
}
