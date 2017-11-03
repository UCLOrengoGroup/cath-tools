/// \file
/// \brief The residue_locating definitions

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

#include "residue_locating.hpp"

#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"

using namespace cath::chop;
using namespace cath::common;

using std::ostream;
using std::string;

/// \brief Generate a string describing the specified residue_locating
///
/// \relates residue_locating
string cath::chop::to_string(const residue_locating &arg_residue_locating ///< The residue_locating to describe
                             ) {
	switch ( arg_residue_locating ) {
		case ( residue_locating::INDEX          ) : { return "res_loc_index"          ; }
		case ( residue_locating::NAME           ) : { return "res_loc_name"           ; }
		case ( residue_locating::NAME_AND_INDEX ) : { return "res_loc_name_and_index" ; }
	}
	BOOST_THROW_EXCEPTION(out_of_range_exception("residue_locating value not recognised"));
}

/// \brief Insert a description of the specified residue_locating into the specified ostream
///
/// \relates residue_locating
ostream & cath::chop::operator<<(ostream                &arg_os,              ///< The ostream into which the description should be inserted
                                 const residue_locating &arg_residue_locating ///< The residue_locating to describe
                                 ) {
	arg_os << to_string( arg_residue_locating );
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates residue_locating
residue_locating cath::chop::make_residue_locating_of_has_name_and_has_index(const bool &arg_has_name, ///< TODOCUMENT
                                                                             const bool &arg_has_index ///< TODOCUMENT
                                                                             ) {
	/// \todo: Come C++14, check out if this can be done with a switch on a using constexpr pair case labels
	///        (eg `switch ( make_pair( has_name, has_index ) ) { case( pair<bool>( true, true) ) [...]`)
	if (   arg_has_name &&   arg_has_index ) {
		return residue_locating::NAME_AND_INDEX;
	}
	if (   arg_has_name && ! arg_has_index ) {
		return residue_locating::NAME;
	}
	if ( ! arg_has_name &&   arg_has_index ) {
		return residue_locating::INDEX;
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot locate residues without using name or index"));
}
