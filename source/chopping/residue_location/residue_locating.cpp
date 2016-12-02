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

#include "exception/invalid_argument_exception.hpp"

using namespace cath::chop;
using namespace cath::common;

/// \brief TODOCUMENT
residue_locating cath::chop::make_residue_locating_of_has_name_and_has_index(const bool &arg_has_name, ///< TODOCUMENT
                                                                             const bool &arg_has_index ///< TODOCUMENT
                                                                             ) {
	/// \todo: Come C++14, check out if this can be done with a switch on a using constexpr pair case labels
	///        (eg `switch ( make_pair( has_name, has_index ) ) { case( pair<bool>( true, true) ) [...]`)
	if      (   arg_has_name &&   arg_has_index ) {
		return residue_locating::NAME_AND_INDEX;
	}
	else if (   arg_has_name && ! arg_has_index ) {
		return residue_locating::NAME;
	}
	else if ( ! arg_has_name &&   arg_has_index ) {
		return residue_locating::INDEX;
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot locate residues without using name or index"));
	return residue_locating::NAME_AND_INDEX; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
}
