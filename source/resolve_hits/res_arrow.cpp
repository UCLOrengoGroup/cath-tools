/// \file
/// \brief The res_arrow class definitions

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

#include "resolve_hits/res_arrow.h"

#include <iostream>
#include <string>

using namespace cath::rslv;

using std::string;
using std::ostream;

static_assert( arrow_before_res( 5 ).res_before() == 4, "" );
static_assert( arrow_before_res( 5 ).res_after () == 5, "" );
static_assert( arrow_after_res ( 5 ).res_before() == 5, "" );
static_assert( arrow_after_res ( 5 ).res_after () == 6, "" );

static_assert( arrow_before_res( 0 ).res_after () == 0, "" );
static_assert( arrow_after_res ( 0 ).res_before() == 0, "" );
static_assert( arrow_after_res ( 0 ).res_after () == 1, "" );

static_assert( start_arrow().res_after() == 0, "" );

static_assert( arrow_after_res( 15 ) - arrow_after_res( 10 ) == 5, "" );

static_assert( arrow_after_res( 10 ) - 5 == arrow_after_res(  5 ), "" );
static_assert( arrow_after_res( 10 ) + 5 == arrow_after_res( 15 ), "" );

/// \brief Generate a string describing the specified res_arrow
///
/// \relates res_arrow
string cath::rslv::to_string(const res_arrow &arg_res_arrow ///< The res_arrow to describe
                             ) {
	return "res_arrow[before residue " + std::to_string( arg_res_arrow.res_after() ) + "]";
}

/// \brief Insert a description of the specified res_arrow into the specified ostream
///
/// \relates res_arrow
ostream & cath::rslv::operator<<(ostream         &arg_os,       ///< The ostream into which the description should be inserted
                                 const res_arrow &arg_res_arrow ///< The res_arrow to describe
                                 ) {
	arg_os << to_string( arg_res_arrow );
	return arg_os;
}
