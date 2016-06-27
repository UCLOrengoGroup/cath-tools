/// \file
/// \brief The hit_arch class definitions

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

#include "hit_arch.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/boost_addenda/range/adaptor/lexical_casted.h"

#include <string>
#include <type_traits>

using namespace cath::common;
using namespace cath::rslv;

using boost::adaptors::transformed;
using boost::algorithm::join;
using std::ostream;
using std::string;

// Come GCC >= 5.0, reinstate this static_assert that should be passing
//static_assert( std::is_nothrow_move_assignable   <hit_arch>::value, "" );
static_assert( std::is_nothrow_move_constructible<hit_arch>::value, "" );

/// \brief TODOCUMENT
///
/// \relates hit_arch
string cath::rslv::to_string(const hit_arch          &arg_hit_arch, ///< TODOCUMENT
                             const hit_output_format &arg_format    ///< TODOCUMENT
                             ) {
	const bool is_jon = ( arg_format == hit_output_format::JON );
	const string prefix    = is_jon ? ""   : "hit_arch[\n\t";
	const string separator = is_jon ? "\n" : "\n\t";
	const string suffix    = is_jon ? "\n" : "\n]";
	return prefix
		+ join(
			arg_hit_arch
				| transformed( [&] (const hit &x) { return to_string( x, arg_format ); } ),
			separator
		)
		+ suffix;
}

/// \brief TODOCUMENT
///
/// \relates hit_arch
ostream & cath::rslv::operator<<(ostream        &arg_ostream, ///< TODOCUMENT
                                 const hit_arch &arg_hit_arch ///< TODOCUMENT
                                 ) {
	arg_ostream << to_string( arg_hit_arch );
	return arg_ostream;
}
