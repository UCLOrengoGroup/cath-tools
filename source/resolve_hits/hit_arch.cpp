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
#include "resolve_hits/full_hit_list.h"

#include <string>
#include <type_traits>

using namespace cath::common;
using namespace cath::rslv;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::optional;
using std::ostream;
using std::string;

// Come GCC >= 5.0, reinstate this static_assert that should be passing
//static_assert( std::is_nothrow_move_assignable   <hit_arch>::value, "" );
static_assert( std::is_nothrow_move_constructible<hit_arch>::value, "" );

/// \brief Generate a string describing the specified hit_arch in the specified format
///
/// This is deliberately separated from a normal to_string() function because
/// the interface requirements of this may change (eg to demand that the client
/// passes the crh_score_spec, trim_spec, hits_boundary_output)
///
/// \relates hit_arch
string cath::rslv::to_output_string(const hit_arch            &arg_hit_arch,     ///< The hit_arch to describe
                                    const full_hit_list       &arg_full_hits,    ///< The list of labels corresponding to the hit
                                    const hit_output_format   &arg_format,       ///< The format in which the hit_arch should be described
                                    const string              &arg_prefix,       ///< Any prefix that should come before the hit in hit_output_format::JON
                                    const optional<trim_spec> &arg_trim_spec_opt ///< An optional trim_spec specifying any trimming that should be performed on the output segments
                                    ) {
	const bool is_jon = ( arg_format == hit_output_format::JON );
	const string prefix    = is_jon ? ""   : "hit_arch[\n\t";
	const string separator = is_jon ? "\n" : "\n\t";
	const string suffix    = is_jon ? "\n" : "\n]";
	return prefix
		+ join(
			arg_hit_arch
				| transformed( [&] (const calc_hit &x) {
					return to_string(
						arg_full_hits[ x.get_label_idx() ],
						arg_format,
						arg_prefix,
						arg_trim_spec_opt
					);
				} ),
			separator
		)
		+ suffix;
}

