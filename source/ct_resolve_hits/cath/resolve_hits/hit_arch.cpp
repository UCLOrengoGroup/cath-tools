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

#include "hit_arch.hpp"

#include "cath/resolve_hits/full_hit.hpp"
#include "cath/resolve_hits/full_hit_list.hpp"
#include "cath/resolve_hits/full_hit_list_fns.hpp"
#include "cath/resolve_hits/options/spec/crh_segment_spec.hpp"

#include <string>
#include <type_traits>

using namespace ::cath::common;
using namespace ::cath::rslv;

using ::std::string;

// Come GCC >= 5.0, reinstate this static_assert that should be passing
//static_assert( std::is_nothrow_move_assignable_v   <hit_arch>, "" );
static_assert( std::is_nothrow_move_constructible_v<hit_arch> );

/// \brief Get a list of the full hits from the specified full_hit_list that correspond to
///        the hits in the specified hit_arch
///
/// \relates hit_arch
full_hit_list cath::rslv::get_full_hits_of_hit_arch(const hit_arch      &prm_hit_arch, ///< The hit_arch whose full hits should be extracted
                                                    const full_hit_list &prm_full_hits ///< The full_hit_list associated with the hit_arch, from which the full_hits should be extracted
                                                    ) {
	return full_hit_list{ transform_build<full_hit_vec>(
		prm_hit_arch,
		[&] (const calc_hit &x) { return prm_full_hits[ x.get_label_idx() ]; }
	) };
}

/// \brief Generate a string describing the specified hit_arch in the specified format
///
/// This is deliberately separated from a normal to_string() function because
/// the interface requirements of this may change (eg to demand that the client
/// passes the crh_score_spec, crh_segment_spec, hits_boundary_output)
///
/// \relates hit_arch
string cath::rslv::to_output_string(const hit_arch            &prm_hit_arch,         ///< The hit_arch to describe
                                    const full_hit_list       &prm_full_hits,        ///< The list of labels corresponding to the hit
                                    const crh_segment_spec    &prm_crh_segment_spec, ///< The crh_segment_spec specifying any trimming that should be performed on the output segments
                                    const hit_output_format   &prm_format,           ///< The format in which the hit_arch should be described
                                    const string              &prm_prefix,           ///< Any prefix that should come before the hit in hit_output_format::JON
                                    const hit_boundary_output &prm_boundary_output   ///< Whether to output the trimmed or original boundaries
                                    ) {
	return to_output_string(
		get_full_hits_of_hit_arch( prm_hit_arch, prm_full_hits ),
		prm_crh_segment_spec,
		prm_format,
		prm_prefix,
		prm_boundary_output
	);
}

