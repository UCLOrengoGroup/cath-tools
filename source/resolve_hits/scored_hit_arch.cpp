/// \file
/// \brief The scored_hit_arch class definitions

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

#include "scored_hit_arch.h"

// #include <boost/algorithm/string/join.hpp>

#include "common/algorithm/transform_build.h"
#include "resolve_hits/hit_list.h"
// #include "common/boost_addenda/range/adaptor/lexical_casted.h"

// #include <string>

using namespace cath::common;
using namespace cath::rslv;

// using boost::algorithm::join;
// using std::ostream;
// using std::string;

// /// \brief TODOCUMENT
// ///
// /// \relates scored_hit_arch
// string cath::rslv::to_string(const scored_hit_arch &arg_scored_hit_arch ///< TODOCUMENT
//                              ) {
// 	return "scored_hit_arch[\n\t"
// 		+ join(
// 			arg_scored_hit_arch | lexical_casted<string>(),
// 			"\n\t"
// 		)
// 		+ "\n]";
// }

// /// \brief TODOCUMENT
// ///
// /// \relates scored_hit_arch
// ostream & cath::rslv::operator<<(ostream   &arg_ostream,      ///< TODOCUMENT
//                                  const scored_hit_arch &arg_scored_hit_arch ///< TODOCUMENT
//                                  ) {
// 	arg_ostream << to_string( arg_scored_hit_arch );
// 	return arg_ostream;
// }

/// \brief TODOCUMENT
///
/// \relates scored_hit_arch
scored_hit_arch cath::rslv::make_scored_hit_arch(const scored_arch_proxy &arg_scored_arch_proxy, ///< TODOCUMENT
                                                 const hit_list          &arg_hit_list           ///< TODOCUMENT
                                                 ) {
	return {
		arg_scored_arch_proxy.get_score(),
		hit_arch{
			transform_build<hit_vec>(
				arg_scored_arch_proxy,
				[&] (const hitidx_t &x) {
					return arg_hit_list[ x ];
				}
			)
		}
	};
}

