/// \file
/// \brief The quad_criteria_result class definitions

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

#include "quad_criteria_result.hpp"

#include "exception/invalid_argument_exception.hpp"

using namespace cath::common;
using namespace cath::scan::detail;
using namespace std;

/// \brief TODOCUMENT
///
/// \relates quad_criteria_result
string cath::scan::detail::to_string(const quad_criteria_result &arg_crit_res ///< TODOCUMENT
                                     ) {
	switch ( arg_crit_res ) {
		case ( quad_criteria_result::PASS                      ) : { return { "PASS"                      }; }
		case ( quad_criteria_result::QUERY_FAILS_SINGLE_CHECKS ) : { return { "QUERY_FAILS_SINGLE_CHECKS" }; }
		case ( quad_criteria_result::INDEX_FAILS_SINGLE_CHECKS ) : { return { "INDEX_FAILS_SINGLE_CHECKS" }; }
		case ( quad_criteria_result::FAILS_VIEW_CHECK          ) : { return { "FAILS_VIEW_CHECK"          }; }
		case ( quad_criteria_result::FAILS_PHI_CHECK           ) : { return { "FAILS_PHI_CHECK"           }; }
		case ( quad_criteria_result::FAILS_PSI_CHECK           ) : { return { "FAILS_PSI_CHECK"           }; }
		case ( quad_criteria_result::FAILS_FRAME_CHECK         ) : { return { "FAILS_FRAME_CHECK"         }; }
		case ( quad_criteria_result::FAILS_QUAD_CHECKS         ) : { return { "FAILS_QUAD_CHECKS"         }; }
		case ( quad_criteria_result::HAS_NO_REP                ) : { return { "HAS_NO_REP"                }; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of quad_criteria_result not recognised whilst converting to_string()"));
}

/// \brief TODOCUMENT
///
/// \relates quad_criteria_result
ostream & cath::scan::detail::operator<<(ostream                    &arg_os,      ///< TODOCUMENT
                                         const quad_criteria_result &arg_crit_res ///< TODOCUMENT
                                         ) {
	arg_os << to_string( arg_crit_res );
	return arg_os;
}
