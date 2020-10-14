/// \file
/// \brief The quad_and_rep_criteria_result class definitions

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

#include "quad_and_rep_criteria_result.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/less_than_helper.hpp"

using namespace ::cath::common;
using namespace ::cath::scan::detail;
using namespace ::std;

/// \brief TODOCUMENT
void quad_and_rep_criteria_result::sanity_check() const {
	if ( get_quad_status() == quad_criteria_result::HAS_NO_REP ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to construct quad_and_rep_criteria_result with quad_status of quad_criteria_result::HAS_NO_REP (which can only apply rep_status)"));
	}
}

/// \brief TODOCUMENT
quad_and_rep_criteria_result::quad_and_rep_criteria_result(const quad_criteria_result &prm_rep_status, ///< TODOCUMENT
                                                           const quad_criteria_result &prm_quad_status ///< TODOCUMENT
                                                           ) : rep_status  ( prm_rep_status  ),
                                                               quad_status ( prm_quad_status ) {
	sanity_check();
}

/// \brief TODOCUMENT
const quad_criteria_result & quad_and_rep_criteria_result::get_rep_status() const {
	return rep_status;
}


/// \brief TODOCUMENT
const quad_criteria_result & quad_and_rep_criteria_result::get_quad_status() const {
	return quad_status;
}

/// \brief TODOCUMENT
///
/// \relates quad_criteria_result
string cath::scan::detail::to_string(const quad_and_rep_criteria_result &prm_crit_res ///< TODOCUMENT
                                     ) {
	return "quad_criteria_result[rep:"
	       + to_string( prm_crit_res.get_rep_status() )
	       + "; quad:"
	       + to_string( prm_crit_res.get_quad_status() )
	       + "]";
}

/// \brief TODOCUMENT
///
/// \relates quad_criteria_result
ostream & cath::scan::detail::operator<<(ostream                            &prm_os,      ///< TODOCUMENT
                                         const quad_and_rep_criteria_result &prm_crit_res ///< TODOCUMENT
                                         ) {
	prm_os << to_string( prm_crit_res );
	return prm_os;
}

/// \brief TODOCUMENT
///
/// \relates quad_criteria_result
bool cath::scan::detail::operator<(const quad_and_rep_criteria_result &prm_criteria_result_a, ///< TODOCUMENT
                                   const quad_and_rep_criteria_result &prm_criteria_result_b  ///< TODOCUMENT
                                   ) {
	auto the_helper = make_less_than_helper( prm_criteria_result_a, prm_criteria_result_b );
	the_helper.register_comparison_field( [] (const quad_and_rep_criteria_result &x) { return x.get_rep_status();  } );
	the_helper.register_comparison_field( [] (const quad_and_rep_criteria_result &x) { return x.get_quad_status(); } );
	return final_less_than_result( the_helper );
}

/// \brief TODOCUMENT
///
/// \relates quad_criteria_result
bool cath::scan::detail::fully_passes(const quad_and_rep_criteria_result &prm_quad_and_rep_criteria_result ///< TODOCUMENT
                                      ) {
	return (
		quad_criteria_result::PASS == prm_quad_and_rep_criteria_result.get_rep_status()
		&&
		quad_criteria_result::PASS == prm_quad_and_rep_criteria_result.get_quad_status()
	);
}
