/// \file
/// \brief The res_pair_index_dirn_criterion definitions

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

#include "res_pair_index_dirn_criterion.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"

#include <iostream>

using namespace ::cath::common;
using namespace ::cath::scan;
using namespace ::std;

/// \brief TODOCUMENT
///
/// \relates res_pair_index_dirn_criterion
ostream & cath::scan::operator<<(ostream                             &prm_os,                           ///< TODOCUMENT
                                 const res_pair_index_dirn_criterion &prm_res_pair_index_dirn_criterion ///< TODOCUMENT
                                 ) {
	switch ( prm_res_pair_index_dirn_criterion ) {
		case ( res_pair_index_dirn_criterion::MUST_MATCH     ) : { prm_os << "MUST_MATCH"     ; return prm_os ; }
		case ( res_pair_index_dirn_criterion::NEED_NOT_MATCH ) : { prm_os << "NEED_NOT_MATCH" ; return prm_os ; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of res_pair_index_dirn_criterion not recognised whilst inserting into an ostream"));
	
}
