/// \file
/// \brief The res_pair_dirn class definitions

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

#include "res_pair_dirn.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath::common;
using namespace ::cath::scan::detail;
using namespace ::std;

/// \brief TODOCUMENT
ostream & cath::scan::detail::operator<<(ostream             &prm_os,           ///< TODOCUMENT
                                         const res_pair_dirn &prm_res_pair_dirn ///< TODOCUMENT
                                         ) {
	switch ( prm_res_pair_dirn ) {
		case ( res_pair_dirn::INCREASE ) : { prm_os << "res_pair_dirn::INCREASE" ; return prm_os ; }
		case ( res_pair_dirn::DECREASE ) : { prm_os << "res_pair_dirn::DECREASE" ; return prm_os ; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of res_pair_dirn not recognised whilst inserting into an ostream"));
}

