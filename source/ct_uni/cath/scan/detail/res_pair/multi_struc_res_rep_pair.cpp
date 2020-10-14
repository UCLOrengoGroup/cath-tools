/// \file
/// \brief The multi_struc_res_rep_pair class definitions

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

#include "multi_struc_res_rep_pair.hpp"

using namespace cath::scan::detail;
using namespace std;

/// \brief TODOCUMENT
ostream & cath::scan::detail::operator<<(ostream                &prm_os,      ///< TODOCUMENT
                                         const multi_struc_res_rep_pair &prm_res_pair ///< TODOCUMENT
                                         ) {
	prm_os << "res_pair[struc: ";
	prm_os << prm_res_pair.get_structure_index();
	prm_os << "; from_rep: ";
	prm_os << prm_res_pair.get_from_res_rep_index();
	prm_os << "; to_rep: ";
	prm_os << prm_res_pair.get_to_res_rep_index();
	prm_os << "; ";
	prm_os << prm_res_pair.get_res_pair_core();
	prm_os << "]";
	return prm_os;
}
