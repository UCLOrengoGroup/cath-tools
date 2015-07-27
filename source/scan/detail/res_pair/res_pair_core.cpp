/// \file
/// \brief The res_pair_core class definitions

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

#include "res_pair_core.h"

using namespace cath::scan::detail;
 using namespace std;

 /// \brief TODOCUMENT
ostream & cath::scan::detail::operator<<(ostream             &arg_os,           ///< TODOCUMENT
                                         const res_pair_core &arg_res_pair_core ///< TODOCUMENT
                                         ) {
	arg_os << "";
	arg_os << "view: [";
	arg_os << arg_res_pair_core.get_view().get<0>();
	arg_os << ", ";
	arg_os << arg_res_pair_core.get_view().get<1>();
	arg_os << ", ";
	arg_os << arg_res_pair_core.get_view().get<2>();
	arg_os << "]; frame: ";
	arg_os << arg_res_pair_core.get_frame();
	arg_os << "; from_phi: ";
	arg_os << arg_res_pair_core.get_from_phi_angle();
	arg_os << "; from_psi: ";
	arg_os << arg_res_pair_core.get_from_psi_angle();
	arg_os << "; to_phi: ";
	arg_os << arg_res_pair_core.get_to_phi_angle();
	arg_os << "; to_psi: ";
	arg_os << arg_res_pair_core.get_to_psi_angle();
	return arg_os;
}
