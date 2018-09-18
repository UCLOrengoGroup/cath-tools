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

#include "res_pair_core.hpp"

using namespace cath::scan::detail;

using std::ostream;

 /// \brief TODOCUMENT
ostream & cath::scan::detail::operator<<(ostream             &prm_os,           ///< TODOCUMENT
                                         const res_pair_core &prm_res_pair_core ///< TODOCUMENT
                                         ) {
	prm_os << "";
	prm_os << "view: [";
	prm_os << prm_res_pair_core.get_view().get<0>();
	prm_os << ", ";
	prm_os << prm_res_pair_core.get_view().get<1>();
	prm_os << ", ";
	prm_os << prm_res_pair_core.get_view().get<2>();
	prm_os << "]; frame: ";
	prm_os << prm_res_pair_core.get_frame();
	prm_os << "; from_phi: ";
	prm_os << prm_res_pair_core.get_from_phi_angle();
	prm_os << "; from_psi: ";
	prm_os << prm_res_pair_core.get_from_psi_angle();
	prm_os << "; to_phi: ";
	prm_os << prm_res_pair_core.get_to_phi_angle();
	prm_os << "; to_psi: ";
	prm_os << prm_res_pair_core.get_to_psi_angle();
	return prm_os;
}
