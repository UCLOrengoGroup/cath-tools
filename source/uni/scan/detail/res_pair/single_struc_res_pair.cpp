/// \file
/// \brief The single_struc_res_pair class definitions

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

#include "single_struc_res_pair.hpp"

using namespace cath::scan;
using namespace cath::scan::detail;

/// \brief TODOCUMENT
constexpr index_type single_struc_res_pair::DUMMY_INDEX_VALUE;

/// \brief TODOCUMENT
std::ostream & cath::scan::detail::operator<<(std::ostream                &arg_os,      ///< TODOCUMENT
                                              const single_struc_res_pair &arg_res_pair ///< TODOCUMENT
                                              ) {
	arg_os << "res_pair[from: ";
	arg_os << arg_res_pair.get_from_res_idx();
	arg_os << "; to: ";
	arg_os << arg_res_pair.get_to_res_idx();
	arg_os << "; ";
	arg_os << arg_res_pair.get_res_pair_core();
	arg_os << "]";
	return arg_os;
}
