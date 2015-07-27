/// \file
/// \brief The std_region_io_spec class definitions

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

#include "std_region_io_spec.h"

using namespace cath::chop;
//using namespace std;

/// \brief Ctor from a chopping_format
std_region_io_spec::std_region_io_spec(const chopping_format &arg_chopping_format ///< TODOCUMENT
                                       ) : chopping_format_ptr( arg_chopping_format.clone() ) {
}
