/// \file
/// \brief The region_reader class definitions

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

#include "region_reader.hpp"

#include "cath/chopping/region/region.hpp"

using namespace ::cath::chop;

using ::std::string;

/// \brief TODOCUMENT
region region_reader::read_region(const string &prm_region_string ///< TODOCUMENT
                                  ) const {
	return do_read_region( prm_region_string );
}
