/// \file
/// \brief The rasmol_style_viewer class definitions

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

#include "rasmol_style_viewer.hpp"

using namespace ::cath;
using namespace ::cath::sup;
using namespace ::std;

/// \brief TODOCUMENT
string rasmol_style_viewer::do_default_file_extension() const {
	return ".spt";
}

/// \brief TODOCUMENT
void rasmol_style_viewer::do_write_alignment_extras(ostream                     &/*prm_os*/,                   ///< TODOCUMENT
                                                    const superposition_context &/*prm_superposition_context*/ ///< TODOCUMENT
                                                    ) const {
}
