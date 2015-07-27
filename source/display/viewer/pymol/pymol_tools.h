/// \file
/// \brief The pymol_tools class header

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

#ifndef PYMOL_TOOLS_H_INCLUDED
#define PYMOL_TOOLS_H_INCLUDED

#include <cstddef>

namespace cath {

	/// \brief TODOCUMENT
	class pymol_tools final {
	private:
		pymol_tools() = delete;
		
	public:
		static double pymol_size(const size_t &,
		                         const double &,
		                         const size_t &,
		                         const double &,
		                         const size_t &);

	};

}

#endif
