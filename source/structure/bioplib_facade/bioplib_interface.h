/// \file
/// \brief The bioplib interface header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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
///
/// Currently, the bioplib subroutines are:
///  - ReadPDB
///  - WritePDB
///  - DupePDB
///  - FreePDB
///  - IndexPDB
///
///  - RenumAtomsPDB
///  - TranslatePDB
///  - ApplyMatrixPDB
///
///  - COOR
///  - PDB
///
///  - matfit()

#ifndef BIOPLIB_INTERFACE_H_INCLUDED
#define BIOPLIB_INTERFACE_H_INCLUDED

extern "C" {
	/// \todo This hits a problem when trying to compile against bioplib >= v3.0
	///       The included files (pdb.h -> general.h -> deprecated.h -> macros.h) define:
	/// ~~~~~
	/// extern int isatty(int);
	/// ~~~~~
	///       but this conflicts with the one with a `throw ()` specification
	///       that gets included into C++ code from unistd.h.
	///       This problem's rather tricky. The include of the `throw ()` version
	///       needn't happen in within the include tree of pdb.h for the problem to occur.
	///       For now, just use bioplib <= v2.1.2
	#include "bioplib/trunk/pdb.h"
}

#include <vector>

namespace cath { namespace file { class bioplib_pdb; } }
namespace cath { namespace geom { class coord; } }
namespace cath { namespace geom { class coord_list; } }
namespace cath { namespace geom { class rotation; } }

namespace cath {
	void translate_pdb(file::bioplib_pdb &,
	                   const geom::coord &);
	void rotate_pdb(file::bioplib_pdb &,
	                const geom::rotation &);

	geom::rotation bioplib_fit(const geom::coord_list &,
	                           const geom::coord_list &);

	/// \brief A class used in the bioplib_fit to temporarily make an array of COOR
	///        from a coord_list
	class bioplib_coord_list_interfacer final {
	private:
		std::vector<COOR> coors;
	public:
		explicit bioplib_coord_list_interfacer(const geom::coord_list &);

		COOR * get_interface();
	};

	COOR COOR_of_coord(const geom::coord &);
	geom::coord coord_of_COOR(const COOR &);
}

#endif
