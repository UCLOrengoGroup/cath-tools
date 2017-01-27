/// \file
/// \brief The dna_atom class definitions

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

#include "amino_acid.hpp"

#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;

using std::string;

/// \brief Convert the specified dna_atom back to the three character string that represents it in PDBs
string cath::to_three_char_str(const dna_atom &arg_dna_atom ///< The dna_atom to convert
                               ) {
	switch ( arg_dna_atom ) {
		case ( dna_atom::DA ) : { return " DA"; }
		case ( dna_atom::DC ) : { return " DC"; }
		case ( dna_atom::DG ) : { return " DG"; }
		case ( dna_atom::DT ) : { return " DT"; }
		case ( dna_atom::DU ) : { return " DU"; }
		case ( dna_atom::A  ) : { return "  A"; }
		case ( dna_atom::C  ) : { return "  C"; }
		case ( dna_atom::G  ) : { return "  G"; }
		case ( dna_atom::N  ) : { return "  N"; }
		case ( dna_atom::U  ) : { return "  U"; }
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of dna_atom not recognised whilst converting to_three_char_str()"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
	return "DNA";
}