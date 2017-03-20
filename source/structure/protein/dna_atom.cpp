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

#include "common/string/char_arr_to_string.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;

using std::string;

/// \brief Convert the specified dna_atom back to the three character string that represents it in PDBs
char_3_arr cath::to_three_char_arr(const dna_atom &arg_dna_atom ///< The dna_atom to convert
                                   ) {
	switch ( arg_dna_atom ) {
		case (dna_atom::A       ) : { return {{ ' ',' ','A' }}; }
		case (dna_atom::C       ) : { return {{ ' ',' ','C' }}; }
		case (dna_atom::G       ) : { return {{ ' ',' ','G' }}; }
		case (dna_atom::I       ) : { return {{ ' ',' ','I' }}; } // eg in 1xnr
		case (dna_atom::N       ) : { return {{ ' ',' ','N' }}; }
		case (dna_atom::T       ) : { return {{ ' ',' ','T' }}; } // eg in 3dpv
		case (dna_atom::U       ) : { return {{ ' ',' ','U' }}; }
		case (dna_atom::DA      ) : { return {{ ' ','D','A' }}; }
		case (dna_atom::DC      ) : { return {{ ' ','D','C' }}; }
		case (dna_atom::DG      ) : { return {{ ' ','D','G' }}; }
		case (dna_atom::DI      ) : { return {{ ' ','D','I' }}; } // eg in 3rzl
		case (dna_atom::DT      ) : { return {{ ' ','D','T' }}; }
		case (dna_atom::DU      ) : { return {{ ' ','D','U' }}; }
		case (dna_atom::PLUS_A  ) : { return {{ ' ','+','A' }}; } // eg in 1hp6
		case (dna_atom::PLUS_C  ) : { return {{ ' ','+','C' }}; } // eg in 356d
		case (dna_atom::PLUS_G  ) : { return {{ ' ','+','G' }}; } // eg in 1gpg
		case (dna_atom::PLUS_U  ) : { return {{ ' ','+','U' }}; } // eg in 1hp6
		case (dna_atom::A_SPACE ) : { return {{ ' ','A',' ' }}; } // eg in 3zvp
		case (dna_atom::C_SPACE ) : { return {{ ' ','C',' ' }}; } // eg in 3zvp
		case (dna_atom::G_SPACE ) : { return {{ ' ','G',' ' }}; } // eg in 4a1d
		case (dna_atom::U_SPACE ) : { return {{ ' ','U',' ' }}; } // eg in 4a1d
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of dna_atom not recognised whilst converting to_three_char_str()"));
		}
	}
	return {{ 'D','N','A' }};
}

/// \brief Convert the specified dna_atom back to the three character string that represents it in PDBs
string cath::to_three_char_str(const dna_atom &arg_dna_atom ///< The dna_atom to convert
                               ) {
	return char_arr_to_string( to_three_char_arr( arg_dna_atom ) );
}