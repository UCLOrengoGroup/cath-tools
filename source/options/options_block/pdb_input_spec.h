/// \file
/// \brief The pdb_input_spec class header

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

#ifndef PDB_INPUT_SPEC_H_INCLUDED
#define PDB_INPUT_SPEC_H_INCLUDED

#include "common/type_aliases.h"

namespace cath {
	namespace opts {

		/// \brief Represent a specification for how PDBs should be read in
		class pdb_input_spec final {
		private:
			/// \brief A list of PDB files that should be read
			path_vec input_files;

			/// \brief Whether to read PDBs from stdin
			bool read_from_stdin = DEFAULT_READ_FROM_STDIN;

		public:
			/// \brief A default value for whether to read PDBs from stdin
			static constexpr bool DEFAULT_READ_FROM_STDIN = false;

			const path_vec & get_input_files() const;
			const bool & get_read_from_stdin() const;

			pdb_input_spec & set_input_files(const path_vec &);
			pdb_input_spec & set_read_from_stdin(const bool &);
		};

		size_t get_num_acquirers(const pdb_input_spec &);

		std::string to_string(const pdb_input_spec &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_input_spec &);

	}
}

#endif
