/// \file
/// \brief The pdb_input_spec class definitions

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

#include "pdb_input_spec.hpp"

#include <boost/algorithm/string/join.hpp>

#include "cath/common/boost_addenda/range/adaptor/lexical_casted.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/string/booled_to_string.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::std::ostream;
using ::std::string;

constexpr bool pdb_input_spec::DEFAULT_READ_FROM_STDIN;

/// \brief Getter for the list of PDB files that should be read
const path_vec & pdb_input_spec::get_input_files() const {
	return input_files;
}

/// \brief Getter for whether to read PDBs from stdin
const bool & pdb_input_spec::get_read_from_stdin() const {
	return read_from_stdin;
}

/// \brief Setter for the list of PDB files that should be read
pdb_input_spec & pdb_input_spec::set_input_files(const path_vec &prm_input_files ///< The list of PDB files that should be read
                                                 ) {
	input_files = prm_input_files;
	return *this;
}

/// \brief Setter for whether to read PDBs from stdin
pdb_input_spec & pdb_input_spec::set_read_from_stdin(const bool &prm_read_from_stdin ///< Whether to read PDBs from stdin
                                                     ) {
	read_from_stdin = prm_read_from_stdin;
	return *this;
}

/// \brief Get the number of pdb_acquirer objects that would be created by get_pdbs_acquirers() on the specified pdb_input_spec
///
/// \relates pdb_input_spec
///
/// \alsorelates pdb_acquirer
size_t cath::opts::get_num_acquirers(const pdb_input_spec &prm_pdb_input_spec ///< The pdb_input_spec to query
                                     ) {
	return
		( ! prm_pdb_input_spec.get_input_files().empty() ? 1_z : 0_z )
		+
		(   prm_pdb_input_spec.get_read_from_stdin()     ? 1_z : 0_z );
}

/// \brief Generate a string describing the specified pdb_input_spec
///
/// \relates pdb_input_spec
string cath::opts::to_string(const pdb_input_spec &prm_pdb_input_spec ///< The pdb_input_spec to describe
                             ) {
	return "pdb_input_spec[input_files: "
		+ join( prm_pdb_input_spec.get_input_files() | lexical_casted<string>(), ", " )
		+ "; read_from_stdin: "
		+ booled_to_string( prm_pdb_input_spec.get_read_from_stdin() )
		+ "]";
}

/// \brief Insert a description of the specified pdb_input_spec into the specified ostream
///
/// \relates pdb_input_spec
ostream & cath::opts::operator<<(ostream              &prm_os,            ///< The ostream into which the description should be inserted
                                 const pdb_input_spec &prm_pdb_input_spec ///< The pdb_input_spec to describe
                                 ) {
	prm_os << to_string( prm_pdb_input_spec );
	return prm_os;
}
