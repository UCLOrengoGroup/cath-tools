/// \file
/// \brief The run_pymol definitions

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

#include "run_pymol.hpp"

#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>

#include "common/command_executer.hpp"

using namespace cath;

using boost::filesystem::path;

/// Run the specified pymol with the specified script file
///
/// \param prm_script_file    The script file to pass to pymol
/// \param prm_pymol_program  The pymol command to call (default: "pymol")
void cath::view::run_pymol( const path &prm_script_file, const path &prm_pymol_program ) {
	const bool pymol_success = command_executer::execute( prm_pymol_program, { prm_script_file.string() } );
	if ( !pymol_success ) {
		BOOST_LOG_TRIVIAL( warning ) << "PyMOL executable " + prm_pymol_program.string() + " did not run/shutdown normally.";
	}
}
