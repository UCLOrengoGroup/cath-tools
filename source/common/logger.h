/// \file
/// \brief The logger class header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef LOGGER_H_INCLUDED
#define LOGGER_H_INCLUDED

#include <string>

namespace cath {

	/// TODOCUMENT
	class logger final {
	public:
		logger() = delete;

		/// TODOCUMENT
		enum class return_code {
			SUCCESS                              =  0, ///< TODOCUMENT
			EXCEPTION_WITHOUT_SPECIFIC_RETCODE   =  1, ///< TODOCUMENT
			TOO_FEW_PDBS_FOR_ALIGNMENT           = 10, ///< TODOCUMENT
			TOO_MANY_PDBS_FOR_ALIGNMENT          = 20, ///< TODOCUMENT
			INSUFFICIENT_RESIDUE_NAME_OVERLAPS   = 30, ///< TODOCUMENT
			NO_SUCH_FILE                         = 40, ///< TODOCUMENT
			MALFORMED_PDB_FILE                   = 50, ///< TODOCUMENT
			UNABLE_TO_LOAD_SSAP_LEGACY_ALIGNMENT = 60, ///< TODOCUMENT
			NO_PDB_FILES_LOADED                  = 70, ///< TODOCUMENT
			MALFORMED_RESOLVE_HITS_INFILE        = 80  ///< TODOCUMENT
		};

		static void log_and_exit(const return_code &,
		                         const std::string &);
	};

}
#endif
