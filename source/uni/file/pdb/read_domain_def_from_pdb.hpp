/// \file
/// \brief The read_domain_def_from_pdb header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_PDB_READ_DOMAIN_DEF_FROM_PDB_HPP
#define _CATH_TOOLS_SOURCE_UNI_FILE_PDB_READ_DOMAIN_DEF_FROM_PDB_HPP

#include <boost/filesystem/path.hpp>

#include "chopping/domain/domain.hpp"

namespace cath { namespace chop { class domain_definition; } }
namespace cath { namespace opts { class data_dirs_spec; } }
namespace cath { namespace file { class pdb; } }

namespace cath {
	namespace file {

		std::pair<boost::filesystem::path, pdb> read_domain_from_pdb(const chop::domain_definition &,
		                                                             const opts::data_dirs_spec &);

		pdb read_domain_from_pdb_file(const boost::filesystem::path &,
		                              const chop::domain &);

	} // namespace file
} // namespace cath

#endif
