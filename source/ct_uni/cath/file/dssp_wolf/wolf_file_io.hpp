/// \file
/// \brief The wolf_file_io header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_WOLF_FILE_IO_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_WOLF_FILE_IO_HPP

#include <filesystem>
#include <vector>

// clang-format off
namespace cath::file { class wolf_file; }
// clang-format on

namespace cath::file {

	wolf_file read_wolf(const ::std::filesystem::path &);

} // namespace cath::file
#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_WOLF_FILE_IO_HPP
