/// \file
/// \brief The dssp_file_io header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_DSSP_FILE_IO_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_DSSP_FILE_IO_HPP

#include <filesystem>
#include <utility>

// clang-format off
namespace cath { class chain_label; }
namespace cath { class dssp_struc_planar_angles; }
namespace cath { class residue; }
namespace cath::file { class dssp_file; }
namespace cath::file { class dssp_file_record; }
// clang-format on

namespace cath::file {

	dssp_file read_dssp_file(const ::std::filesystem::path &);
	dssp_file read_dssp(std::istream &);

	using size_residue_pair = std::pair<size_t, residue>;
	size_residue_pair parse_dssp_residue_line(const std::string &);

} // namespace cath::file

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_DSSP_FILE_IO_HPP
