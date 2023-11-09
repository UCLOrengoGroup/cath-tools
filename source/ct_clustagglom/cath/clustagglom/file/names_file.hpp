/// \file
/// \brief The names_file class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE_NAMES_FILE_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE_NAMES_FILE_HPP

#include <filesystem>
#include <iosfwd>

#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::common { class id_of_str_bidirnl; }
// clang-format on

namespace cath::clust {

	doub_vec parse_names(std::istream &,
	                     common::id_of_str_bidirnl &);

	doub_vec parse_names(const std::string &,
	                     common::id_of_str_bidirnl &);

	doub_vec parse_names(const ::std::filesystem::path &,
	                     common::id_of_str_bidirnl &);

	std::pair<doub_vec, common::id_of_str_bidirnl> parse_names(std::istream &);

	std::pair<doub_vec, common::id_of_str_bidirnl> parse_names(const std::string &);

	std::pair<doub_vec, common::id_of_str_bidirnl> parse_names(const ::std::filesystem::path &);

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE_NAMES_FILE_HPP
