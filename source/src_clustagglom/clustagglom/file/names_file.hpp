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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_FILE_NAMES_FILE_H
#define _CATH_TOOLS_SOURCE_CLUSTER_FILE_NAMES_FILE_H

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.hpp"

#include <iosfwd>

namespace cath { namespace common { class id_of_str_bidirnl; } }

namespace cath {
	namespace clust {

		doub_vec parse_names(std::istream &,
		                     common::id_of_str_bidirnl &);

		doub_vec parse_names(const std::string &,
		                     common::id_of_str_bidirnl &);

		doub_vec parse_names(const boost::filesystem::path &,
		                     common::id_of_str_bidirnl &);

		std::pair<doub_vec, common::id_of_str_bidirnl> parse_names(std::istream &);

		std::pair<doub_vec, common::id_of_str_bidirnl> parse_names(const std::string &);

		std::pair<doub_vec, common::id_of_str_bidirnl> parse_names(const boost::filesystem::path &);


	} // namespace clust
} // namespace cath

#endif
