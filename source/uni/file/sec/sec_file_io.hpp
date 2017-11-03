/// \file
/// \brief The sec_file_io header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_SEC_SEC_FILE_IO_H
#define _CATH_TOOLS_SOURCE_FILE_SEC_SEC_FILE_IO_H

#include <boost/filesystem/path.hpp>

#include "structure/structure_type_aliases.hpp"

#include <iosfwd>
#include <vector>

namespace cath { namespace file { class sec_file_record; } }
namespace cath { namespace file { class sec_file; } }
namespace cath { class sec_struc_planar_angles; }

namespace cath {
	namespace file {
		sec_file read_sec(const boost::filesystem::path &);
		sec_file read_sec(std::istream &);

		namespace detail {
			std::pair<size_t, sec_file_record> parse_sec_main_line(const std::string &);
		} // namespace detail

		sec_struc_planar_angles_vec parse_sec_angles_line(const std::string &);

	} // namespace file
} // namespace cath
#endif
