/// \file
/// \brief The sec_file class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_SEC_SEC_FILE_H
#define _CATH_TOOLS_SOURCE_FILE_SEC_SEC_FILE_H

#include "structure/structure_type_aliases.hpp"

#include <cstddef>
#include <vector>

namespace cath { namespace file { class sec_file_record; } }
namespace cath { class sec_struc; }
namespace cath { class sec_struc_planar_angles; }

namespace cath {
	namespace file {

		/// \brief Represent the data parsed out of a sec file
		class sec_file final {
		private:
			using record_list_type = std::vector<sec_file_record>;

			/// \brief The list of sec_file_records
			record_list_type records;

			/// \brief The list of planar_angle_lists
			///        (as in the sec file, each entry contains the angles between that secondary structure and
			///         all of the following secondary structures)
			sec_struc_planar_angles_vec_vec inter_planar_angles;

		public:
			sec_file(const std::vector<sec_file_record> &,
			         const sec_struc_planar_angles_vec_vec &);

			using iterator       = record_list_type::const_iterator;
			using const_iterator = record_list_type::const_iterator;
			using size_type      = record_list_type::size_type;

			size_type size() const;
			const sec_struc_planar_angles & get_planar_angles_of_indices(const size_t &,
			                                                             const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		sec_struc_vec make_sec_struc_list(const sec_file &);

	} // namespace file
} // namespace cath

#endif
