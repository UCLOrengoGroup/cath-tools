/// \file
/// \brief The sec_file_record class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_SEC_SEC_FILE_RECORD_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_SEC_SEC_FILE_RECORD_HPP

#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/rotation.hpp"
#include "cath/structure/protein/sec_struc_type.hpp"

// clang-format off
namespace cath { class sec_struc; }
// clang-format on

namespace cath::file {

	/// \brief Represent the data for a single secondary structure, as parsed out of the main line of a sec file
	class sec_file_record final {
	private:
		/// \brief The start residue of the secondary structure (\todo sequential number? offset?)
		size_t         start_residue_num;
		/// \brief The stop residue of the secondary structure (\todo sequential number? offset?)
		size_t         stop_residue_num;
		/// \brief The type of secondary structure: sec_struc_type::ALPHA_HELIX or sec_struc_type::BETA_STRAND (represented as 'H' or 'S' respectively in sec files)
		sec_struc_type type;
		/// \brief The coordinates of the secondary structure's midpoint
		geom::coord    midpoint;
		/// \brief A unit vector along the secondary structure (\todo towards the start of the secondary structure?)
		geom::coord    unit_dirn;

		void check_unit_dirn_length() const;
		void check_sec_struc_type() const;

	public:
		sec_file_record(const size_t &,
		                const size_t &,
		                const sec_struc_type &,
		                geom::coord,
		                geom::coord);

		[[nodiscard]] size_t         get_start_residue_num() const;
		[[nodiscard]] size_t         get_stop_residue_num() const;
		[[nodiscard]] sec_struc_type get_type() const;
		[[nodiscard]] geom::coord    get_midpoint() const;
		[[nodiscard]] geom::coord    get_unit_dirn() const;
	};

	bool operator==(const sec_file_record &,
	                const sec_file_record &);

	std::string to_string(const sec_file_record &);
	std::ostream & operator<<(std::ostream &,
	                          const sec_file_record &);

	sec_struc make_sec_struc(const sec_file_record &);

} // namespace cath::file
#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_SEC_SEC_FILE_RECORD_HPP
