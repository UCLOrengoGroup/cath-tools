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

#ifndef _CATH_TOOLS_SOURCE_FILE_SEC_SEC_FILE_RECORD_H
#define _CATH_TOOLS_SOURCE_FILE_SEC_SEC_FILE_RECORD_H

#include "structure/geometry/coord.h"
#include "structure/geometry/rotation.h"
#include "structure/protein/sec_struc_type.h"

namespace cath { class sec_struc; }

namespace cath {
	namespace file {

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

			void check_unit_dirn_length();
			void check_sec_struc_type() const;

		public:
			sec_file_record(const size_t &,
			                const size_t &,
			                const sec_struc_type &,
			                const geom::coord &,
			                const geom::coord &);

			size_t         get_start_residue_num() const;
			size_t         get_stop_residue_num() const;
			sec_struc_type get_type() const;
			geom::coord    get_midpoint() const;
			geom::coord    get_unit_dirn() const;
		};

		sec_struc make_sec_struc(const sec_file_record &);

	} // namespace file
} // namespace cath
#endif
