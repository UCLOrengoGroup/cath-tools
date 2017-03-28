/// \file
/// \brief The alnd_rgn class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_ALND_RGN_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_ALND_RGN_H

#include "resolve_hits/res_arrow.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath {
	namespace rslv {

		/// \brief Represent a continuous region of (sequence) residues that are aligned to each other
		struct alnd_rgn {
			/// \brief The start of the aligned region in the first  sequence
			res_arrow start_res_a;

			/// \brief The start of the aligned region in the second sequence
			res_arrow start_res_b;

			/// \brief The length of the aligned region
			residx_t length;

		public:
			alnd_rgn(res_arrow,
			         res_arrow,
			         const residx_t &) noexcept;

			const res_arrow & get_start_res_a() const;
			const res_arrow & get_start_res_b() const;
			const residx_t & get_length() const;
		};

		std::string to_string(const alnd_rgn &);
		std::string to_string(const alnd_rgn_vec &);

	} // namespace rslv
} // namespace cath

#endif
