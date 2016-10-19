/// \file
/// \brief The plate_scan class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_DETAIL_PLATE_PLATE_SCAN_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_DETAIL_PLATE_PLATE_SCAN_H

#include <cstddef>

namespace cath {
	namespace index {
		namespace detail {

			/// \brief TODOCUMENT
			///
			/// \todo This should be an iterator
			class plate_scan final {
			private:
				/// \brief TODOCUMENT
				size_t from_a;

				/// \brief TODOCUMENT
				size_t to_a;

				/// \brief TODOCUMENT
				size_t size_a;


				/// \brief TODOCUMENT
				size_t from_b;

				/// \brief TODOCUMENT
				size_t to_b;

				/// \brief TODOCUMENT
				size_t size_b;


				/// \brief TODOCUMENT
				size_t from_step_size;

				/// \brief TODOCUMENT
				size_t to_step_size;

			public:
				bool operator()() const;

				// operator++();
				// operator*();
			};

		} // namespace detail
	} // namespace index
} // namespace cath

#endif
