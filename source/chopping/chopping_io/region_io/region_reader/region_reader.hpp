/// \file
/// \brief The region_reader class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER_REGION_READER_H
#define _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER_REGION_READER_H

#include <string>

namespace cath { namespace chop { class region; } }

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class region_reader {
		private:
			virtual region do_read_region(const std::string &) const = 0;

		public:
			region_reader() = default;
			virtual ~region_reader() noexcept = default;

			region_reader(const region_reader &) = default;
			region_reader(region_reader &&) noexcept = default;
			region_reader & operator=(const region_reader &) = default;
			region_reader & operator=(region_reader &&) noexcept = default;

			region read_region(const std::string &) const;
		};

	} // namespace chop
} // namespace cath

#endif
