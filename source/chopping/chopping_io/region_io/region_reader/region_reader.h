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

#ifndef REGION_READER_H_INCLUDED
#define REGION_READER_H_INCLUDED

#include <string>

namespace cath { namespace chop { class region; } }

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class region_reader {
		private:
			virtual region do_read_region(const std::string &) const = 0;

		public:
			virtual ~region_reader() noexcept = default;

			region read_region(const std::string &) const;
		};

	}
}

#endif
