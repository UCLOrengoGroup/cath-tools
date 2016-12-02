/// \file
/// \brief The std_region_reader class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER_STD_REGION_READER_H
#define _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_IO_REGION_IO_REGION_READER_STD_REGION_READER_H

#include "chopping/chopping_io/region_io/region_reader/region_reader.hpp"
#include "chopping/chopping_io/region_io/std_region_io_spec.hpp"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class std_region_reader final : public region_reader {
		private:
			/// \brief TODOCUMENT
			std_region_io_spec region_io_spec;

			virtual region do_read_region(const std::string &) const override final;

		public:
			explicit std_region_reader(const std_region_io_spec &);
		};

	} // namespace chop
} // namespace cath

#endif
