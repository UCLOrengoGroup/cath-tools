/// \file
/// \brief The std_region_io_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_IO_REGION_IO_STD_REGION_IO_SPEC_H
#define _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_IO_REGION_IO_STD_REGION_IO_SPEC_H

#include "chopping/chopping_format/chopping_format.h"
#include "common/clone/clone_ptr.h"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class std_region_io_spec final {
		private:
			common::clone_ptr<chopping_format> chopping_format_ptr;

		public:
			explicit std_region_io_spec(const chopping_format &);
		};

	} // namespace chop
} // namespace cath

#endif
