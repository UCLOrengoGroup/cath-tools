/// \file
/// \brief The std_region_writer class header

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

#ifndef STD_REGION_WRITER_H_INCLUDED
#define STD_REGION_WRITER_H_INCLUDED

#include "chopping/chopping_io/region_io/region_writer/region_writer.h"
#include "chopping/chopping_io/region_io/std_region_io_spec.h"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class std_region_writer final : public region_writer {
		private:
			/// \brief TODOCUMENT
			std_region_io_spec region_io_spec;

			virtual std::string do_write_region(const region &) const override final;

		public:
			explicit std_region_writer(const std_region_io_spec &);
		};

	}
}

#endif
