/// \file
/// \brief The ostream_superposition_outputter class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef OSTREAM_SUPERPOSITION_OUTPUTTER_H_INCLUDED
#define OSTREAM_SUPERPOSITION_OUTPUTTER_H_INCLUDED

#include "options/outputter/superposition_outputter/superposition_outputter.h"

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class ostream_superposition_outputter final : public superposition_outputter {
		private:
			virtual std::unique_ptr<superposition_outputter> do_clone() const override final;
			virtual void do_output_superposition(const sup::superposition_context &,
			                                     std::ostream &) const override final;
			virtual bool do_involves_display_spec() const override final;

		public:
			virtual ~ostream_superposition_outputter() noexcept = default;
		};

	}
}

#endif
