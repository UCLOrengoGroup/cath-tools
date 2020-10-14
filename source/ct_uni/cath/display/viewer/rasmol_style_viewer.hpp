/// \file
/// \brief The rasmol_style_viewer class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_DISPLAY_VIEWER_RASMOL_STYLE_VIEWER_HPP
#define _CATH_TOOLS_SOURCE_UNI_DISPLAY_VIEWER_RASMOL_STYLE_VIEWER_HPP

#include "cath/display/viewer/viewer.hpp"

namespace cath {

	/// \brief TODOCUMENT
	class rasmol_style_viewer : public viewer {
	private:
		std::string do_default_executable() const override = 0;
		std::string do_default_file_extension() const final;
		void do_write_alignment_extras(std::ostream &,
		                               const sup::superposition_context &) const override = 0;
	};

} // namespace cath

#endif
