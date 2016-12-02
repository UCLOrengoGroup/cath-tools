/// \file
/// \brief The alignment_free_display_colourer class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_ALIGNMENT_FREE_DISPLAY_COLOURER_H
#define _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_ALIGNMENT_FREE_DISPLAY_COLOURER_H

#include "display/display_colour_spec/broad_display_colour_spec.hpp"
#include "display/display_colourer/display_colourer.hpp"

namespace cath {

	/// \brief ABC for display_colourers that can colour without knowledge of the alignment,
	///        only the number of entries
	class alignment_free_display_colourer : public display_colourer {
	private:
		using super = display_colourer;

		/// \brief Pure virtual method with which each concrete alignment_free_display_colourer must define how to colour based on the number of entries
		virtual broad_display_colour_spec do_get_colour_spec_from_num_entries(const size_t &) const = 0;

		virtual display_colour_spec do_get_colour_spec(const align::alignment_context &) const override final;

	public:
		alignment_free_display_colourer() noexcept = default;
		explicit alignment_free_display_colourer(const detail::score_colour_handler &);

		broad_display_colour_spec get_colour_spec_from_num_entries(const size_t &) const;
	};

} // namespace cath

#endif
