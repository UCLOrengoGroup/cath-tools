/// \file
/// \brief The display_colourer_consecutive class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_CONSECUTIVE_H
#define _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_CONSECUTIVE_H

#include "display/display_colourer/detail/score_colour_handler.h"
#include "display/display_colourer/score_adjusted_display_colourer.h"
#include "display_colour/display_colour_list.h"

namespace cath {

	/// \brief TODOCUMENT
	///
	/// \todo Does display_colourer_alignment need to use score_adjusted_display_colourer through inheritance rather than aggregation?
	///       If this is developed further, it might be better to make score_adjusted_display_colourer a utility class.
	///
	/// \todo It would be more flexible gradient colouring algorithms to be used for opacities and vice versa
	///       by having them both be stored in classes that generate a "display_value_spec", a sister to display_colour_spec
	///       that stores values rather than colours. This can then be translated to colours via a gradient or translated
	///       to opacities and applied to a display_colour_spec.
	class display_colourer_consecutive final : public score_adjusted_display_colourer {
	private:
		/// \brief A useful type alias for the parent class
		using super = score_adjusted_display_colourer;

		/// \brief The list of colours with which the structures should be coloured
		display_colour_list colours;

		virtual std::unique_ptr<display_colourer> do_clone() const override final;

		virtual display_colour_spec do_get_colour_spec_before_scoring(const align::alignment_context &) const override final;

		const display_colour_list & get_colours() const;

	public:
		display_colourer_consecutive(const detail::score_colour_handler &,
		                             const display_colour_list &);
		virtual ~display_colourer_consecutive() noexcept = default;
	};

} // namespace cath

#endif
