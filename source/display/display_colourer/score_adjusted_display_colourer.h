/// \file
/// \brief The score_adjusted_display_colourer class header

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

#ifndef SCORE_ADJUSTED_DISPLAY_COLOURER_H_INCLUDED
#define SCORE_ADJUSTED_DISPLAY_COLOURER_H_INCLUDED

#include "display/display_colourer/detail/score_colour_handler.h"
#include "display/display_colourer/display_colourer.h"

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
	class score_adjusted_display_colourer : public display_colourer {
	private:
		/// \brief TODOCUMENT
		detail::score_colour_handler score_handler;

		/// \brief TODOCUMENT
		virtual display_colour_spec do_get_colour_spec_before_scoring(const align::alignment_context &) const = 0;

		virtual display_colour_spec do_get_colour_spec(const align::alignment_context &) const override final;

		const detail::score_colour_handler & get_score_handler() const;

	public:
		score_adjusted_display_colourer(const detail::score_colour_handler &);
		virtual ~score_adjusted_display_colourer() noexcept = default;
	};

}

#endif
