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

#ifndef _CATH_TOOLS_SOURCE_UNI_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_CONSECUTIVE_HPP
#define _CATH_TOOLS_SOURCE_UNI_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_CONSECUTIVE_HPP

#include "cath/display/display_colourer/alignment_free_display_colourer.hpp"
#include "cath/display/display_colourer/detail/score_colour_handler.hpp"
#include "cath/display_colour/display_colour_list.hpp"

namespace cath {

	/// \brief TODOCUMENT
	///
	/// \todo It would be more flexible gradient colouring algorithms to be used for opacities and vice versa
	///       by having them both be stored in classes that generate a "display_value_spec", a sister to display_colour_spec
	///       that stores values rather than colours. This can then be translated to colours via a gradient or translated
	///       to opacities and applied to a display_colour_spec.
	class display_colourer_consecutive final : public alignment_free_display_colourer {
	private:
		/// \brief A useful type alias for the parent class
		using super = alignment_free_display_colourer;

		/// \brief The list of colours with which the structures should be coloured
		display_colour_list colours;

		std::unique_ptr<display_colourer> do_clone() const final;

		broad_display_colour_spec do_get_colour_spec_from_regions(const chop::region_vec_opt_vec &) const final;

		std::string do_get_label() const final;

		const display_colour_list & get_colours() const;

	public:
		explicit display_colourer_consecutive(display_colour_list);
		display_colourer_consecutive(display_colour_list,
		                             const detail::score_colour_handler &);
	};

} // namespace cath

#endif
