/// \file
/// \brief The display_colourer_score class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_SCORE_H
#define _CATH_TOOLS_SOURCE_UNI_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_SCORE_H

#include "display/display_colourer/detail/score_colour_handler.hpp"
#include "display/display_colourer/display_colourer.hpp"
#include "display_colour/display_colour_gradient.hpp"

namespace cath {

	/// \brief Concrete display_colourer that colours based on the strength of the score
	///        stretching it in proportion to the alignment's residue scores if they're present
	///
	/// If the alignment doesn't have residue_scores, then it is performed linearly.
	///
	/// The approach to making the progression through the gradient aims to:
	///  * use the two extremes of the gradient for the two ends of the alignment
	///  * be symmetric (ie reversing the alignment wouldn't affect the proportion of progress any two residues),
	///    which is achieved by making the amount of progress between two residues be the mean of their two scores
	///
	/// \todo It would be more flexible gradient colouring algorithms to be used for opacities and vice versa
	///       by having them both be stored in classes that generate a "display_value_spec", a sister to display_colour_spec
	///       that stores values rather than colours. This can then be translated to colours via a gradient or translated
	///       to opacities and applied to a display_colour_spec.
	class display_colourer_score final : public display_colourer {
	private:
		/// \brief A useful type alias for the parent class
		using super = display_colourer;

		/// \brief The gradient with which this display_colourer_score should colour alignments
		display_colour_gradient gradient;

		std::unique_ptr<display_colourer> do_clone() const final;

		display_colour_spec do_get_colour_spec(const align::alignment_context &) const final;

		std::string do_get_label() const final;

		const display_colour_gradient & get_gradient() const;

	public:
		explicit display_colourer_score(display_colour_gradient);
		display_colourer_score(display_colour_gradient,
		                       const detail::score_colour_handler &);
	};

} // namespace cath

#endif
