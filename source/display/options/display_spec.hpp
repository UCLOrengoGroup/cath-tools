/// \file
/// \brief The display_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_OPTIONS_DISPLAY_SPEC_H
#define _CATH_TOOLS_SOURCE_DISPLAY_OPTIONS_DISPLAY_SPEC_H

#include "common/type_aliases.hpp"

#include <memory>
#include <string>

namespace cath { class display_colour_list; }

namespace cath {

	/// \brief TODOCUMENT
	class display_spec final {
	private:
		static const std::string COLOURS_UNSPECIFIED;

		/// \brief A string describing the colours to use
		std::string display_colours_string = COLOURS_UNSPECIFIED;

		/// \brief Whether to display a gradient of colours
		bool gradient_colour_alignment = DEFAULT_GRADIENT_COLOUR_ALIGNMENT;

		/// \brief Whether to use colour to indicate scores (if they're present)
		bool show_scores_if_present    = DEFAULT_SHOW_SCORES_IF_PRESENT;

		/// \brief Whether to colour based on scores to the *present* equivalent positions
		bool scores_to_equivs          = DEFAULT_SCORES_TO_EQUIVS;

		/// \brief Whether to colour based on scores normalised across the alignment, rather than absolute scores
		bool normalise_scores          = DEFAULT_NORMALISE_SCORES;

	public:
		/// \brief The default value for whether to display a gradient of colours
		static constexpr bool DEFAULT_GRADIENT_COLOUR_ALIGNMENT = false;

		/// \brief The default value for whether to use colour to indicate scores (if they're present)
		static constexpr bool DEFAULT_SHOW_SCORES_IF_PRESENT    = false;

		/// \brief The default value for whether to colour based on scores to the *present* equivalent positions
		static constexpr bool DEFAULT_SCORES_TO_EQUIVS          = false;

		/// \brief The default value for whether to colour based on scores normalised across the alignment, rather than absolute scores
		static constexpr bool DEFAULT_NORMALISE_SCORES          = false;

		display_spec() = default;
		display_spec(const std::string &,
		             const bool &,
		             const bool &,
		             const bool &,
		             const bool &);

		str_opt get_display_colours_string() const;

		bool get_gradient_colour_alignment() const;
		bool get_show_scores_if_present() const;
		bool get_scores_to_equivs() const;
		bool get_normalise_scores() const;

		void set_display_colours_string(const std::string &);
		void set_gradient_colour_alignment(const bool &);
		void set_show_scores_if_present(const bool &);
		void set_scores_to_equivs(const bool &);
		void set_normalise_scores(const bool &);
	};

	bool requires_alignment(const display_spec &);
	bool is_consecutive(const display_spec &);

	bool has_display_colours_string(const display_spec &);
	str_opt invalid_string(const display_spec &);
	display_colour_list get_colour_list(const display_spec &);
} // namespace cath
#endif
