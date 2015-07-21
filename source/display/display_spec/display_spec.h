/// \file
/// \brief The display_spec class header

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

#ifndef DISPLAY_SPEC_H_INCLUDED
#define DISPLAY_SPEC_H_INCLUDED

#include "display/display_spec/display_spec.h"
#include "common/type_aliases.h"
#include "display/display_colour/display_colour_gradient.h"

#include <memory>
#include <string>

namespace cath { class display_colour_list; }
namespace cath { class display_colourer; }
namespace cath { namespace opts { class display_options_block; } }
namespace display_spec_test_suite { struct basic; }

namespace cath {

	/// \brief TODOCUMENT
	///
	/// This makes display_options_block into a friend and is a specification that
	/// gives display_options_block special access to configure it.
	class display_spec final {
	private:
		friend class opts::display_options_block;
		friend struct display_spec_test_suite::basic;

		static const std::string COLOURS_UNSPECIFIED;

		/// \brief TODOCUMENT
		std::string display_colours_string = COLOURS_UNSPECIFIED;

		/// \brief TODOCUMENT
		bool gradient_colour_alignment = false;

		/// \brief TODOCUMENT
		bool show_scores_if_present = false;

		/// \brief TODOCUMENT
		bool scores_to_equivs = false;

		/// \brief TODOCUMENT
		bool normalise_scores = false;

		display_spec() = default;

		const std::string & get_display_colours_string_ref() const;
		std::string & get_display_colours_string_ref();
		bool & get_gradient_colour_alignment_ref();
		bool & get_show_scores_if_present_ref();
		bool & get_scores_to_equivs_ref();
		bool & get_normalise_scores_ref();

		bool has_display_colours_string() const;
		const std::string & get_display_colours_string_or_default() const;

		opt_str invalid_string() const;

		display_colour_list get_colour_list() const;
		bool get_gradient_colour_alignment() const;
		bool get_show_scores_if_present() const;
		bool get_scores_to_equivs() const;
		bool get_normalise_scores() const;

	public:
		display_spec(const std::string &,
		             const bool &,
		             const bool &,
		             const bool &,
		             const bool &);

		std::unique_ptr<const display_colourer> get_display_colourer(const display_colour_gradient &arg_colour_gradient = make_default_colour_gradient() ) const;
	};

}
#endif
