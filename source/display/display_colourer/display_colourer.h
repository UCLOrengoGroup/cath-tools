/// \file
/// \brief The display_colourer class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_H
#define _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOURER_H

#include <boost/optional.hpp>

#include "common/type_aliases.h"
#include "display/display_colourer/detail/score_colour_handler.h"
#include "display_colour/display_colour_gradient.h"

#include <iosfwd>
#include <memory>

namespace cath { class alignment_free_display_colourer; }
namespace cath { class display_colour_spec; }
namespace cath { class display_spec; }
namespace cath { class viewer; }
namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { class alignment_context; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {

	/// \brief TODOCUMENT
	class display_colourer {
	private:
		/// \brief Optional specification for post-modifying the colouring based on scores
		detail::score_colour_handler_opt the_score_colour_handler{ boost::none };

		/// \brief Pure virtual method with which each concrete display_colourer must define how to create a clone of itself
		virtual std::unique_ptr<display_colourer> do_clone() const = 0;

		virtual display_colour_spec do_get_colour_spec(const align::alignment_context &) const = 0;

	public:
		display_colourer() noexcept = default;
		explicit display_colourer(const detail::score_colour_handler &);
		virtual ~display_colourer() noexcept = default;

		std::unique_ptr<display_colourer> clone() const;

		const detail::score_colour_handler_opt & get_score_colour_handler_opt() const;
		display_colour_spec get_colour_spec(const align::alignment_context &) const;
	};

	bool has_score_colour_handler(const display_colourer &);
	const detail::score_colour_handler & get_score_colour_handler(const display_colourer &);

	std::unique_ptr<const display_colourer> get_display_colourer(const display_spec &,
	                                                             const display_colour_gradient & = make_default_colour_gradient() );

	display_colour_spec get_colour_spec(const display_colourer &,
	                                    const file::pdb_list &,
	                                    const str_vec &,
	                                    const align::alignment &);

	void colour_viewer(const display_colourer &,
	                   std::ostream &,
	                   const viewer &,
	                   const align::alignment_context &);

	void colour_viewer(const alignment_free_display_colourer &,
                       std::ostream &,
                       const viewer &,
                       const str_vec &);

} // namespace cath

#endif
