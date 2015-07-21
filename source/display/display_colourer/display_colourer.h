/// \file
/// \brief The display_colourer class header

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

#ifndef DISPLAY_COLOURER_H_INCLUDED
#define DISPLAY_COLOURER_H_INCLUDED

#include "common/type_aliases.h"

#include <iosfwd>
#include <memory>

namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { class alignment_context; } }
namespace cath { class display_colour_spec; }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace sup { class superposition_context; } }
namespace cath { class viewer; }

namespace cath {

	/// \brief TODOCUMENT
	class display_colourer {
	private:

		/// \brief Pure virtual method with which each concrete display_colourer must define how to create a clone of itself
		virtual std::unique_ptr<display_colourer> do_clone() const = 0;

		virtual display_colour_spec do_get_colour_spec(const align::alignment_context &) const = 0;

	public:
		virtual ~display_colourer() noexcept = default;

		std::unique_ptr<display_colourer> clone() const;

		display_colour_spec get_colour_spec(const align::alignment_context &) const;
	};

	display_colour_spec get_colour_spec(const display_colourer &,
	                                    const file::pdb_list &,
	                                    const str_vec &,
	                                    const align::alignment &);

	void colour_viewer(const display_colourer &,
	                   std::ostream &,
	                   const viewer &,
	                   const align::alignment_context &);

}

#endif
