/// \file
/// \brief The broad_display_colour_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOUR_SPEC_BROAD_DISPLAY_COLOUR_SPEC_H
#define _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOUR_SPEC_BROAD_DISPLAY_COLOUR_SPEC_H

#include <boost/optional.hpp>

#include "display_colour/display_colour_type_aliases.h"
#include "display_colour/display_colour.h"

namespace cath { class viewer; }

namespace cath {

	/// \brief Represent the colouring of structures at a broad (ie not residue-specific) level
	class broad_display_colour_spec final {
	private:
		/// \brief The main overall colour to use
		display_colour_opt           base_clr;

		/// \brief The colours to use for the individual PDBs
		size_display_colour_map      clr_of_pdb;

	public:
		void colour_base(const display_colour &,
		                 const bool & = false);

		void colour_pdb(const size_t &,
		                const display_colour &,
		                const bool & = false);

		const display_colour_opt & get_base_clr() const;

		const size_display_colour_map & get_clr_of_pdb() const;
	};

	display_colour_opt get_clr_of_pdb_index(const broad_display_colour_spec &,
	                                        const size_t &);

	display_colour_vec get_pdb_colours(const broad_display_colour_spec &);

	size_vec get_pdbs_of_colour(const broad_display_colour_spec &,
	                            const display_colour &);

	str_vec generate_colour_names(const size_t &);

	namespace detail {
		void colour_pdbs_impl(const display_colour_vec &,
		                      const broad_display_colour_spec &,
		                      const viewer &,
		                      const str_vec &,
		                      std::ostream &);
	} // namespace detail

	void colour_viewer_with_spec(const broad_display_colour_spec &,
	                             const viewer &,
	                             const str_vec &,
	                             std::ostream &);

} // namespace cath

#endif
