/// \file
/// \brief The display_colour_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOUR_SPEC_H
#define _CATH_TOOLS_SOURCE_DISPLAY_DISPLAY_COLOURER_DISPLAY_COLOUR_SPEC_H

#include <boost/optional.hpp>

#include "display_colour/display_colour_type_aliases.h"
#include "display_colour/display_colour.h"

namespace cath { namespace align { class alignment; } }
namespace cath { namespace align { class alignment_context; } } 
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace sup { class superposition_context; } }
namespace cath { class viewer; }

namespace cath {

	/// \brief TODOCUMENT
	class display_colour_spec final {
	private:
		/// \brief TODOCUMENT
		display_colour_opt           base_clr;

		/// \brief TODOCUMENT
		size_display_colour_map      clr_of_pdb;

		/// \brief TODOCUMENT
		size_size_display_colour_map clr_of_pdb_and_res;

	public:
		void colour_base(const display_colour &,
		                 const bool &arg_overwrite = false);

		void colour_pdb(const size_t &,
		                const display_colour &,
		                const bool &arg_overwrite = false);

		void colour_pdb_residue(const size_t &,
		                        const size_t &,
		                        const display_colour &,
		                        const bool &arg_overwrite = false);

		const display_colour_opt & get_base_clr() const;

		const size_display_colour_map & get_clr_of_pdb() const;

		const size_size_display_colour_map & get_clr_of_pdb_and_res() const;
	};

	display_colour_opt get_clr_of_pdb_index(const display_colour_spec &,
	                                        const size_t &);

	display_colour_opt get_clr_of_pdb_and_res_indices(const display_colour_spec &,
	                                                  const size_t &,
	                                                  const size_t &);

	display_colour_vec get_pdb_colours(const display_colour_spec &);

	display_colour_vec get_residue_colours(const display_colour_spec &);

	display_colour_vec get_all_colours(const display_colour_spec &);

	size_vec get_pdbs_of_colour(const display_colour_spec &,
	                            const display_colour &);

	size_size_vec_map get_residues_of_colour(const display_colour_spec &,
	                                         const display_colour &);

	str_vec generate_colour_names(const size_t &);

	void colour_viewer_with_spec(const display_colour_spec &,
	                             const viewer &,
	                             const align::alignment_context &,
	                             std::ostream &);

//	void colour_alignment_with_spec(const display_colour_spec &,
//	                                const align::alignment &,
//	                                const pdb_list &,
//	                                const str_vec &,
//	                                std::ostream &);
} // namespace cath

#endif
