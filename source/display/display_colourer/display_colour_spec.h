/// \file
/// \brief The display_colour_spec class header

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

#ifndef DISPLAY_COLOUR_SPEC_H_INCLUDED
#define DISPLAY_COLOUR_SPEC_H_INCLUDED

#include <boost/optional.hpp>

#include "display/display_colour/display_colour.h"
#include "display/display_type_aliases.h"

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
		opt_display_colour           base_clr;

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

		const opt_display_colour & get_base_clr() const;

		const size_display_colour_map & get_clr_of_pdb() const;

		const size_size_display_colour_map & get_clr_of_pdb_and_res() const;
	};

	opt_display_colour get_clr_of_pdb_index(const display_colour_spec &,
	                                        const size_t &);

	opt_display_colour get_clr_of_pdb_and_res_indices(const display_colour_spec &,
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
}

#endif
