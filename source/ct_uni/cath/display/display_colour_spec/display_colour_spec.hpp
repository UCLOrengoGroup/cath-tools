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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC_DISPLAY_COLOUR_SPEC_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC_DISPLAY_COLOUR_SPEC_HPP

#include "cath/display/display_colour_spec/broad_display_colour_spec.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"

// clang-format off
namespace cath::align { class alignment_context; }
// clang-format on

namespace cath {

	/// \brief Represent the colouring of structures, possibly down to the residue-specific level
	class display_colour_spec final {
	private:
		/// \brief The broad level colouring
		broad_display_colour_spec    the_broad_spec;

		/// \brief Any residue-specific colouring
		size_size_display_colour_map clr_of_pdb_and_res;

	public:
		display_colour_spec() = default;
		explicit display_colour_spec(broad_display_colour_spec);

		void colour_base(const display_colour &,
		                 const bool & = false);

		void colour_pdb(const size_t &,
		                const display_colour &,
		                const bool & = false);

		void colour_pdb_residue(const size_t &,
		                        const size_t &,
		                        const display_colour &,
		                        const bool & = false);

		[[nodiscard]] const broad_display_colour_spec &get_broad_spec() const;

		[[nodiscard]] const size_size_display_colour_map &get_clr_of_pdb_and_res() const;
	};

	const display_colour_opt & get_base_clr(const display_colour_spec &);

	const size_display_colour_map & get_clr_of_pdb(const display_colour_spec &);

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

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC_DISPLAY_COLOUR_SPEC_HPP
