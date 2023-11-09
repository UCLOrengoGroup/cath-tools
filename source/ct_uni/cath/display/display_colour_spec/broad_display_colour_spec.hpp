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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC_BROAD_DISPLAY_COLOUR_SPEC_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC_BROAD_DISPLAY_COLOUR_SPEC_HPP

#include "cath/display/colour_category.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"

// clang-format off
namespace cath { class viewer; }
namespace cath::file { class pdb_list; }
// clang-format on

namespace cath {

	// /// \brief TODOCUMENT
	// using region_display_colour_pair              = std::pair<region, display_colour>;

	// /// \brief TODOCUMENT
	// using region_display_colour_pair_vec          = std::vector<region_display_colour_pair>;

	// /// \brief TODOCUMENT
	// using size_region_display_colour_pair_vec_map = std::map<size_t, region_display_colour_pair_vec>;


	/// \brief Represent the colouring of structures at a broad (ie not residue-specific) level
	class broad_display_colour_spec final {
	private:
		/// \brief The main overall colour to use
		display_colour_opt                      base_clr;

		/// \brief The colours to use for the individual PDBs
		size_display_colour_map                 clr_of_pdb;

		// /// \brief The colours to use for regions of individual PDBs
		// size_region_display_colour_pair_vec_map clr_of_region;

	public:
		broad_display_colour_spec & colour_base(const display_colour &,
		                                        const bool & = false);

		broad_display_colour_spec & colour_pdb(const size_t &,
		                                       const display_colour &,
		                                       const bool & = false);

		// broad_display_colour_spec & colour_region(const size_t &,
		//                                           const region &,
		//                                           const display_colour &);

		[[nodiscard]] const display_colour_opt &get_base_clr() const;

		[[nodiscard]] const size_display_colour_map &get_clr_of_pdb() const;

		// const size_region_display_colour_pair_vec_map & get_clr_of_regions() const;
	};

	bool has_base_colour(const broad_display_colour_spec &);
	
	display_colour get_base_colour(const broad_display_colour_spec &);

	display_colour_opt get_clr_of_pdb_index(const broad_display_colour_spec &,
	                                        const size_t &);

	display_colour_vec get_pdb_colours(const broad_display_colour_spec &);

	size_vec get_pdbs_of_colour(const broad_display_colour_spec &,
	                            const display_colour &);

	void define_all_colours( const display_colour_vec &,
	                         const viewer &,
	                         ::std::ostream &,
	                         const colour_category & );

	namespace detail {
		void colour_base_and_pdbs_impl(const display_colour_vec &,
		                               const broad_display_colour_spec &,
		                               const viewer &,
		                               const file::pdb_list &,
		                               const str_vec &,
		                               const colour_category &,
		                               std::ostream &);
	} // namespace detail

	void colour_viewer_with_spec(const broad_display_colour_spec &,
	                             const viewer &,
	                             const file::pdb_list &,
	                             const str_vec &,
	                             std::ostream &);

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_DISPLAY_COLOUR_SPEC_BROAD_DISPLAY_COLOUR_SPEC_HPP
