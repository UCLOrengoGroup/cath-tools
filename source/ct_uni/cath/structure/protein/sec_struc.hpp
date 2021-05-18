/// \file
/// \brief The sec_struc class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_SEC_STRUC_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_SEC_STRUC_HPP

#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/protein/sec_struc_type.hpp"
#include "cath/structure/structure_type_aliases.hpp"

#include <cstddef>

namespace cath {
	class sec_struc_planar_angles;

	/// \todo Move planar_angles out of this into a separate class that is stored in the protein because
	///       they aren't really properties of a sec_struc but of the relationships between sec_strucs.
	///
	/// \todo In fact, do planar_angles need to be stored at all or could they just be calculated dynamically?

	/// \brief Represent a secondary structure as used in the SSAP code
	class sec_struc final {
	private:
		/// \brief TODOCUMENT
		size_t         start_residue_num;

		/// \brief TODOCUMENT
		size_t         stop_residue_num;

		/// \brief TODOCUMENT
		sec_struc_type type;   ///< \brief Alpha-helix ('H') or beta-strand ('S').

		/// \brief TODOCUMENT
		geom::coord    midpoint;

		/// \brief TODOCUMENT
		geom::coord    unit_dirn;

		/// \brief The planar angles between this sec_struc and the others as parsed from the sec file
		///
		/// Note that the angles to preceding sec_strucs are just stored as copies of the angles from those
		/// preceding sec_struc - no negation is applied.
		///
		/// The entry for this sec_struc is stored with all zeroes.
		sec_struc_planar_angles_vec planar_angles;

	private:
		void check_planar_angles_index_is_valid(const size_t &) const;
		void check_sec_struc_type() const;

	public:
		sec_struc(const size_t &,
		          const size_t &,
		          const sec_struc_type &,
		          geom::coord,
		          geom::coord);

		void set_planar_angles(const sec_struc_planar_angles_vec &);

		[[nodiscard]] size_t                      get_start_residue_num() const;
		[[nodiscard]] size_t                      get_stop_residue_num() const;
		[[nodiscard]] sec_struc_type              get_type() const;
		[[nodiscard]] geom::coord                 get_midpoint() const;
		[[nodiscard]] geom::coord                 get_unit_dirn() const;
		[[nodiscard]] sec_struc_planar_angles_vec get_planar_angles() const;

		[[nodiscard]] const sec_struc_planar_angles &get_planar_angles_of_index( const size_t & ) const;
		[[nodiscard]] size_t                         get_num_planar_angles() const;

		// /// \todo Sort out indexing and then get rid of this because a blank secondary structure doesn't really make any sense
		// static const sec_struc NULL_SEC_STRUC;
	};

	geom::coord calculate_inter_sec_struc_vector(const sec_struc &,
	                                             const sec_struc &,
	                                             const sec_struc &);
} // namespace cath
#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_SEC_STRUC_HPP
