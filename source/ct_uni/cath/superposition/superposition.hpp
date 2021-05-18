/// \file
/// \brief The superposition class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSITION_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSITION_HPP

#include <array>
#include <iosfwd>

#include <boost/operators.hpp>
#include <boost/range/size.hpp>
#include <boost/tuple/tuple.hpp>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/rotation.hpp"

namespace cath { namespace file { class pdb; } }
namespace cath { namespace geom { class coord_list; } }

namespace cath {
	namespace sup {

		/// \brief The centres of gravity and the rotation matrix arising from superposing structures (coord_lists)
		///
		/// For now the object is read-only after construction and the coord_list objects
		/// must be provided at construction. This keeps things simpler.
		///
		/// For each structure, the superposition stores a translation and rotation that should be performed on that structure
		/// to move it into the desired location.
		///
		/// Note that the translation must be performed before the rotation.
		class superposition final : private boost::equality_comparable<superposition>  {
		private:
			/// \brief TODOCUMENT
			geom::coord_vec    translations;

			/// \brief TODOCUMENT
			geom::rotation_vec rotations;

			// Find the translation and rotation that will make the first coord_list best fit the second in the seconds current position
			static geom::coord_rot_pair fit_second_to_first(const geom::coord_list &,
			                                                const geom::coord_list &);

		public:
			using indices_and_coord_lists_type = std::tuple<size_t, geom::coord_list, size_t, geom::coord_list>;

			explicit superposition(const std::vector<indices_and_coord_lists_type> &,
			                       const size_t         & = 0,
			                       const geom::coord    & = geom::ORIGIN_COORD,
			                       const geom::rotation & = geom::rotation::IDENTITY_ROTATION());
			superposition(geom::coord_vec,
			              geom::rotation_vec);

			[[nodiscard]] size_t                get_num_entries() const;
			[[nodiscard]] const geom::coord &   get_translation_of_index( const size_t & ) const;
			[[nodiscard]] const geom::rotation &get_rotation_of_index( const size_t & ) const;

			superposition & post_translate(const geom::coord &);
			superposition & post_rotate(const geom::rotation &);

			static constexpr size_t          NUM_ENTRIES_IN_PAIRWISE_SUPERPOSITION     = 2;
			static constexpr size_t          INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION  = 0;
			static constexpr size_t          INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION = 1;
			static constexpr auto            SUPERPOSITION_CHAIN_CHARS = common::make_array(
				'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
				'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'
			);
			static const     chain_label_vec SUPERPOSITION_CHAIN_LABELS;
		};

		void post_translate_and_rotate(superposition &,
		                               const geom::coord &,
		                               const geom::rotation & = geom::rotation::IDENTITY_ROTATION() );

		superposition post_translate_and_rotate_copy(superposition,
		                                             const geom::coord &,
		                                             const geom::rotation & = geom::rotation::IDENTITY_ROTATION() );

		void transform(const superposition &,
		               const size_t &,
		               geom::coord &);

		geom::coord transform_copy(const superposition &,
		                           const size_t &,
		                           geom::coord);

		void transform(const superposition &,
		               const size_t &,
		               geom::coord_list &);

		geom::coord_list transform_copy(const superposition &,
		                                const size_t &,
		                                geom::coord_list);

		void transform(const superposition &,
		               const size_t &,
		               geom::coord_list_vec &);

		geom::coord_list_vec transform_copy(const superposition &,
		                                    const size_t &,
		                                    geom::coord_list_vec);

		double superposed_distance(const superposition &,
		                           const size_t &,
		                           geom::coord,
		                           const size_t &,
		                           geom::coord);

		double calc_rmsd_between_superposed_entries(const superposition &,
		                                            const size_t &,
		                                            geom::coord_list,
		                                            const size_t &,
		                                            geom::coord_list);

		bool operator==(const superposition &,
		                const superposition &);

		std::ostream & operator<<(std::ostream &,
		                          const superposition &);

		bool are_close(const superposition &,
		               const superposition &);

		void write_superposition(std::ostream &,
		                         const superposition &);

		superposition read_superposition(std::istream &);

		superposition create_pairwise_superposition(const geom::coord_list &,
		                                            const geom::coord_list &,
		                                            const bool           & = true,
		                                            const geom::coord    & = geom::ORIGIN_COORD,
		                                            const geom::rotation & = geom::rotation::IDENTITY_ROTATION());

		void superpose_second_coords_to_first(const geom::coord_list &,
		                                      geom::coord_list &);

		geom::coord_list superpose_copy_second_coords_to_first(const geom::coord_list &,
		                                                       geom::coord_list);

		void superpose_second_coords_to_first(const geom::coord_list_vec &,
		                                      geom::coord_list_vec &);

		geom::coord_list_vec superpose_copy_second_coords_to_first(const geom::coord_list_vec &,
		                                                           geom::coord_list_vec);

		double calc_pairwise_superposition_rmsd(const geom::coord_list &,
		                                        const geom::coord_list &);

		void check_superposition_is_pairwise(const superposition &);
//		double get_rmsd_from_pairwise_superposition(const superposition &);

		superposition make_identity_superposition(const size_t &);

		/// Make an identity superpostion with a number of entries to match the size() of the specified argument
		template <typename T>
		superposition make_identity_superposition_of(T &&arg ///< The object whose size should be used for the number of entries
		                                             ) {
			return make_identity_superposition( ::boost::size( ::std::forward<T>( arg ) ) );
		}

	} // namespace sup
} // namespace cath
#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSITION_HPP
