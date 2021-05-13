/// \file
/// \brief The proximity_calculator class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PROXIMITY_CALCULATOR_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PROXIMITY_CALCULATOR_HPP

#include <functional>
#include <optional>
#include <vector>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/file/file_type_aliases.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/coord_linkage.hpp"
#include "cath/structure/structure_type_aliases.hpp"

namespace cath { namespace file { class pdb; } }
namespace cath { namespace file { class pdb_residue; } }

namespace cath {
	namespace file {
		namespace detail {

			/// \brief Represent a residue's index, key coord and maximum distance from the key coord to
			///        any of the other atoms
			///
			/// This struct is used in the proximity calculator
			struct res_index_key_coord_and_dist final {
				/// \brief The index of the residue in the source pdb
				size_t res_index;

				/// \brief The key coord (the CA atom where present or the first atom otherwise)
				geom::coord key_coord;

				/// \brief The furthest distance from the key_coord to any of the other atoms
				double furthest_atom_dist;

				/// \brief Standard ctor
				///
				/// This is defined (rather than just using aggregate initialisation) to
				/// allow the use of emplace_back in a vector of res_index_key_coord_and_dist
				res_index_key_coord_and_dist(const size_t &prm_res_index,         ///< The index of the residue in the source pdb
				                             geom::coord   prm_key_coord,         ///< The key coord (the CA atom where present or the first atom otherwise)
				                             const double &prm_furthest_atom_dist ///< The furthest distance from the key_coord to any of the other atoms
				                             ) : res_index          { prm_res_index              },
				                                 key_coord          { std::move( prm_key_coord ) },
				                                 furthest_atom_dist { prm_furthest_atom_dist     } {}
			};

			/// \brief Type alias for a vector of res_index_key_coord_and_dist values
			using res_index_key_coord_and_dist_vec = std::vector<res_index_key_coord_and_dist>;

		} // namespace detail


		/// \brief Index a PDB by storing a key coordinate for each residue and the maximum
		///        distance from that point to the other atoms in the residue to make it a bit
		///        faster to identify whether a point is close to the PDB
		///
		/// If necessary, it should be possible to further optimise this sort of calculation,
		/// eg if the proximity_calculator returns a distance rather than a bool, one large distance
		/// value could be used to rule out large chunks of the data that's spatially close
		/// to that query point (using spatial cells)
		class proximity_calculator final {
		private:
			/// \brief A reference to the source PDB
			std::reference_wrapper<const pdb> source_pdb;

			/// \brief The res_index_key_coord_and_dist values calculated for the PDB's residues
			detail::res_index_key_coord_and_dist_vec obj_res_indices_key_coords_and_dists;

			static detail::res_index_key_coord_and_dist_vec get_res_indices_key_coords_and_dists(const pdb &,
			                                                                                     const chop::region_vec_opt &);

		public:
			proximity_calculator(const pdb &,
			                     const chop::region_vec_opt & = ::std::nullopt);

			proximity_calculator(const proximity_calculator &) = default;
			proximity_calculator(proximity_calculator &&) noexcept = default;
			proximity_calculator & operator=(const proximity_calculator &) = default;
			proximity_calculator & operator=(proximity_calculator &&) noexcept = default;

			bool is_within_distance(const geom::coord &,
			                        const double &) const;
		};

		bool is_within_distance(const proximity_calculator &,
		                        const pdb_residue &,
		                        const double &);

		chain_label_set nearby_dna_rna_chain_labels(const pdb &,
		                                            const proximity_calculator &,
		                                            const double &);

		chain_label_set nearby_dna_rna_chain_labels(const pdb &,
		                                            const proximity_calculator &,
		                                            const doub_opt &);

		void restrict_to_linkage_proximate(geom::coord_coord_linkage_pair_vec &,
		                                   const proximity_calculator &,
		                                   const double &,
		                                   const double &);

		geom::coord_coord_linkage_pair_vec restrict_to_linkage_proximate_copy(geom::coord_coord_linkage_pair_vec,
		                                                                      const proximity_calculator &,
		                                                                      const double &,
		                                                                      const double &);

		pdb_residue_vec get_nearby_post_ter_res_atoms(const pdb &,
		                                              const proximity_calculator &,
		                                              const double &,
		                                              const double &);

		pdb_residue_vec get_nearby_post_ter_res_atoms(const pdb &,
		                                              const proximity_calculator &,
		                                              const doub_opt &,
		                                              const double &);

	} // namespace file
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PROXIMITY_CALCULATOR_HPP
