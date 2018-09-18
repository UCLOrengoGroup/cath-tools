/// \file
/// \brief The multi_struc_res_rep_pair class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_RES_PAIR_MULTI_STRUC_RES_REP_PAIR_H
#define _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_RES_PAIR_MULTI_STRUC_RES_REP_PAIR_H

#include "scan/detail/res_pair/res_pair_core.hpp"
#include "scan/detail/res_pair_dirn/res_pair_dirn.hpp"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Store data on a from/to pair of residue representatives from one of
			///        potentially multiple structures for the purpose of fast scanning
			///
			/// The indices that appear as direct members refer to the corresponding
			/// scan_multi_structure_data and scan_structure_data and mean nothing outside
			/// that context.
			class multi_struc_res_rep_pair final {
			private:
				/// \brief The core of data that scanning code uses for all res_pairs
				res_pair_core the_core;

				/// \brief The index of the structure from which this from/to residue representative pair comes
				///
				/// This is the index in the corresponding scan_multi_structure_data
				index_type structure_index;

				/// \brief The index of the from residue representative in the corresponding scan_structure_data
				///
				/// NOTE: This is not the same as the residue's index in the structure,
				///       which is typically much higher.
				res_rep_index_type from_res_rep_index;

				/// \brief The index of the to residue representative in the corresponding scan_structure_data
				///
				/// NOTE: This is not the same as the residue's index in the structure,
				///       which is typically much higher.
				res_rep_index_type to_res_rep_index;

			public:
				multi_struc_res_rep_pair(res_pair_core,
				                         const index_type &,
				                         const res_rep_index_type &,
				                         const res_rep_index_type &);

				const res_pair_core & get_res_pair_core() const;
				const index_type & get_structure_index() const;
				const res_rep_index_type & get_from_res_rep_index() const;
				const res_rep_index_type & get_to_res_rep_index() const;
			};

			const view_base_type & get_view_x(const multi_struc_res_rep_pair &);
			const view_base_type & get_view_y(const multi_struc_res_rep_pair &);
			const view_base_type & get_view_z(const multi_struc_res_rep_pair &);

			/// \brief Ctor from a res_pair_core, the index of the structure and the indices of the from/to residue reps
			inline multi_struc_res_rep_pair::multi_struc_res_rep_pair(res_pair_core             prm_core,               ///< The core properties of the res_pair
			                                                          const index_type         &prm_structure_index,    ///< The index of the structure in the scan_multi_structure_data (name?)
			                                                          const res_rep_index_type &prm_from_res_rep_index, ///< The rep index of the from rep-residue
			                                                          const res_rep_index_type &prm_to_res_rep_index    ///< The rep index of the to   rep-residue
			                                                          ) : the_core           { std::move( prm_core )  },
			                                                              structure_index    { prm_structure_index    },
			                                                              from_res_rep_index { prm_from_res_rep_index },
			                                                              to_res_rep_index   { prm_to_res_rep_index   } {
			}

			/// \brief TODOCUMENT
			inline const res_pair_core & multi_struc_res_rep_pair::get_res_pair_core() const {
				return the_core;
			}

			/// \brief TODOCUMENT
			inline const index_type & multi_struc_res_rep_pair::get_structure_index() const {
				return structure_index;
			}

			/// \brief TODOCUMENT
			inline const res_rep_index_type & multi_struc_res_rep_pair::get_from_res_rep_index() const {
				return from_res_rep_index;
			}

			/// \brief TODOCUMENT
			inline const res_rep_index_type & multi_struc_res_rep_pair::get_to_res_rep_index() const {
				return to_res_rep_index;
			}

			/// \brief Whether a res_pair's from-residue comes before (INCREASE) or after (DECREASE) its to-residue
			///
			/// \relates multi_struc_res_rep_pair
			inline res_pair_dirn direction(const multi_struc_res_rep_pair &prm_res_pair ///< The single_struc_res_pair to query
			                               ) {
#ifndef NDEBUG
				if ( prm_res_pair.get_from_res_rep_index() == prm_res_pair.get_to_res_rep_index() ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("direction() cannot process res_pairs with matching to/from rep residues"));
				}
#endif
				return ( prm_res_pair.get_from_res_rep_index() < prm_res_pair.get_to_res_rep_index() ) ? res_pair_dirn::INCREASE
				                                                                                       : res_pair_dirn::DECREASE;
			}

			/// \brief Return whether two res_pairs both have the same direction
			///        (ie both have from-residue before to-residue or both have from-residue after to-residue)
			///
			/// \relates multi_struc_res_rep_pair
			inline bool same_direction(const multi_struc_res_rep_pair &prm_res_pair_a, ///< The first  res_pair to compare
			                           const multi_struc_res_rep_pair &prm_res_pair_b  ///< The second res_pair to compare
			                           ) {
				return direction( prm_res_pair_a ) == direction( prm_res_pair_b );
			}

			/// \brief TODOCUMENT
			///
			/// \todo Add a version that takes a protein, a structure index, a pair of striders and a pair of indices
			///
			/// \relates multi_struc_res_rep_pair
			inline multi_struc_res_rep_pair make_multi_struc_res_rep_pair(const residue            &prm_from_residue,       ///< TODOCUMENT
			                                                              const residue            &prm_to_residue,         ///< TODOCUMENT
			                                                              const index_type         &prm_structure_index,    ///< TODOCUMENT
			                                                              const res_rep_index_type &prm_from_res_rep_index, ///< TODOCUMENT
			                                                              const res_rep_index_type &prm_to_res_rep_index    ///< TODOCUMENT
			                                                              ) {
				return {
					make_res_pair_core( prm_from_residue, prm_to_residue ),
					prm_structure_index,
					prm_from_res_rep_index,
					prm_to_res_rep_index
				};
			}

			/// \brief TODOCUMENT
			inline const view_base_type & get_view_x(const multi_struc_res_rep_pair &prm_res_pair ///< TODOCUMENT
			                                         ) {
				return get_view_x( prm_res_pair.get_res_pair_core() );
			}

			/// \brief TODOCUMENT
			inline const view_base_type & get_view_y(const multi_struc_res_rep_pair &prm_res_pair ///< TODOCUMENT
			                                         ) {
				return get_view_y( prm_res_pair.get_res_pair_core() );
			}

			/// \brief TODOCUMENT
			inline const view_base_type & get_view_z(const multi_struc_res_rep_pair &prm_res_pair ///< TODOCUMENT
			                                         ) {
				return get_view_z( prm_res_pair.get_res_pair_core() );
			}

			std::ostream & operator<<(std::ostream &,
			                          const multi_struc_res_rep_pair &);
		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
