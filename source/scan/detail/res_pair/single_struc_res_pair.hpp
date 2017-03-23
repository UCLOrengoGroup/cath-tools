/// \file
/// \brief The single_struc_res_pair class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_DETAIL_RES_PAIR_SINGLE_STRUC_RES_PAIR_H
#define _CATH_TOOLS_SOURCE_SCAN_DETAIL_RES_PAIR_SINGLE_STRUC_RES_PAIR_H

#include "scan/detail/res_pair/res_pair_core.hpp"
#include "scan/detail/res_pair/functions/res_pair_core_functions.hpp"
#include "scan/detail/res_pair_dirn/res_pair_dirn.hpp"

#include <iosfwd>

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Store data on a from/to pair of residues from a known single structure
			///        for the purpose of fast scanning
			///
			/// This is useful for representing the list of residue pairs that neighbour a pair of rep residues
			class single_struc_res_pair final {
			private:
				/// \brief The core properties of the res_pair
				res_pair_core the_core;

				/// \brief The from-residue index (the absolute index in the source structure)
				index_type    from_res_idx;

				/// \brief The to-residue index (the absolute index in the source structure)
				index_type    to_res_idx;

			public:
				single_struc_res_pair();
				single_struc_res_pair(res_pair_core,
				                      const index_type &,
				                      const index_type &);

				single_struc_res_pair(const view_type &,
				                      const frame_quat_rot &,
				                      const angle_type &,
				                      const angle_type &,
				                      const angle_type &,
				                      const angle_type &,
				                      const index_type &,
				                      const index_type &);

				const res_pair_core & get_res_pair_core() const;
				const index_type & get_from_res_idx() const;
				const index_type & get_to_res_idx() const;

				static constexpr index_type DUMMY_INDEX_VALUE = std::numeric_limits<index_type>::max();
			};

			/// \brief Ctor to build a dummy res_pair
			inline single_struc_res_pair::single_struc_res_pair() : the_core(),
			                                                        from_res_idx( DUMMY_INDEX_VALUE ),
			                                                        to_res_idx  ( DUMMY_INDEX_VALUE ) {
			}

			/// \brief Ctor from a res_pair_core and the indices of the from/to residues
			inline single_struc_res_pair::single_struc_res_pair(res_pair_core     arg_res_pair_core, ///< The core properties of the res_pair
			                                                    const index_type &arg_from_res_idx,  ///< The from-residue index
			                                                    const index_type &arg_to_res_idx     ///< The to-residue   index
			                                                    ) : the_core     { std::move( arg_res_pair_core ) },
			                                                        from_res_idx { arg_from_res_idx               },
			                                                        to_res_idx   { arg_to_res_idx                 } {
			}

			/// \brief Ctor from all the parts
			inline single_struc_res_pair::single_struc_res_pair(const view_type      &arg_view,           ///< The view of the to_residue from the from_residue
			                                                    const frame_quat_rot &arg_frame,          ///< The coordinate frame of the from_residue, as determined by its core atoms
			                                                    const angle_type     &arg_from_phi_angle, ///< The phi angle of the from_residue
			                                                    const angle_type     &arg_from_psi_angle, ///< The psi angle of the from_residue
			                                                    const angle_type     &arg_to_phi_angle,   ///< The phi angle of the to_residue
			                                                    const angle_type     &arg_to_psi_angle,   ///< The psi angle of the to_residue
			                                                    const index_type     &arg_from_res_idx,   ///< The from-residue index
			                                                    const index_type     &arg_to_res_idx      ///< The to-residue   index
			                                                    ) : single_struc_res_pair(
			                                                        	res_pair_core(
			                                                        		arg_view,
			                                                        		arg_frame,
			                                                        		arg_from_phi_angle,
			                                                        		arg_from_psi_angle,
			                                                        		arg_to_phi_angle,
			                                                        		arg_to_psi_angle
			                                                        	),
			                                                        	arg_from_res_idx,
			                                                        	arg_to_res_idx
			                                                        ) {
			}

			/// \brief Getter for the res_pair core
			inline const res_pair_core & single_struc_res_pair::get_res_pair_core() const {
				return the_core;
			}

			/// \brief Getter for the from-residue index
			inline const index_type & single_struc_res_pair::get_from_res_idx() const {
				return from_res_idx;
			}

			/// \brief Getter for the to-residue index
			inline const index_type & single_struc_res_pair::get_to_res_idx() const {
				return to_res_idx;
			}

			/// \brief Whether a res_pair's from-residue comes before (INCREASE) or after (DECREASE) its to-residue
			///
			/// \relates single_struc_res_pair
			inline res_pair_dirn direction(const single_struc_res_pair &arg_res_pair ///< The single_struc_res_pair to query
			                               ) {
#ifndef NDEBUG
				if ( arg_res_pair.get_from_res_idx() == arg_res_pair.get_to_res_idx() ) {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("direction() cannot process res_pairs with matching to/from residues"));
				}
#endif
				return ( arg_res_pair.get_from_res_idx() < arg_res_pair.get_to_res_idx() ) ? res_pair_dirn::INCREASE
				                                                                           : res_pair_dirn::DECREASE;
			}

			/// \brief Return whether two res_pairs both have the same direction
			///        (ie both have from-residue before to-residue or both have from-residue after to-residue)
			///
			/// \relates single_struc_res_pair
			inline bool same_direction(const single_struc_res_pair &arg_res_pair_a, ///< The first  res_pair to compare
			                           const single_struc_res_pair &arg_res_pair_b  ///< The second res_pair to compare
			                           ) {
				return direction( arg_res_pair_a ) == direction( arg_res_pair_b );
			}

			/// \brief TODOCUMENT
			///
			/// \relates single_struc_res_pair
			inline bool is_dummy(const single_struc_res_pair &arg_res_pair ///< TODOCUMENT
			                     ) {
				return (
					arg_res_pair.get_from_res_idx() == single_struc_res_pair::DUMMY_INDEX_VALUE
					&&
					arg_res_pair.get_to_res_idx()   == single_struc_res_pair::DUMMY_INDEX_VALUE
				);
			}

			/// \brief TODOCUMENT
			///
			/// \relates single_struc_res_pair
			inline single_struc_res_pair make_single_res_pair(const residue    &arg_from_residue, ///< TODOCUMENT
			                                                  const residue    &arg_to_residue,   ///< TODOCUMENT
			                                                  const index_type &arg_from_index,   ///< TODOCUMENT
			                                                  const index_type &arg_to_index      ///< TODOCUMENT
			                                                  ) {
				return {
					make_res_pair_core( arg_from_residue, arg_to_residue ),
					arg_from_index,
					arg_to_index
				};
			}

			/// \brief TODOCUMENT
			///
			/// \relates single_struc_res_pair
			inline single_struc_res_pair make_single_res_pair(const protein    &arg_protein,    ///< TODOCUMENT
			                                                  const index_type &arg_from_index, ///< TODOCUMENT
			                                                  const index_type &arg_to_index    ///< TODOCUMENT
			                                                  ) {
				return make_single_res_pair(
					arg_protein.get_residue_ref_of_index( arg_from_index ),
					arg_protein.get_residue_ref_of_index( arg_to_index   ),
					arg_from_index,
					arg_to_index
				);
			}

			/// \brief Calculate the squared distance between the views of the two residue pairs
			///
			/// Each view is the location of the to_residue's carbon-beta atom as seen from
			/// the coordinate frame of the from_residue
			///
			/// \relates single_struc_res_pair
			inline double squared_distance(const single_struc_res_pair &arg_res_pair_a, ///< The first  res_pair
			                               const single_struc_res_pair &arg_res_pair_b  ///< The second res_pair
			                               ) {
				return squared_distance(
					arg_res_pair_a.get_res_pair_core(),
					arg_res_pair_b.get_res_pair_core()
				);
			}

			std::ostream & operator<<(std::ostream &,
			                          const single_struc_res_pair &);

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
