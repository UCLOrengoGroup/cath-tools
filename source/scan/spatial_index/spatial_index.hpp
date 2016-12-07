/// \file
/// \brief The spatial_index class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_SPATIAL_INDEX_SPATIAL_INDEX_H
#define _CATH_TOOLS_SOURCE_SCAN_SPATIAL_INDEX_SPATIAL_INDEX_H

#include "common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/size_t_literal.hpp"
#include "file/pdb/pdb.hpp"
#include "scan/detail/scan_index_store/scan_index_lattice_store.hpp"
#include "scan/detail/scan_type_aliases.hpp"
#include "scan/res_pair_keyer/res_pair_keyer.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_x_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_y_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_z_keyer_part.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/protein/residue.hpp"

using namespace cath::common::literals;

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		struct simple_locn_index final {
			using view_t = detail::view_base_type;

			/// \brief TODOCUMENT
			detail::view_base_type view_x;

			/// \brief TODOCUMENT
			detail::view_base_type view_y;

			/// \brief TODOCUMENT
			detail::view_base_type view_z;

			/// \brief TODOCUMENT
			unsigned int index;
		};

		/// \brief Convenience function to get the x component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const detail::view_base_type & get_view_x(const simple_locn_index &arg_locn_index ///< The simple_locn_index to query
		                                                           ) {
		        return arg_locn_index.view_x;
		}

		/// \brief Convenience function to get the y component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const detail::view_base_type & get_view_y(const simple_locn_index &arg_locn_index ///< The simple_locn_index to query
		                                                           ) {
		        return arg_locn_index.view_y;
		}

		/// \brief Convenience function to get the z component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const detail::view_base_type & get_view_z(const simple_locn_index &arg_locn_index ///< The simple_locn_index to query
		                                                           ) {
		        return arg_locn_index.view_z;
		}

		inline detail::view_base_type get_squared_distance(const simple_locn_index &arg_locn_index_a, ///< The simple_locn_index to query
		                                                   const simple_locn_index &arg_locn_index_b  ///< The simple_locn_index to query
		                                                   ) {
			return (
				(
					( get_view_x( arg_locn_index_a ) - get_view_x( arg_locn_index_b ) )
					*
					( get_view_x( arg_locn_index_a ) - get_view_x( arg_locn_index_b ) )
				)
				+
				(
					( get_view_y( arg_locn_index_a ) - get_view_y( arg_locn_index_b ) )
					*
					( get_view_y( arg_locn_index_a ) - get_view_y( arg_locn_index_b ) )
				)
				+
				(
					( get_view_z( arg_locn_index_a ) - get_view_z( arg_locn_index_b ) )
					*
					( get_view_z( arg_locn_index_a ) - get_view_z( arg_locn_index_b ) )
				)
			);
		}

		inline bool are_within_distance(const simple_locn_index &arg_locn_index_a, ///< The simple_locn_index to query
		                                const simple_locn_index &arg_locn_index_b, ///< The simple_locn_index to query
		                                const float &arg_max_dist,                 ///< TODOCUMENT
		                                const float &arg_max_squared_dist          ///< TODOCUMENT
		                                ) {
			const auto dist_x = get_view_x( arg_locn_index_a ) - get_view_x( arg_locn_index_b );
			if ( dist_x > arg_max_dist ) {
				return false;
			}
			const auto dist_y = get_view_y( arg_locn_index_a ) - get_view_y( arg_locn_index_b );
			if ( dist_y > arg_max_dist ) {
				return false;
			}
			const auto dist_z = get_view_z( arg_locn_index_a ) - get_view_z( arg_locn_index_b );
			if ( dist_z > arg_max_dist ) {
				return false;
			}
			return (
				( dist_x * dist_x )
				+
				( dist_y * dist_y )
				+
				( dist_z * dist_z )
			) < arg_max_squared_dist;
		}

		/// \brief TODOCUMENT
		inline simple_locn_index make_simple_locn_index(const geom::coord  &arg_coord, ///< TODOCUMENT
		                                                const unsigned int &arg_index  ///< TODOCUMENT
		                                                ) {
			return {
				debug_numeric_cast< detail::view_base_type>( arg_coord.get_x() ),
				debug_numeric_cast< detail::view_base_type>( arg_coord.get_y() ),
				debug_numeric_cast< detail::view_base_type>( arg_coord.get_z() ),
				arg_index
			};
		}

		/// \brief TODOCUMENT
		inline simple_locn_index make_simple_locn_index_of_ca(const residue      &arg_res,   ///< TODOCUMENT
		                                                      const unsigned int &arg_index  ///< TODOCUMENT
		                                                      ) {
			return make_simple_locn_index(
				arg_res.get_carbon_alpha_coord(),
				arg_index
			);
		}

		/// \brief TODOCUMENT
		inline simple_locn_index make_simple_locn_index_of_ca(const file::pdb_residue &arg_res,   ///< TODOCUMENT
		                                                      const unsigned int      &arg_index  ///< TODOCUMENT
		                                                      ) {
			return make_simple_locn_index(
				get_carbon_alpha_coord( arg_res ),
				arg_index
			);
		}

		/// \brief TODOCUMENT
		struct simple_locn_crit final {

			/// \brief TODOCUMENT
			detail::view_base_type maximum_squared_distance;
		};

		inline std::string to_string(const simple_locn_index &arg_simple_locn_index
		                             ) {
			return "simple_locn_index[ ("
				+ ::std::to_string( arg_simple_locn_index.view_x )
				+ ", "
				+ ::std::to_string( arg_simple_locn_index.view_y )
				+ ", "
				+ ::std::to_string( arg_simple_locn_index.view_z )
				+ "), "
				+ ::std::to_string( arg_simple_locn_index.index )
				+ "]";
		}

		/// \brief TODOCUMENT
		inline detail::view_base_type get_maximum_distance(const simple_locn_crit &arg_crit ///< TODOCUMENT
		                                                   ) {
			return std::sqrt( arg_crit.maximum_squared_distance );
		}

		/// \brief TODOCUMENT
		using locn_index_store = detail::scan_index_lattice_store<std::tuple<short, short, short>, std::vector<simple_locn_index> >;

		/// \brief TODOCUMENT
		using simple_locn_x_keyer_part = detail::axis_keyer_part< detail::res_pair_view_x_keyer_part_spec< simple_locn_index, simple_locn_crit > >;

		/// \brief TODOCUMENT
		using simple_locn_y_keyer_part = detail::axis_keyer_part< detail::res_pair_view_y_keyer_part_spec< simple_locn_index, simple_locn_crit > >;

		/// \brief TODOCUMENT
		using simple_locn_z_keyer_part = detail::axis_keyer_part< detail::res_pair_view_z_keyer_part_spec< simple_locn_index, simple_locn_crit > >;

		template <typename Fn>
		void scan_sparse_lattice(const locn_index_store &arg_store,     ///< TODOCUMENT
		                         const protein          &arg_protein,   ///< TODOCUMENT
		                         const float            &arg_cell_size, ///< TODOCUMENT
		                         const float            &arg_max_dist,  ///< TODOCUMENT
		                         Fn                      arg_fn         ///< TODOCUMENT
		                         ) {
			const auto the_keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ arg_cell_size },
				simple_locn_y_keyer_part{ arg_cell_size },
				simple_locn_z_keyer_part{ arg_cell_size }
			);

			const float max_squared_dist = arg_max_dist * arg_max_dist;
			const simple_locn_crit the_crit{ max_squared_dist };

			for (const size_t &the_res_idx : boost::irange( 0_z, arg_protein.get_length() ) ) {
				const auto &the_res = arg_protein.get_residue_ref_of_index( the_res_idx );
				const auto  data    = make_simple_locn_index_of_ca( the_res, debug_numeric_cast<unsigned int>( the_res_idx ) );
				for (const auto &key : common::cross( the_keyer.make_close_keys( data, the_crit ) ) ) {
					if ( arg_store.has_matches( key ) ) {
						for (const simple_locn_index &eg : arg_store.find_matches( key ) ) {
							// if ( are_within_distance( eg, data, arg_max_dist, max_squared_dist ) ) {
							if ( get_squared_distance( eg, data ) < max_squared_dist ) {
								arg_fn( data, eg );
							}
						}
					}
				}
			}
		}

		template <typename Fn>
		void scan_sparse_lattice(const locn_index_store &arg_store,     ///< TODOCUMENT
		                         const file::pdb        &arg_pdb,       ///< TODOCUMENT
		                         const float            &arg_cell_size, ///< TODOCUMENT
		                         const float            &arg_max_dist,  ///< TODOCUMENT
		                         Fn                      arg_fn         ///< TODOCUMENT
		                         ) {
			const auto the_keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ arg_cell_size },
				simple_locn_y_keyer_part{ arg_cell_size },
				simple_locn_z_keyer_part{ arg_cell_size }
			);

			const float max_squared_dist = arg_max_dist * arg_max_dist;
			const simple_locn_crit the_crit{ max_squared_dist };

			for (const size_t &the_res_idx : boost::irange( 0_z, arg_pdb.get_num_residues() ) ) {
				const auto &the_res = arg_pdb.get_residue_cref_of_index__backbone_unchecked( the_res_idx );
				const auto  data    = make_simple_locn_index_of_ca( the_res, debug_numeric_cast<unsigned int>( the_res_idx ) );
				for (const auto &key : common::cross( the_keyer.make_close_keys( data, the_crit ) ) ) {
					if ( arg_store.has_matches( key ) ) {
						for (const simple_locn_index &eg : arg_store.find_matches( key ) ) {
							// if ( are_within_distance( eg, data, arg_max_dist, max_squared_dist ) ) {
							if ( get_squared_distance( eg, data ) < max_squared_dist ) {
								arg_fn( data, eg );
							}
						}
					}
				}
			}
		}

		template <typename Fn>
		void scan_dense_lattice(const locn_index_store &arg_store,     ///< TODOCUMENT
		                        const protein          &arg_protein,   ///< TODOCUMENT
		                        const float            &arg_cell_size, ///< TODOCUMENT
		                        const float            &arg_max_dist,  ///< TODOCUMENT
		                        Fn                      arg_fn         ///< TODOCUMENT
		                        ) {
			const auto the_keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ arg_cell_size },
				simple_locn_y_keyer_part{ arg_cell_size },
				simple_locn_z_keyer_part{ arg_cell_size }
			);

			const float max_squared_dist = arg_max_dist * arg_max_dist;
			// const simple_locn_crit the_crit{ max_squared_dist };

			for (const size_t &the_res_idx : boost::irange( 0_z, arg_protein.get_length() ) ) {
				const auto &the_res = arg_protein.get_residue_ref_of_index( the_res_idx );
				const auto  data    = make_simple_locn_index_of_ca( the_res, debug_numeric_cast<unsigned int>( the_res_idx ) );
				const auto &key = the_keyer.make_key( data );
				for (const simple_locn_index &eg : arg_store.find_matches( key ) ) {
					// if ( are_within_distance( eg, data, arg_max_dist, max_squared_dist ) ) {
					if ( get_squared_distance( eg, data ) < max_squared_dist ) {
						arg_fn( data, eg );
					}
				}
			}
		}

		locn_index_store make_sparse_lattice(const protein &,
		                                     const float &);

		locn_index_store make_dense_lattice(const protein &,
		                                    const float &,
		                                    const float &);


		locn_index_store make_sparse_lattice(const file::pdb &,
		                                     const float &);

		locn_index_store make_dense_lattice(const file::pdb &,
		                                    const float &,
		                                    const float &);


	} // namespace scan
} // namespace cath

#endif
