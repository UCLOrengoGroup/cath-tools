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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_SPATIAL_INDEX_SPATIAL_INDEX_H
#define _CATH_TOOLS_SOURCE_UNI_SCAN_SPATIAL_INDEX_SPATIAL_INDEX_H

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "cath/common/constexpr_ignore_unused.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_lattice_store.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/res_pair_from_to_index_keyer_part.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_x_keyer_part.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_y_keyer_part.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_z_keyer_part.hpp"
#include "cath/scan/spatial_index/simple_locn_index.hpp"
#include "cath/structure/protein/residue.hpp"

using namespace ::cath::common::literals;

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		inline simple_locn_index make_simple_locn_index_of_ca(const residue      &prm_res,   ///< TODOCUMENT
		                                                      const unsigned int &prm_index  ///< TODOCUMENT
		                                                      ) {
			return make_simple_locn_index(
				prm_res.get_carbon_alpha_coord(),
				prm_index
			);
		}

		/// \brief TODOCUMENT
		inline simple_locn_index make_simple_locn_index_of_ca(const file::pdb_residue &prm_res,   ///< TODOCUMENT
		                                                      const unsigned int      &prm_index  ///< TODOCUMENT
		                                                      ) {
			return make_simple_locn_index(
				get_carbon_alpha_coord( prm_res ),
				prm_index
			);
		}

		/// \brief TODOCUMENT
		using locn_index_store = detail::scan_index_lattice_store<std::tuple<short, short, short>, std::vector<simple_locn_index> >;

		// temporarily, from simple_locn_index, need:
		//    * key: tpl_elmnt_skip_t{}, value: index

		// then, from fuller_locn_index, need:
		//    * key:  tpl_elmnt_skip_t{}, value: from_index // vec<stores>
		//    * key:  tpl_elmnt_skip_t{}, value: to_index   // vec<stores>
		//    * key:  to_index,           value: from_index // store<>
		//    * key:  from_index,         value: to_index   // store<>

		/// \brief TODOCUMENT
		using simple_locn_x_keyer_part = detail::axis_keyer_part< detail::res_pair_view_x_keyer_part_spec< simple_locn_index, simple_locn_crit > >;

		/// \brief TODOCUMENT
		using simple_locn_y_keyer_part = detail::axis_keyer_part< detail::res_pair_view_y_keyer_part_spec< simple_locn_index, simple_locn_crit > >;

		/// \brief TODOCUMENT
		using simple_locn_z_keyer_part = detail::axis_keyer_part< detail::res_pair_view_z_keyer_part_spec< simple_locn_index, simple_locn_crit > >;

		template <typename Fn>
		void scan_sparse_lattice(const locn_index_store &prm_store,     ///< TODOCUMENT
		                         const protein          &prm_protein,   ///< TODOCUMENT
		                         const float            &prm_cell_size, ///< TODOCUMENT
		                         const float            &prm_max_dist,  ///< TODOCUMENT
		                         Fn                      prm_fn         ///< TODOCUMENT
		                         ) {
			const auto the_keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ prm_cell_size },
				simple_locn_y_keyer_part{ prm_cell_size },
				simple_locn_z_keyer_part{ prm_cell_size },
				res_pair_from_to_index_keyer_part{}
			);

			const float max_squared_dist = prm_max_dist * prm_max_dist;
			const simple_locn_crit the_crit{ max_squared_dist };

			for (const size_t &the_res_idx : common::indices( prm_protein.get_length() ) ) {
				const auto &the_res = prm_protein.get_residue_ref_of_index( the_res_idx );
				const auto  data    = make_simple_locn_index_of_ca( the_res, debug_numeric_cast<unsigned int>( the_res_idx ) );
				for (const auto &key : common::cross( the_keyer.make_close_keys( data, the_crit ) ) ) {
					if ( prm_store.has_matches( key ) ) {
						for (const simple_locn_index &eg : prm_store.find_matches( key ) ) {
							// if ( are_within_distance( eg, data, prm_max_dist, max_squared_dist ) ) {
							if ( get_squared_distance( eg, data ) < max_squared_dist ) {
								prm_fn( data, eg );
							}
						}
					}
				}
			}
		}

		template <typename Fn>
		void scan_sparse_lattice(const locn_index_store &prm_store,     ///< TODOCUMENT
		                         const file::pdb        &prm_pdb,       ///< TODOCUMENT
		                         const float            &prm_cell_size, ///< TODOCUMENT
		                         const float            &prm_max_dist,  ///< TODOCUMENT
		                         Fn                      prm_fn         ///< TODOCUMENT
		                         ) {
			const auto the_keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ prm_cell_size },
				simple_locn_y_keyer_part{ prm_cell_size },
				simple_locn_z_keyer_part{ prm_cell_size },
				res_pair_from_to_index_keyer_part{}
			);

			const float max_squared_dist = prm_max_dist * prm_max_dist;
			const simple_locn_crit the_crit{ max_squared_dist };

			for (const size_t &the_res_idx : common::indices( prm_pdb.get_num_residues() ) ) {
				const auto &the_res = prm_pdb.get_residue_of_index__backbone_unchecked( the_res_idx );
				const auto  data    = make_simple_locn_index_of_ca( the_res, debug_numeric_cast<unsigned int>( the_res_idx ) );
				for (const auto &key : common::cross( the_keyer.make_close_keys( data, the_crit ) ) ) {
					if ( prm_store.has_matches( key ) ) {
						for (const simple_locn_index &eg : prm_store.find_matches( key ) ) {
							// if ( are_within_distance( eg, data, prm_max_dist, max_squared_dist ) ) {
							if ( get_squared_distance( eg, data ) < max_squared_dist ) {
								prm_fn( data, eg );
							}
						}
					}
				}
			}
		}

		template <typename Fn>
		void scan_dense_lattice(const locn_index_store &prm_store,     ///< TODOCUMENT
		                        const protein          &prm_protein,   ///< TODOCUMENT
		                        const float            &prm_cell_size, ///< TODOCUMENT
		                        const float            &prm_max_dist,  ///< TODOCUMENT
		                        Fn                      prm_fn         ///< TODOCUMENT
		                        ) {
			const auto the_keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ prm_cell_size },
				simple_locn_y_keyer_part{ prm_cell_size },
				simple_locn_z_keyer_part{ prm_cell_size }
			);

			const float max_squared_dist = prm_max_dist * prm_max_dist;
			// const simple_locn_crit the_crit{ max_squared_dist };

			for (const size_t &the_res_idx : common::indices( prm_protein.get_length() ) ) {
				const auto &the_res = prm_protein.get_residue_ref_of_index( the_res_idx );
				const auto  data    = make_simple_locn_index_of_ca( the_res, debug_numeric_cast<unsigned int>( the_res_idx ) );
				const auto &key = the_keyer.make_key( data );
				for (const simple_locn_index &eg : prm_store.find_matches( key ) ) {
					// if ( are_within_distance( eg, data, prm_max_dist, max_squared_dist ) ) {
					if ( get_squared_distance( eg, data ) < max_squared_dist ) {
						prm_fn( data, eg );
					}
				}
			}
		}

		locn_index_store make_sparse_lattice(const protein &,
		                                     const float &,
		                                     const float &);

		locn_index_store make_dense_lattice(const protein &,
		                                    const float &,
		                                    const float &);


		locn_index_store make_sparse_lattice(const file::pdb &,
		                                     const float &,
		                                     const float &);

		locn_index_store make_dense_lattice(const file::pdb &,
		                                    const float &,
		                                    const float &);


		/// \brief Key part generator for simple_locn_index that extracts the index
		class simple_locn_index_keyer_part final {
		public:
			/// \brief The type for the value extracted from the simple_locn_index
			using value_t           = decltype( simple_locn_index::index );

			/// \brief The type for the cell index
			using cell_index_t      = unsigned int;

			/// \brief The type for the range of indices of close cells
			using cell_index_list_t = std::array<cell_index_t, 1>;

			/// \brief The type for the search radius
			using search_radius_t   = unsigned int;

			/// \brief Get a short name that describes this key part
			std::string get_name() const {
				return "index";
			}

			/// \brief Extract the relevant value from the specified simple_locn_index
			constexpr value_t get_value(const simple_locn_index &prm_locn_index ///< TODOCUMENT
			                            ) const {
				return prm_locn_index.index;
			}

			/// \brief Extract the search radius from the specified criteria
			constexpr search_radius_t get_search_radius(const simple_locn_crit &/*prm_criteria*/ ///< The criteria defining what is considered a match
			                                            ) const {
				return 0;
			}

			/// \brief Generate the key part for the specified value
			constexpr cell_index_t key_part(const value_t &prm_value ///< The value for which the key_part should be extracted
			                                ) const {
				return prm_value;
			}

			/// \brief Generate a list of all key parts for all conceivable simple_locn_indexs that would match the specified value
			///        within the specified search radius
			constexpr cell_index_list_t close_key_parts(const value_t         &prm_value,        ///< The value for which the key_part should be extracted
			                                            const search_radius_t &prm_search_radius ///< The search radius defining what is considered a match
			                                            ) const {
				return
#ifndef NDEBUG
					( prm_search_radius == 0 )
#else
					common::constexpr_ignore_unused( prm_search_radius )
#endif
						? cell_index_list_t{ { prm_value } }
						: throw std::invalid_argument("simple_locn_index_keyer_part currently requires that the search radius is 0 (ie requires matching indices)");
			}
		};


	} // namespace scan
} // namespace cath

#endif
