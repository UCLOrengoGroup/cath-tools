/// \file
/// \brief The indexed_refiner class header

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

#ifndef _CATH_TOOLS_SOURCE_INDEXED_REFINER_INDEXED_REFINER_H
#define _CATH_TOOLS_SOURCE_INDEXED_REFINER_INDEXED_REFINER_H

#include "alignment/align_type_aliases.hpp"
#include "alignment/pair_alignment.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "common/type_aliases.hpp"
#include "scan/detail/scan_index_store/scan_index_hash_store.hpp"
#include "scan/detail/scan_index_store/scan_index_lattice_store.hpp"
#include "scan/detail/scan_index_store/scan_index_store_helper.hpp"
#include "scan/detail/scan_index_store/scan_index_vector_store.hpp"
#include "scan/spatial_index/spatial_index.hpp"

// namespace cath { namespace index { class view_cache_list; } }
namespace cath { class protein_list; }
namespace cath { namespace align { class alignment; } }
// namespace cath { namespace align { namespace detail { class alignment_split; } } }
// namespace cath { namespace align { namespace detail { class alignment_split_list; } } }
// namespace cath { namespace align { namespace gap { class gap_penalty; } } }

namespace cath {
	namespace align {

		enum class fot {
			FROM,
			TO
		};

		// namespace detail {
		// 	template <typename Idx,
		// 	          typename Keyer>
		// 	Idx build_gubbins(const protein   &arg_protein,         ///< TODOCUMENT
		// 	                  const scan::sod &arg_sparse_or_dense, ///< TODOCUMENT
		// 	                  const fot       &arg_from_or_to,      ///< TODOCUMENT
		// 	                  const float     &arg_cell_size,       ///< TODOCUMENT
		// 	                  const Keyer     &arg_keyer            ///< TODOCUMENT
		// 	                  ) {
		// 		constexpr float MAX_DIST = 7.0; // ??????????????????????
		// 		constexpr scan::simple_locn_crit the_crit{ MAX_DIST * MAX_DIST };

		// 		const auto the_range = boost::irange( 0_z, arg_protein.get_length() )
		// 			| transformed( [&] (const size_t &x) {
		// 				return scan::make_simple_locn_index_of_ca( arg_protein.get_residue_ref_of_index( x ), debug_numeric_cast<unsigned int>( x ) );
		// 			} );

		// 		auto the_store = make_dense_empty_lattice_store( the_range, arg_keyer, the_crit );
		// 		for (const auto &value : the_range) {
		// 			const auto close_keys = arg_keyer.make_close_keys( value, the_crit );
		// 			for (const auto &the_key : common::cross( close_keys ) ) {
		// 				the_store.push_back_entry_to_cell( the_key, value );
		// 			}
		// 		}
		// 		return the_store;


		// 		return make_dense_lattice_store(
		// 			,
		// 			arg_keyer,
		// 			,
		// 			arg_cell_size,
		// 			MAX_DIST
		// 		);

		// 		// make_dense_lattice_store_impl(
		// 		// 	irange( 0_z, arg_protein.get_length() )
		// 		// 		| transformed( [&] (const size_t &x) {
		// 		// 			return make_simple_locn_index_of_ca( arg_protein.get_residue_ref_of_index( x ), debug_numeric_cast<unsigned int>( x ) );
		// 		// 		} ),
		// 		// 	arg_cell_size,
		// 		// 	arg_max_dist
		// 		// );
		// 	}
		// }

		/// \brief TODOCUMENT
		inline geom::coord view_vector_of_rotns_coords_and_indices(const std::vector<geom::rotation> &arg_rotns,      ///< TODOCUMENT
		                                                           const std::vector<geom::coord>    &arg_coords,     ///< TODOCUMENT
		                                                           const size_t                      &arg_from_index, ///< TODOCUMENT
		                                                           const size_t                      &arg_to_index    ///< TODOCUMENT
		                                                           ) {
			return rotate_copy(
				arg_rotns[ arg_from_index ],
				arg_coords[ arg_to_index ] - arg_coords[ arg_from_index ]
			);
		}

		/// \brief TODOCUMENT
		inline auto make_range_of_index(const std::vector<geom::rotation> &arg_rotns,      ///< TODOCUMENT
		                                const std::vector<geom::coord>    &arg_coords,     ///< TODOCUMENT
		                                const fot                         &arg_from_or_to, ///< TODOCUMENT
		                                const size_t                      &arg_index       ///< TODOCUMENT
		                                ) {
			return boost::irange( 0u, debug_numeric_cast<unsigned int>( arg_rotns.size() ) )
				| boost::adaptors::transformed( [&] (const unsigned int &x) {
					return scan::make_simple_locn_index(
						view_vector_of_rotns_coords_and_indices(
							arg_rotns,
							arg_coords,
							( ( arg_from_or_to == fot::FROM ) ? arg_index : x         ),
							( ( arg_from_or_to == fot::FROM ) ? x         : arg_index )
						),
						x
					);
				} );
		}

		/// \brief TODOCUMENT
		inline auto make_range_of_all(const std::vector<geom::rotation> &arg_rotns,     ///< TODOCUMENT
		                              const std::vector<geom::coord>    &arg_coords,    ///< TODOCUMENT
		                              const fot                         &arg_from_or_to ///< TODOCUMENT
		                              ) {
			const auto index_range = boost::irange( 0u, debug_numeric_cast<unsigned int>( arg_rotns.size() ) );
			return common::cross( index_range, index_range )
				| boost::adaptors::transformed( [&] (const std::tuple<unsigned int, unsigned int> &x) {
					return scan::make_simple_locn_index(
						view_vector_of_rotns_coords_and_indices(
							arg_rotns,
							arg_coords,
							( ( arg_from_or_to == fot::FROM ) ? std::get<0>( x ) : std::get<1>( x ) ),
							( ( arg_from_or_to == fot::FROM ) ? std::get<1>( x ) : std::get<0>( x ) )
						),
						std::get<1>( x )
					);
				} );
		}

		/// \brief TODOCUMENT
		inline auto make_rotns_and_coords(const protein &arg_protein ///< TODOCUMENT
		                                  ) {
			std::pair<std::vector<geom::rotation>, std::vector<geom::coord>> result;
			for (const size_t &x : boost::irange( 0_z, arg_protein.get_length() ) ) {
				result . first  . push_back( arg_protein.get_residue_ref_of_index( x ).get_frame()             );
				result . second . push_back( arg_protein.get_residue_ref_of_index( x ).get_carbon_beta_coord() );
			}
			return result;
		}

		/// \brief TODOCUMENT
		inline auto make_range_of_index(const protein &arg_protein,    ///< TODOCUMENT
		                                const fot     &arg_from_or_to, ///< TODOCUMENT
		                                const size_t  &arg_index       ///< TODOCUMENT
		                                ) {
			const auto rotns_and_coords = make_rotns_and_coords( arg_protein );
			return make_range_of_index(
				rotns_and_coords.first,
				rotns_and_coords.second,
				arg_from_or_to,
				arg_index
			);
		}

		/// \brief TODOCUMENT
		inline auto make_range_of_all(const protein &arg_protein,   ///< TODOCUMENT
		                              const fot     &arg_from_or_to ///< TODOCUMENT
		                              ) {
			const auto rotns_and_coords = make_rotns_and_coords( arg_protein );
			return make_range_of_all(
				rotns_and_coords.first,
				rotns_and_coords.second,
				arg_from_or_to
			);
		}

		// Expect Cell to be something like std::vector<scan::simple_locn_index> or a small-size-optimized version of that
		template <typename Cell>
		struct vector_refine_index final {

			/// \brief TODOCUMENT
			using store_type = scan::detail::scan_index_vector_store<std::tuple<int, int, int, int>, Cell>;

			/// \brief TODOCUMENT
			scan::res_pair_keyer<scan::simple_locn_x_keyer_part,
			                     scan::simple_locn_y_keyer_part,
			                     scan::simple_locn_z_keyer_part,
			                     scan::res_pair_from_to_index_keyer_part> the_keyer;

			/// \brief TODOCUMENT
			store_type the_store;

			/// \brief TODOCUMENT
			static constexpr scan::detail::view_base_type MAX_DIST = 7.0;

			/// \brief TODOCUMENT
			template <scan::sod Sod>
			vector_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
			                    const protein                                &arg_protein,    ///< TODOCUMENT
			                    const fot                                    &arg_from_or_to, ///< TODOCUMENT
			                    const float                                  &arg_cell_size   ///< TODOCUMENT
			                    ) : the_keyer{
			                        	make_res_pair_keyer(
			                        		scan::simple_locn_x_keyer_part{ arg_cell_size },
			                        		scan::simple_locn_y_keyer_part{ arg_cell_size },
			                        		scan::simple_locn_z_keyer_part{ arg_cell_size },
			                        		scan::res_pair_from_to_index_keyer_part{}
			                        	)
			                        },
			                        the_store{ scan::detail::store_maker<Sod, store_type>{}(
			                        	make_range_of_all( arg_protein, arg_from_or_to ),
			                        	the_keyer,
			                        	scan::simple_locn_crit{ MAX_DIST * MAX_DIST }
			                        ) } {
			}

			/// \brief TODOCUMENT
			scan::info_quantity get_info_size() const {
				return the_store.get_info_size();
			}
		};

		/// \brief TODOCUMENT
		template <typename Cell>
		struct vec_of_vectors_refine_index final {

			/// \brief TODOCUMENT
			using store_type = scan::detail::scan_index_vector_store<std::tuple<int, int, int>, Cell>;

			/// \brief TODOCUMENT
			scan::res_pair_keyer<scan::simple_locn_x_keyer_part,
			                     scan::simple_locn_y_keyer_part,
			                     scan::simple_locn_z_keyer_part,
			                     scan::res_pair_from_to_index_keyer_part> the_keyer;

			/// \brief TODOCUMENT
			std::vector<store_type> the_store;

			/// \brief TODOCUMENT
			static constexpr scan::detail::view_base_type MAX_DIST = 7.0;

			/// \brief TODOCUMENT
			template <scan::sod Sod>
			vec_of_vectors_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
			                            const protein                                &arg_protein,    ///< TODOCUMENT
			                            const fot                                    &arg_from_or_to, ///< TODOCUMENT
			                            const float                                  &arg_cell_size   ///< TODOCUMENT
			                            ) : the_keyer{
			                                	make_res_pair_keyer(
			                                		scan::simple_locn_x_keyer_part{ arg_cell_size },
			                                		scan::simple_locn_y_keyer_part{ arg_cell_size },
			                                		scan::simple_locn_z_keyer_part{ arg_cell_size },
			                                		scan::res_pair_from_to_index_keyer_part{}
			                                	)
			                                },
			                                the_store{ common::transform_build<std::vector<store_type>>(
			                                	boost::irange( 0_z, arg_protein.get_length() ),
			                                	[&] (const size_t &x) {
			                                		return scan::detail::store_maker<Sod, store_type>{}(
			                                			make_range_of_index(
			                                				arg_protein,
			                                				arg_from_or_to,
			                                				x
			                                			),
			                                			the_keyer,
			                                			scan::simple_locn_crit{ MAX_DIST * MAX_DIST }
			                                		);
			                                	}
			                                ) } {
			}

			/// \brief TODOCUMENT
			scan::info_quantity get_info_size() const {
				return accumulate(
					the_store
						| boost::adaptors::transformed( [] (const store_type &x) {
							return x.get_info_size();
						} ),
					scan::info_quantity{ 0 }
				)
				+ ( sizeof( decltype( the_store ) ) * boost::units::information::bytes );
			}

			void store_sparse(const int &,
			                  const scan::simple_locn_index &);
		};

		template <typename Cell>
		struct hash_refine_index final {

			/// \brief TODOCUMENT
			using store_type = scan::detail::scan_index_hash_store<std::tuple<int, int, int, int>, Cell>;

			/// \brief TODOCUMENT
			scan::res_pair_keyer<scan::simple_locn_x_keyer_part,
			                     scan::simple_locn_y_keyer_part,
			                     scan::simple_locn_z_keyer_part,
			                     scan::res_pair_from_to_index_keyer_part> the_keyer;

			/// \brief TODOCUMENT
			store_type the_store;

			/// \brief TODOCUMENT
			static constexpr scan::detail::view_base_type MAX_DIST = 7.0;

			/// \brief TODOCUMENT
			template <scan::sod Sod>
			hash_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
			                  const protein                                &arg_protein,    ///< TODOCUMENT
			                  const fot                                    &arg_from_or_to, ///< TODOCUMENT
			                  const float                                  &arg_cell_size   ///< TODOCUMENT
			                  ) : the_keyer{
			                      	make_res_pair_keyer(
			                      		scan::simple_locn_x_keyer_part{ arg_cell_size },
			                      		scan::simple_locn_y_keyer_part{ arg_cell_size },
			                      		scan::simple_locn_z_keyer_part{ arg_cell_size },
			                      		scan::res_pair_from_to_index_keyer_part{}
			                      	)
			                      },
			                      the_store{ scan::detail::store_maker<Sod, store_type>{}(
			                      	make_range_of_all( arg_protein, arg_from_or_to ),
			                      	the_keyer,
			                      	scan::simple_locn_crit{ MAX_DIST * MAX_DIST }
			                      ) } {
			}

			/// \brief TODOCUMENT
			scan::info_quantity get_info_size() const {
				return the_store.get_info_size();
			}

			void store_sparse(const int &,
			                  const scan::simple_locn_index &);
		};

		template <typename Cell>
		struct vec_of_hashes_refine_index final {

			/// \brief TODOCUMENT
			using store_type = scan::detail::scan_index_hash_store<std::tuple<int, int, int>, Cell>;

			/// \brief TODOCUMENT
			scan::res_pair_keyer<scan::simple_locn_x_keyer_part,
			                     scan::simple_locn_y_keyer_part,
			                     scan::simple_locn_z_keyer_part,
			                     scan::res_pair_from_to_index_keyer_part> the_keyer;

			/// \brief TODOCUMENT
			std::vector< store_type > the_store;

			/// \brief TODOCUMENT
			static constexpr scan::detail::view_base_type MAX_DIST = 7.0;

			/// \brief TODOCUMENT
			template <scan::sod Sod>
			vec_of_hashes_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
			                           const protein                                &arg_protein,    ///< TODOCUMENT
			                           const fot                                    &arg_from_or_to, ///< TODOCUMENT
			                           const float                                  &arg_cell_size   ///< TODOCUMENT
			                           ) : the_keyer{
			                               	make_res_pair_keyer(
			                               		scan::simple_locn_x_keyer_part{ arg_cell_size },
			                               		scan::simple_locn_y_keyer_part{ arg_cell_size },
			                               		scan::simple_locn_z_keyer_part{ arg_cell_size },
			                               		scan::res_pair_from_to_index_keyer_part{}
			                               	)
			                               },
			                               the_store{ common::transform_build<std::vector<store_type>>(
			                               	boost::irange( 0_z, arg_protein.get_length() ),
			                               	[&] (const size_t &x) {
			                               		return scan::detail::store_maker<Sod, store_type>{}(
			                               			make_range_of_index(
			                               				arg_protein,
			                               				arg_from_or_to,
			                               				x
			                               			),
			                               			the_keyer,
			                               			scan::simple_locn_crit{ MAX_DIST * MAX_DIST }
			                               		);
			                               	}
			                               ) } {
			}

			/// \brief TODOCUMENT
			scan::info_quantity get_info_size() const {
				return accumulate(
					the_store
						| boost::adaptors::transformed( [] (const store_type &x) {
							return x.get_info_size();
						} ),
					scan::info_quantity{ 0 }
				)
				+ ( sizeof( decltype( the_store ) ) * boost::units::information::bytes );
			}

			void store_sparse(const int &,
			                  const scan::simple_locn_index &);
		};

		template <typename Cell>
		struct lattice_refine_index final {

			/// \brief TODOCUMENT
			using store_type = scan::detail::scan_index_lattice_store<std::tuple<int, int, int, int>, Cell>;

			/// \brief TODOCUMENT
			scan::res_pair_keyer<scan::simple_locn_x_keyer_part,
			                     scan::simple_locn_y_keyer_part,
			                     scan::simple_locn_z_keyer_part,
			                     scan::res_pair_from_to_index_keyer_part> the_keyer;

			/// \brief TODOCUMENT
			store_type the_store;

			/// \brief TODOCUMENT
			static constexpr scan::detail::view_base_type MAX_DIST = 7.0;

			/// \brief TODOCUMENT
			template <scan::sod Sod>
			lattice_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
			                     const protein                                &arg_protein,    ///< TODOCUMENT
			                     const fot                                    &arg_from_or_to, ///< TODOCUMENT
			                     const float                                  &arg_cell_size   ///< TODOCUMENT
			                     ) : the_keyer{
			                         	make_res_pair_keyer(
			                         		scan::simple_locn_x_keyer_part{ arg_cell_size },
			                         		scan::simple_locn_y_keyer_part{ arg_cell_size },
			                         		scan::simple_locn_z_keyer_part{ arg_cell_size },
			                         		scan::res_pair_from_to_index_keyer_part{}
			                         	)
			                         },
			                         the_store{ scan::detail::store_maker<Sod, store_type>{}(
			                         	make_range_of_all( arg_protein, arg_from_or_to ),
			                         	the_keyer,
			                         	scan::simple_locn_crit{ MAX_DIST * MAX_DIST }
			                         ) } {
				// if ( arg_sparse_or_dense == scan::sod::SPARSE ) {
				// 	// locn_index_store cath::scan::
				// 	make_sparse_lattice(
				// 		arg_protein,
				// 		arg_cell_size,
				// 		arg_max_dist
				// 	);
				// }
				// else {
				// 	// locn_index_store cath::scan::
				// 	make_dense_lattice(const protein &arg_protein,   ///< TODOCUMENT
				// 	                   const float   &arg_cell_size, ///< TODOCUMENT
				// 	                   const float   &arg_max_dist   ///< TODOCUMENT
				// 	                   );
				// }
			}

			/// \brief TODOCUMENT
			scan::info_quantity get_info_size() const {
				return the_store.get_info_size();
			}

			void store_sparse(const int &,
			                  const scan::simple_locn_index &);
		};

		template <typename Cell>
		struct vec_of_lattices_refine_index final {

			/// \brief TODOCUMENT
			using store_type = scan::detail::scan_index_lattice_store<std::tuple<int, int, int>, Cell>;

			/// \brief TODOCUMENT
			scan::res_pair_keyer<scan::simple_locn_x_keyer_part,
			                     scan::simple_locn_y_keyer_part,
			                     scan::simple_locn_z_keyer_part,
			                     scan::res_pair_from_to_index_keyer_part> the_keyer;

			/// \brief TODOCUMENT
			std::vector< store_type > the_store;

			/// \brief TODOCUMENT
			static constexpr scan::detail::view_base_type MAX_DIST = 7.0;

			/// \brief TODOCUMENT
			template <scan::sod Sod>
			vec_of_lattices_refine_index(const std::integral_constant<scan::sod, Sod> &, ///< TODOCUMENT
			                             const protein &arg_protein,                     ///< TODOCUMENT
			                             const fot     &arg_from_or_to,                  ///< TODOCUMENT
			                             const float   &arg_cell_size                    ///< TODOCUMENT
			                             ) : the_keyer{
			                                 	make_res_pair_keyer(
			                                 		scan::simple_locn_x_keyer_part{ arg_cell_size },
			                                 		scan::simple_locn_y_keyer_part{ arg_cell_size },
			                                 		scan::simple_locn_z_keyer_part{ arg_cell_size },
			                                 		scan::res_pair_from_to_index_keyer_part{}
			                                 	)
			                                 },
			                                 the_store{ common::transform_build<std::vector<store_type>>(
			                                 	boost::irange( 0_z, arg_protein.get_length() ),
			                                 	[&] (const size_t &x) {
			                                 		return scan::detail::store_maker<Sod, store_type>{}(
			                                 			make_range_of_index(
			                                 				arg_protein,
			                                 				arg_from_or_to,
			                                 				x
			                                 			),
			                                 			the_keyer,
			                                 			scan::simple_locn_crit{ MAX_DIST * MAX_DIST }
			                                 		);
			                                 	}
			                                 ) } {
			}

			/// \brief TODOCUMENT
			scan::info_quantity get_info_size() const {
				return accumulate(
					the_store
						| boost::adaptors::transformed( [] (const store_type &x) {
							return x.get_info_size();
						} ),
					scan::info_quantity{ 0 }
				)
				+ ( sizeof( decltype( the_store ) ) * boost::units::information::bytes );
			}
		};

		// /// \brief TODOCUMENT
		// class indexed_refiner final {
		// private:
		// 	/// \brief TODOCUMENT
		// 	float_score_vec_vec from_and_to_alignment_scores;

		// 	/// \brief TODOCUMENT
		// 	scan::detail::scan_index_lattice_store< std::tuple<int, int, int, int>, std::vector<scan::simple_locn_index> > barry;

		// 	// /// \brief TODOCUMENT
		// 	// float_score_vec_vec to_alignment_scores;

		// 	detail::bool_aln_pair iterate_step(const alignment &,
		// 	                                   const protein_list &,
		// 	                                   // const index::view_cache_list &,
		// 	                                   const gap::gap_penalty &);

		// 	detail::bool_aln_pair iterate_step_for_alignment_split_list(const alignment &,
		// 	                                                            const protein_list &,
		// 	                                                            // const index::view_cache_list &,
		// 	                                                            const gap::gap_penalty &,
		// 	                                                            const detail::alignment_split_list &);

		// 	detail::bool_aln_pair iterate_step_for_alignment_split(const alignment &,
		// 	                                                       const protein_list &,
		// 	                                                       // const index::view_cache_list &,
		// 	                                                       const gap::gap_penalty &,
		// 	                                                       const detail::alignment_split &);

		// public:
		// 	alignment iterate(const alignment &,
		// 	                  const protein_list &,
		// 	                  const gap::gap_penalty &);

		// };

		void do_some_gubbins(const protein &,
		                     const protein &,
		                     const alignment &);

	} // namespace align

} // namespace cath

#endif
