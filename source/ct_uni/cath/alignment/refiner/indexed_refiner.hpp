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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_INDEXED_REFINER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_INDEXED_REFINER_HPP

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/alignment/pair_alignment.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_hash_store.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_lattice_store.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_store_helper.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_vector_store.hpp"
#include "cath/scan/spatial_index/spatial_index.hpp"

// clang-format off
namespace cath { class protein_list; }
namespace cath::align { class alignment; }
// clang-format on

namespace cath::align {

	enum class fot : bool {
		FROM,
		TO
	};

	struct indexed_refiner_constants final {
		static constexpr float MAX_DIST    = 7.0;
		static constexpr float MAX_SQ_DIST = MAX_DIST * MAX_DIST;
	};

	// namespace detail {
	// 	template <typename Idx,
	// 	          typename Keyer>
	// 	Idx build_gubbins(const protein   &prm_protein,         ///< TODOCUMENT
	// 	                  const scan::sod &prm_sparse_or_dense, ///< TODOCUMENT
	// 	                  const fot       &prm_from_or_to,      ///< TODOCUMENT
	// 	                  const float     &prm_cell_size,       ///< TODOCUMENT
	// 	                  const Keyer     &prm_keyer            ///< TODOCUMENT
	// 	                  ) {
	// 		constexpr float MAX_DIST = 7.0; // ??????????????????????
	// 		constexpr scan::simple_locn_crit the_crit{ indexed_refiner_constants::MAX_SQ_DIST };

	// 		const auto the_range = common::indices( prm_protein.get_length() )
	// 			| transformed( [&] (const size_t &x) {
	// 				return scan::make_simple_locn_index_of_ca( prm_protein.get_residue_ref_of_index( x ), debug_numeric_cast<unsigned int>( x ) );
	// 			} );

	// 		auto the_store = make_dense_empty_lattice_store( the_range, prm_keyer, the_crit );
	// 		for (const auto &value : the_range) {
	// 			const auto close_keys = prm_keyer.make_close_keys( value, the_crit );
	// 			for (const auto &the_key : common::cross( close_keys ) ) {
	// 				the_store.push_back_entry_to_cell( the_key, value );
	// 			}
	// 		}
	// 		return the_store;


	// 		return make_dense_lattice_store(
	// 			,
	// 			prm_keyer,
	// 			,
	// 			prm_cell_size,
	// 			MAX_DIST
	// 		);

	// 		// make_dense_lattice_store_impl(
	// 		// 	common::indices( prm_protein.get_length() )
	// 		// 		| transformed( [&] (const size_t &x) {
	// 		// 			return make_simple_locn_index_of_ca( prm_protein.get_residue_ref_of_index( x ), debug_numeric_cast<unsigned int>( x ) );
	// 		// 		} ),
	// 		// 	prm_cell_size,
	// 		// 	prm_max_dist
	// 		// );
	// 	}
	// }

	/// \brief TODOCUMENT
	inline geom::coord view_vector_of_rotns_coords_and_indices(const std::vector<geom::rotation> &prm_rotns,      ///< TODOCUMENT
	                                                           const std::vector<geom::coord>    &prm_coords,     ///< TODOCUMENT
	                                                           const size_t                      &prm_from_index, ///< TODOCUMENT
	                                                           const size_t                      &prm_to_index    ///< TODOCUMENT
	                                                           ) {
		return rotate_copy(
			prm_rotns[ prm_from_index ],
			prm_coords[ prm_to_index ] - prm_coords[ prm_from_index ]
		);
	}

	/// \brief TODOCUMENT
	///
	/// Note that the range this returns uses references to all of the arguments so it
	/// mustn't be used beyond the end of any of their lifetimes.
	inline auto make_range_of_index(const std::vector<geom::rotation> &prm_rotns,      ///< TODOCUMENT
	                                const std::vector<geom::coord>    &prm_coords,     ///< TODOCUMENT
	                                const fot                         &prm_from_or_to, ///< TODOCUMENT
	                                const size_t                      &prm_index       ///< TODOCUMENT
	                                ) {
		return common::indices( debug_numeric_cast<unsigned int>( prm_rotns.size() ) )
			| boost::adaptors::transformed( [&] (const unsigned int &x) {
				return scan::make_simple_locn_index(
					view_vector_of_rotns_coords_and_indices(
						prm_rotns,
						prm_coords,
						( ( prm_from_or_to == fot::FROM ) ? prm_index : x         ),
						( ( prm_from_or_to == fot::FROM ) ? x         : prm_index )
					),
					x
				);
			} );
	}

	/// \brief TODOCUMENT
	inline auto make_range_of_all(const std::vector<geom::rotation> &prm_rotns,     ///< TODOCUMENT
	                              const std::vector<geom::coord>    &prm_coords,    ///< TODOCUMENT
	                              const fot                         &prm_from_or_to ///< TODOCUMENT
	                              ) {
		const auto index_range = common::indices( debug_numeric_cast<unsigned int>( prm_rotns.size() ) );
		return common::cross( index_range, index_range )
			| boost::adaptors::transformed( [&] (const std::tuple<unsigned int, unsigned int> &x) {
				return scan::make_simple_locn_index(
					view_vector_of_rotns_coords_and_indices(
						prm_rotns,
						prm_coords,
						( ( prm_from_or_to == fot::FROM ) ? std::get<0>( x ) : std::get<1>( x ) ),
						( ( prm_from_or_to == fot::FROM ) ? std::get<1>( x ) : std::get<0>( x ) )
					),
					std::get<1>( x )
				);
			} );
	}

	/// \brief TODOCUMENT
	inline auto make_rotns_and_coords(const protein &prm_protein ///< TODOCUMENT
	                                  ) {
		std::pair<std::vector<geom::rotation>, std::vector<geom::coord>> result;
		for (const size_t &x : common::indices( prm_protein.get_length() ) ) {
			result . first  . push_back( prm_protein.get_residue_ref_of_index( x ).get_frame()             );
			result . second . push_back( prm_protein.get_residue_ref_of_index( x ).get_carbon_beta_coord() );
		}
		return result;
	}

	/// \brief TODOCUMENT
	inline auto make_range_of_all(const protein &prm_protein,   ///< TODOCUMENT
	                              const fot     &prm_from_or_to ///< TODOCUMENT
	                              ) {
		const auto rotns_and_coords = make_rotns_and_coords( prm_protein );
		return make_range_of_all(
			rotns_and_coords.first,
			rotns_and_coords.second,
			prm_from_or_to
		);
	}

	/// \brief TODOCUMENT
	template <scan::sod Sod,
	          typename T,
	          typename Keyer>
	std::vector<T> prepare_vec_based_simple_locn_index_store(const protein &prm_protein,    ///< TODOCUMENT
	                                                         const fot     &prm_from_or_to, ///< TODOCUMENT
	                                                         const Keyer   &prm_keyer       ///< TODOCUMENT
	                                                         ) {
		const auto    rotns_and_coords = make_rotns_and_coords( prm_protein );
		const size_t &length           = prm_protein.get_length();

		std::vector<T> result;
		result.reserve( length );

		std::vector<scan::simple_locn_index> lists;
		lists.reserve( length );

		for (const size_t &outer_index : common::indices( length ) ) {
			lists.clear();
			for (const unsigned int &inner_index : common::indices( debug_numeric_cast<unsigned int>( length ) ) ) {
				const auto the_view = view_vector_of_rotns_coords_and_indices(
					rotns_and_coords.first,
					rotns_and_coords.second,
					( ( prm_from_or_to == fot::FROM ) ? outer_index : inner_index ),
					( ( prm_from_or_to == fot::FROM ) ? inner_index : outer_index )
				);
				lists.emplace_back(
					debug_numeric_cast< scan::detail::view_base_type>( the_view.get_x() ),
					debug_numeric_cast< scan::detail::view_base_type>( the_view.get_y() ),
					debug_numeric_cast< scan::detail::view_base_type>( the_view.get_z() ),
					inner_index
				);
			}

			result.push_back(
				scan::detail::store_maker<Sod, T>{}(
					lists,
					prm_keyer,
					scan::simple_locn_crit{ indexed_refiner_constants::MAX_SQ_DIST }
				)
			);
		}

		return result;
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
		template <scan::sod Sod>
		vector_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
		                    const protein                                &prm_protein,    ///< TODOCUMENT
		                    const fot                                    &prm_from_or_to, ///< TODOCUMENT
		                    const float                                  &prm_cell_size   ///< TODOCUMENT
		                    ) : the_keyer{
		                        	make_res_pair_keyer(
		                        		scan::simple_locn_x_keyer_part{ prm_cell_size },
		                        		scan::simple_locn_y_keyer_part{ prm_cell_size },
		                        		scan::simple_locn_z_keyer_part{ prm_cell_size },
		                        		scan::res_pair_from_to_index_keyer_part{}
		                        	)
		                        },
		                        the_store{ scan::detail::store_maker<Sod, store_type>{}(
		                        	make_range_of_all( prm_protein, prm_from_or_to ),
		                        	the_keyer,
		                        	scan::simple_locn_crit{ indexed_refiner_constants::MAX_SQ_DIST }
		                        ) } {
		}

		/// \brief TODOCUMENT
		[[nodiscard]] scan::info_quantity get_info_size() const {
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
		template <scan::sod Sod>
		vec_of_vectors_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
		                            const protein                                &prm_protein,    ///< TODOCUMENT
		                            const fot                                    &prm_from_or_to, ///< TODOCUMENT
		                            const float                                  &prm_cell_size   ///< TODOCUMENT
		                            ) : the_keyer{
		                                	make_res_pair_keyer(
		                                		scan::simple_locn_x_keyer_part{ prm_cell_size },
		                                		scan::simple_locn_y_keyer_part{ prm_cell_size },
		                                		scan::simple_locn_z_keyer_part{ prm_cell_size },
		                                		scan::res_pair_from_to_index_keyer_part{}
		                                	)
		                                },
		                                the_store{
		                                	prepare_vec_based_simple_locn_index_store<Sod, store_type>(
		                                		prm_protein,
		                                		prm_from_or_to,
		                                		the_keyer
		                                	)
		                                } {
		}

		/// \brief TODOCUMENT
		[[nodiscard]] scan::info_quantity get_info_size() const {
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
		template <scan::sod Sod>
		hash_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
		                  const protein                                &prm_protein,    ///< TODOCUMENT
		                  const fot                                    &prm_from_or_to, ///< TODOCUMENT
		                  const float                                  &prm_cell_size   ///< TODOCUMENT
		                  ) : the_keyer{
		                      	make_res_pair_keyer(
		                      		scan::simple_locn_x_keyer_part{ prm_cell_size },
		                      		scan::simple_locn_y_keyer_part{ prm_cell_size },
		                      		scan::simple_locn_z_keyer_part{ prm_cell_size },
		                      		scan::res_pair_from_to_index_keyer_part{}
		                      	)
		                      },
		                      the_store{ scan::detail::store_maker<Sod, store_type>{}(
		                      	make_range_of_all( prm_protein, prm_from_or_to ),
		                      	the_keyer,
		                      	scan::simple_locn_crit{ indexed_refiner_constants::MAX_SQ_DIST }
		                      ) } {
		}

		/// \brief TODOCUMENT
		[[nodiscard]] scan::info_quantity get_info_size() const {
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
		template <scan::sod Sod>
		vec_of_hashes_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
		                           const protein                                &prm_protein,    ///< TODOCUMENT
		                           const fot                                    &prm_from_or_to, ///< TODOCUMENT
		                           const float                                  &prm_cell_size   ///< TODOCUMENT
		                           ) : the_keyer{
		                               	make_res_pair_keyer(
		                               		scan::simple_locn_x_keyer_part{ prm_cell_size },
		                               		scan::simple_locn_y_keyer_part{ prm_cell_size },
		                               		scan::simple_locn_z_keyer_part{ prm_cell_size },
		                               		scan::res_pair_from_to_index_keyer_part{}
		                               	)
		                               },
		                               the_store{
		                               	prepare_vec_based_simple_locn_index_store<Sod, store_type>(
		                               		prm_protein,
		                               		prm_from_or_to,
		                               		the_keyer
		                               	)
		                               } {
		}

		/// \brief TODOCUMENT
		[[nodiscard]] scan::info_quantity get_info_size() const {
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
		template <scan::sod Sod>
		lattice_refine_index(const std::integral_constant<scan::sod, Sod> &,               ///< TODOCUMENT
		                     const protein                                &prm_protein,    ///< TODOCUMENT
		                     const fot                                    &prm_from_or_to, ///< TODOCUMENT
		                     const float                                  &prm_cell_size   ///< TODOCUMENT
		                     ) : the_keyer{
		                         	make_res_pair_keyer(
		                         		scan::simple_locn_x_keyer_part{ prm_cell_size },
		                         		scan::simple_locn_y_keyer_part{ prm_cell_size },
		                         		scan::simple_locn_z_keyer_part{ prm_cell_size },
		                         		scan::res_pair_from_to_index_keyer_part{}
		                         	)
		                         },
		                         the_store{ scan::detail::store_maker<Sod, store_type>{}(
		                         	make_range_of_all( prm_protein, prm_from_or_to ),
		                         	the_keyer,
		                         	scan::simple_locn_crit{ indexed_refiner_constants::MAX_SQ_DIST }
		                         ) } {
			// if ( prm_sparse_or_dense == scan::sod::SPARSE ) {
			// 	// locn_index_store cath::scan::
			// 	make_sparse_lattice(
			// 		prm_protein,
			// 		prm_cell_size,
			// 		prm_max_dist
			// 	);
			// }
			// else {
			// 	// locn_index_store cath::scan::
			// 	make_dense_lattice(const protein &prm_protein,   ///< TODOCUMENT
			// 	                   const float   &prm_cell_size, ///< TODOCUMENT
			// 	                   const float   &prm_max_dist   ///< TODOCUMENT
			// 	                   );
			// }
		}

		/// \brief TODOCUMENT
		[[nodiscard]] scan::info_quantity get_info_size() const {
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
		template <scan::sod Sod>
		vec_of_lattices_refine_index(const std::integral_constant<scan::sod, Sod> &, ///< TODOCUMENT
		                             const protein &prm_protein,                     ///< TODOCUMENT
		                             const fot     &prm_from_or_to,                  ///< TODOCUMENT
		                             const float   &prm_cell_size                    ///< TODOCUMENT
		                             ) : the_keyer{
		                                 	make_res_pair_keyer(
		                                 		scan::simple_locn_x_keyer_part{ prm_cell_size },
		                                 		scan::simple_locn_y_keyer_part{ prm_cell_size },
		                                 		scan::simple_locn_z_keyer_part{ prm_cell_size },
		                                 		scan::res_pair_from_to_index_keyer_part{}
		                                 	)
		                                 },
		                                 the_store{
		                                 	prepare_vec_based_simple_locn_index_store<Sod, store_type>(
		                                 		prm_protein,
		                                 		prm_from_or_to,
		                                 		the_keyer
		                                 	)
		                                 } {
		}

		/// \brief TODOCUMENT
		[[nodiscard]] scan::info_quantity get_info_size() const {
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

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_REFINER_INDEXED_REFINER_HPP
