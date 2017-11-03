/// \file
/// \brief The scan_index_store helper header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_STORE_HELPER_H
#define _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_STORE_HELPER_H

#include "common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "scan/detail/scan_index_store/scan_index_lattice_store.hpp"
#include "scan/detail/scan_role.hpp"
#include "scan/detail/scan_type_aliases.hpp"
#include "scan/detail/stride/roled_scan_stride.hpp"
#include "scan/scan_policy.hpp"
#include "structure/protein/protein.hpp"

using namespace std::literals::string_literals;

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		enum class sod : bool {
			SPARSE,
			DENSE
		};

		/// \brief TODOCUMENT
		inline std::string to_string(const sod &arg_sod ///< TODOCUMENT
		                             ) {
			switch ( arg_sod ) {
				case ( sod::SPARSE ) : { return "sparse"s; }
				case ( sod::DENSE  ) : { return "dense"s;  }
			}
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of sod not recognised whilst converting to_string()"));
		}

		namespace detail {

			/// \brief TODOCUMENT
			template <sod Sod,
			          typename Store>
			struct empty_store_maker final {
				template <typename Rng,
				          typename Keyer,
				          typename Crit>
				Store operator()(const Rng   &/*arg_rng*/,   ///< TODOCUMENT
				                 const Keyer &/*arg_keyer*/, ///< TODOCUMENT
				                 const Crit  &/*arg_crit*/   ///< TODOCUMENT
				                 ) {
					return Store{};
				}
			};

			// template <typename T> class TD;

			/// \brief TODOCUMENT
			template <typename Key,
			          typename Cell>
			struct empty_store_maker<sod::SPARSE, scan_index_lattice_store<Key, Cell>> final {

				template <typename Rng,
				          typename Keyer,
				          typename Crit>
				auto operator()(const Rng   &arg_rng,     ///< TODOCUMENT
				                const Keyer &arg_keyer,   ///< TODOCUMENT
				                const Crit  &/*arg_crit*/ ///< TODOCUMENT
				                ) {
					using value_t = common::range_value_t< Rng >;

					if ( boost::empty( arg_rng ) ) {
						BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot make a scan_index_store for an empty range"));
					}

					const auto mins_maxs = common::tuple_mins_maxs_element(
						arg_rng
							| boost::adaptors::transformed( [&] (const value_t &value) { return arg_keyer.make_key( value ); } )
					);

					// TD< value_t > value_t_obj;
					// TD< decltype( mins_maxs.first  ) > first_obj;

					// static_assert( std::is_same< decltype( mins_maxs.first  ), value_t >::value, "" );
					// static_assert( std::is_same< decltype( mins_maxs.second ), value_t >::value, "" );

					return scan_index_lattice_store<Key, Cell>( mins_maxs.first, mins_maxs.second );
				}
			};

			/// \brief TODOCUMENT
			template <typename Key,
			          typename Cell>
			struct empty_store_maker<sod::DENSE, scan_index_lattice_store<Key, Cell>> final {

				template <typename Rng,
				          typename Keyer,
				          typename Crit>
				auto operator()(const Rng   &arg_rng,   ///< TODOCUMENT
				                const Keyer &arg_keyer, ///< TODOCUMENT
				                const Crit  &arg_crit   ///< TODOCUMENT
				                ) {
					using value_t = common::range_value_t< Rng >;

					if ( boost::empty( arg_rng ) ) {
						BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot make a scan_index_store for an empty range"));
					}
					const auto mins_maxs = common::mins_maxs_tuple_pair_mins_maxs_element(
						arg_rng
							| boost::adaptors::transformed( [&] (const value_t &value) {
								return std::make_pair(
									arg_keyer.make_min_close_key( value, arg_crit ),
									arg_keyer.make_max_close_key( value, arg_crit )
								);
							} )
					);

					// static_assert(std::is_same< decltype( mins_maxs.first  ), value_t >::value, "" );
					// static_assert(std::is_same< decltype( mins_maxs.second ), value_t >::value, "" );

					return scan_index_lattice_store<Key, Cell>( mins_maxs.first, mins_maxs.second );
				}
			};



			/// \brief TODOCUMENT
			template <sod Sod,
			          typename Store>
			struct store_maker final {};

			/// \brief TODOCUMENT
			template <typename Store>
			struct store_maker<sod::SPARSE, Store> final {

				template <typename Rng,
				          typename Keyer,
				          typename Crit>
				/// \brief TODOCUMENT
				auto operator()(const Rng   &arg_rng,   ///< TODOCUMENT
				                const Keyer &arg_keyer, ///< TODOCUMENT
				                const Crit  &arg_crit   ///< TODOCUMENT
				                ) {
					Store the_store = empty_store_maker<sod::SPARSE, Store>{}( arg_rng, arg_keyer, arg_crit );
					for (const auto &data : arg_rng) {
						arg_keyer.store_emplace_value(
							the_store,
							arg_keyer.make_key( data ),
							data
						);
					}
					return the_store;
				}
			};

			/// \brief TODOCUMENT
			template <typename Store>
			struct store_maker<sod::DENSE, Store> final {

				/// \brief TODOCUMENT
				template <typename Rng,
				          typename Keyer,
				          typename Crit>
				auto operator()(const Rng   &arg_rng,   ///< TODOCUMENT
				                const Keyer &arg_keyer, ///< TODOCUMENT
				                const Crit  &arg_crit   ///< TODOCUMENT
				                ) {
					Store the_store = empty_store_maker<sod::DENSE, Store>{}( arg_rng, arg_keyer, arg_crit );
					for (const auto &data : arg_rng) {
						const auto close_keys = arg_keyer.make_close_keys( data, arg_crit );
						for (const auto &the_key : common::cross( close_keys ) ) {
							arg_keyer.store_emplace_value(
								the_store,
								the_key,
								data
							);
						}
					}
					return the_store;
				}
			};

			/// \brief TODOCUMENT
			///
			/// \todo Eliminate redundancy between add_structure_to_store() and dense_add_structure_to_store()
			template <typename T, typename... KPs>
			void add_structure_to_store(T                         &arg_store,           ///< TODOCUMENT
			                            const index_type          &arg_structure_index, ///< TODOCUMENT
			                            const protein             &arg_protein,         ///< TODOCUMENT
			                            const scan_policy<KPs...> &arg_scan_policy,     ///< TODOCUMENT
			                            const scan_role           &arg_scan_role        ///< TODOCUMENT
			                            ) {
				const auto &the_keyer        = arg_scan_policy.get_keyer();
//				BOOST_LOG_TRIVIAL( warning ) << "Keyer is : " << the_keyer;
				const auto  roled_stride     = roled_scan_stride{ arg_scan_role, arg_scan_policy.get_scan_stride() };
				const auto  num_residues     = debug_unwarned_numeric_cast<index_type>( arg_protein.get_length() );
				const auto  from_rep_strider = get_this_from_strider( roled_stride );
				const auto  to_rep_strider   = get_this_to_strider  ( roled_stride );
				for (const auto &from_rep_index : get_rep_indices_range( from_rep_strider, num_residues ) ) {
//					 BOOST_LOG_TRIVIAL( warning ) << "\tFrom rep " << from_rep_index;
					for (const auto &to_rep_index : get_rep_indices_range( to_rep_strider, num_residues ) ) {

						if ( from_rep_index != to_rep_index ) {
							const auto the_res_pair = make_multi_struc_res_rep_pair(
								arg_protein.get_residue_ref_of_index( get_index_of_rep_index( from_rep_strider, from_rep_index ) ),
								arg_protein.get_residue_ref_of_index( get_index_of_rep_index( to_rep_strider,   to_rep_index   ) ),
								arg_structure_index,
								from_rep_index,
								to_rep_index
							);
//							BOOST_LOG_TRIVIAL( warning ) << "\t\tTo rep " << to_rep_index << " - the rep : " << the_res_pair;
							arg_store.push_back_entry_to_cell( the_keyer.make_key( the_res_pair ), the_res_pair );
						}
					}
				}
			}

			/// \brief TODOCUMENT
			///
			/// \todo Eliminate redundancy between add_structure_to_store() and dense_add_structure_to_store()
			template <typename T, typename... KPs>
			void dense_add_structure_to_store(T                         &arg_store,           ///< TODOCUMENT
			                                  const index_type          &arg_structure_index, ///< TODOCUMENT
			                                  const protein             &arg_protein,         ///< TODOCUMENT
			                                  const scan_policy<KPs...> &arg_scan_policy,     ///< TODOCUMENT
			                                  const scan_role           &arg_scan_role        ///< TODOCUMENT
			                                  ) {
				const auto &the_criteria     = arg_scan_policy.get_criteria();
				const auto &the_keyer        = arg_scan_policy.get_keyer();
//				BOOST_LOG_TRIVIAL( warning ) << "Keyer is : " << the_keyer;
				const auto  roled_stride     = roled_scan_stride{ arg_scan_role, arg_scan_policy.get_scan_stride() };
				const auto  num_residues     = debug_unwarned_numeric_cast<index_type>( arg_protein.get_length() );
				const auto  from_rep_strider = get_this_from_strider( roled_stride );
				const auto  to_rep_strider   = get_this_to_strider  ( roled_stride );
				for (const auto &from_rep_index : get_rep_indices_range( from_rep_strider, num_residues ) ) {
//					BOOST_LOG_TRIVIAL( warning ) << "\tFrom rep " << from_rep_index;
					for (const auto &to_rep_index : get_rep_indices_range( to_rep_strider, num_residues ) ) {
//						BOOST_LOG_TRIVIAL( warning ) << "\t\tTo rep " << to_rep_index;
						if ( from_rep_index != to_rep_index ) {
							const auto the_res_pair = make_multi_struc_res_rep_pair(
								arg_protein.get_residue_ref_of_index( get_index_of_rep_index( from_rep_strider, from_rep_index ) ),
								arg_protein.get_residue_ref_of_index( get_index_of_rep_index( to_rep_strider, to_rep_index )   ),
								arg_structure_index,
								from_rep_index,
								to_rep_index
							);
//							BOOST_LOG_TRIVIAL( warning ) << "\t\tMade multi_struc_res_rep_pair : " << the_res_pair;
//							const auto the_key    = the_keyer.make_key       ( the_res_pair );
//							BOOST_LOG_TRIVIAL( warning ) << "\t\tThe actual key itself is      : " << output_key(the_key );

							const auto close_keys = the_keyer.make_close_keys( the_res_pair, the_criteria );
//							BOOST_LOG_TRIVIAL( warning ) << "\t\t\t(from_phi : "
//							                             << the_res_pair.get_res_pair_core().get_from_phi_angle()
//							                             << ") There are "
//							                             << boost::distance( common::cross( close_keys ) )
//							                             << " keys close to "
//							                             << output_key( the_keyer.make_key( the_res_pair ) );

							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 0 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<0>( close_keys )
							// 			| boost::adaptors::transformed( [](const uint8_t &x) { return std::to_string( static_cast<size_t>( x) ); } ),
							// 		"\", \""
							// 	)
							// 	<< "\"";
							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 1 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<1>( close_keys )
							// 			| boost::adaptors::transformed( [](const uint8_t &x) { return std::to_string( static_cast<size_t>( x) ); } ),
							// 		"\", \""
							// 	)
							// 	<< "\"";
							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 2 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<2>( close_keys )
							// 			| boost::adaptors::transformed( [](const uint8_t &x) { return std::to_string( static_cast<size_t>( x) ); } ),
							// 		"\", \""
							// 	)
							// 	<< "\"";
							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 3 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<3>( close_keys )
							// 			| boost::adaptors::transformed( [](const uint8_t &x) { return std::to_string( static_cast<size_t>( x) ); } ),
							// 		"\", \""
							// 	)
							// 	<< "\"";
							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 4 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<4>( close_keys )
							// 			| common::lexical_casted<std::string>(),
							// 		"\", \""
							// 	)
							// 	<< "\"";
							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 5 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<5>( close_keys )
							// 			| common::lexical_casted<std::string>(),
							// 		"\", \""
							// 	)
							// 	<< "\"";
							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 6 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<6>( close_keys )
							// 			| common::lexical_casted<std::string>(),
							// 		"\", \""
							// 	)
							// 	<< "\"";
							// BOOST_LOG_TRIVIAL( warning )
							// 	<< "\t\tRange 7 : \""
							// 	<< boost::algorithm::join (
							// 		std::get<7>( close_keys )
							// 			| common::lexical_casted<std::string>(),
							// 		"\", \""
							// 	)
							// 	<< "\"";

//							KEYER_PARTS: tuple<
//								res_pair_phi_psi_angle_keyer_part<res_pair_from_phi_keyer_part_spec>,
//								res_pair_phi_psi_angle_keyer_part<res_pair_from_psi_keyer_part_spec>,
//								res_pair_phi_psi_angle_keyer_part<res_pair_to_phi_keyer_part_spec>,
//								res_pair_phi_psi_angle_keyer_part<res_pair_to_psi_keyer_part_spec>,
//								res_pair_index_dirn_keyer_part,
//								res_pair_view_axis_keyer_part<res_pair_view_x_keyer_part_spec>,
//								res_pair_view_axis_keyer_part<res_pair_view_y_keyer_part_spec>,
//								res_pair_view_axis_keyer_part<res_pair_view_z_keyer_part_spec>
//							>
//							KEY: tuple<unsigned char, unsigned char, unsigned char, unsigned char, res_pair_dirn, short, short, short>
//							CLOSE_KEYS: tuple<
//								joined_range<const integer_range<unsigned char>, const integer_range<unsigned char> >,
//								joined_range<const integer_range<unsigned char>, const integer_range<unsigned char> >,
//								joined_range<const integer_range<unsigned char>, const integer_range<unsigned char> >,
//								joined_range<const integer_range<unsigned char>, const integer_range<unsigned char> >,
//								array<res_pair_dirn, 1u>,
//								integer_range<short int>,
//								integer_range<short int>,
//								integer_range<short int>
//							>;
//							ELEMENT_OF_CLOSE_KEYS: tuple<unsigned char, unsigned char, unsigned char, unsigned char, const res_pair_dirn &, short, short, short>

//							const bool contains_key = common::contains( common::cross( close_keys ), the_keyer.make_key( the_res_pair ) );
//							if ( contains_key ) {
//								BOOST_LOG_TRIVIAL( warning ) << "\t\t\tDoes contain the central key";
//							}
//							else {
//								BOOST_LOG_TRIVIAL( warning ) << "\t\t\t***** Doesn't contain the central key";
//							}

//							size_t counter = 0;
							for (const auto &x : common::cross( close_keys ) ) {
//								if ( ! contains_key ) {
//									BOOST_LOG_TRIVIAL( warning ) << "\t\t\t\t[" << counter << "] : Adding entry for close key " << output_key( x );
//								}
								arg_store.push_back_entry_to_cell( x, the_res_pair );
//								++counter;
							}
						}
					}
				}
				arg_store.summarize();
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
