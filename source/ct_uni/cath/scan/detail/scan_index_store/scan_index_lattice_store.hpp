/// \file
/// \brief The scan_index_lattice_store class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_LATTICE_STORE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_LATTICE_STORE_HPP

#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/join.hpp>
#include <boost/range/numeric.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "cath/common/cpp17/apply.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/tuple/mins_maxs_tuple_pair_mins_maxs_element.hpp"
#include "cath/common/tuple/tuple_increment.hpp"
#include "cath/common/tuple/tuple_lattice_index.hpp"
#include "cath/common/tuple/tuple_mins_maxs_element.hpp"
#include "cath/common/tuple/tuple_multiply_args.hpp"
#include "cath/common/tuple/tuple_subtract.hpp"
#include "cath/common/tuple/tuple_within_range.hpp"
#include "cath/common/type_traits.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"

#include <utility>

using namespace ::cath::common::literals;

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			class scan_index_lattice_store final {
			private:
				/// \brief TODOCUMENT
				using key_cell_pair = std::pair<Key, Cell>;

				/// \brief TODOCUMENT
				using key_cell_pair_vec = std::vector<key_cell_pair>;

				/// \brief TODOCUMENT
				using value_t = common::range_value_t<Cell>;

				/// \brief TODOCUMENT
				Key mins_key;

				/// \brief TODOCUMENT
				Key nums_of_cells_key;

				/// \brief TODOCUMENT
				key_cell_pair_vec the_store;

				// /// \brief TODOCUMENT
				// Cell empty_cell;

				/// \brief const-agnostic implementation of find_cell()
				///
				/// See GSL rule: Pro.Type.3: Don't use const_cast to cast away const (i.e., at all)
				/// (https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Pro-type-constcast)
				template <typename Store>
				static auto find_cell_impl(Store     &prm_store, ///< TODOCUMENT
				                           const Key &prm_key    ///< TODOCUMENT
				                           ) -> decltype( prm_store.find_cell( prm_key ) ) {
					return prm_store.the_store[
						debug_numeric_cast<size_t>(
							common::tuple_lattice_index(
								common::tuple_subtract( prm_key, prm_store.mins_key ),
								prm_store.nums_of_cells_key
							)
						)
					].second;
				}

				Cell & find_cell(const Key &);
				const Cell & find_cell(const Key &) const;

			public:
				/// \brief TODOCUMENT
				using const_iterator = typename key_cell_pair_vec::const_iterator;

				scan_index_lattice_store(const Key &,
				                         const Key &);

				void push_back_entry_to_cell(const Key &,
				                             const value_t &);

				template <typename... Ts>
				void emplace_back_entry_to_cell(const Key &,
				                                Ts &&...);

				bool has_matches(const Key &) const;

				const Cell & find_matches(const Key &) const;

				const_iterator begin() const;
				const_iterator end() const;

				[[nodiscard]] info_quantity get_info_size() const;
			};

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_lattice_store<Key, Cell>::find_cell(const Key &prm_key ///< TODOCUMENT
			                                                           ) -> Cell & {
				return find_cell_impl( *this, prm_key );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_lattice_store<Key, Cell>::find_cell(const Key &prm_key ///< TODOCUMENT
			                                                           ) const -> const Cell & {
				return find_cell_impl( *this, prm_key );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline scan_index_lattice_store<Key, Cell>::scan_index_lattice_store(const Key &prm_mins_key, ///< TODOCUMENT
			                                                                     const Key &prm_maxs_key  ///< TODOCUMENT
			                                                                     ) : mins_key         { prm_mins_key                                        },
			                                                                         nums_of_cells_key{
			                                                                         	common::tuple_increment(
			                                                                         		common::tuple_subtract( prm_maxs_key, prm_mins_key )
			                                                                         	)
			                                                                         },
			                                                                         the_store        {
			                                                                         	boost::numeric_cast<size_t>( common::tuple_multiply_args( nums_of_cells_key ) )
			                                                                         } {
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline void scan_index_lattice_store<Key, Cell>::push_back_entry_to_cell(const Key     &prm_key, ///< TODOCUMENT
			                                                                         const value_t &prm_data ///< TODOCUMENT
			                                                                         ) {
				find_cell( prm_key ).push_back( prm_data );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			template <typename... Ts>
			inline void scan_index_lattice_store<Key, Cell>::emplace_back_entry_to_cell(const Key  &    prm_key, ///< TODOCUMENT
			                                                                            Ts        &&... prm_data ///< TODOCUMENT
			                                                                            ) {
				find_cell( prm_key ).emplace_back( std::forward<Ts>( prm_data )... );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline bool scan_index_lattice_store<Key, Cell>::has_matches(const Key &prm_key ///< TODOCUMENT
			                                                             ) const {
				return common::tuple_within_range( common::tuple_subtract( prm_key, mins_key ), nums_of_cells_key );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_lattice_store<Key, Cell>::find_matches(const Key &prm_key ///< TODOCUMENT
			                                                              ) const -> const Cell & {
				return find_cell( prm_key );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_lattice_store<Key, Cell>::begin() const -> const_iterator {
				return ::std::cbegin( the_store );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_lattice_store<Key, Cell>::end() const -> const_iterator {
				return ::std::cend( the_store );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			info_quantity scan_index_lattice_store<Key, Cell>::get_info_size() const {
				const auto num_bytes =
					  sizeof( common::remove_cvref_t< decltype( *this ) > )
					+ sizeof( Cell ) * the_store.size()
					+ sizeof( value_t ) * boost::accumulate(
						the_store
							| boost::adaptors::map_values
							| boost::adaptors::transformed( [] (const Cell &x) { return x.size(); } ),
						0_z
					);
				return num_bytes * boost::units::information::bytes;
			}

			// /// \brief TODOCUMENT
			// template <typename Cell, typename Key>
			// scan_index_lattice_store<Key, Cell> make_scan_index_lattice_store(const Key &prm_mins_key, ///< TODOCUMENT
			//                                                                   const Key &prm_maxs_key  ///< TODOCUMENT
			//                                                                   ) {
			// 	return { prm_mins_key, prm_maxs_key };
			// }

			// /// \brief TODOCUMENT
			// template <typename Rng, typename Keyer>
			// auto make_sparse_empty_lattice_store(const Rng   &prm_rng,  ///< TODOCUMENT
			//                                      const Keyer &prm_keyer ///< TODOCUMENT
			//                                      ) {
			// 	using value_t = common::range_value_t< Rng >;

			// 	const auto mins_maxs = common::tuple_mins_maxs_element(
			// 		prm_rng
			// 			| boost::adaptors::transformed( [&] (const value_t &value) { return prm_keyer.make_key( value ); } )
			// 	);

			// 	return make_scan_index_lattice_store< std::vector<value_t> >( mins_maxs.first, mins_maxs.second );
			// }

			// /// \brief TODOCUMENT
			// template <typename Rng, typename Keyer, typename Crit>
			// auto make_dense_empty_lattice_store(const Rng   &prm_rng,   ///< TODOCUMENT
			//                                     const Keyer &prm_keyer, ///< TODOCUMENT
			//                                     const Crit  &prm_crit   ///< TODOCUMENT
			//                                     ) {
			// 	using value_t = common::range_value_t< Rng >;

			// 	const auto mins_maxs = common::mins_maxs_tuple_pair_mins_maxs_element(
			// 		prm_rng
			// 			| boost::adaptors::transformed( [&] (const value_t &value) {
			// 				return std::make_pair(
			// 					prm_keyer.make_min_close_key( value, prm_crit ),
			// 					prm_keyer.make_max_close_key( value, prm_crit )
			// 				);
			// 			} )
			// 	);

			// 	return make_scan_index_lattice_store< std::vector<value_t> >( mins_maxs.first, mins_maxs.second );
			// }

			// /// \brief TODOCUMENT
			// template <typename Rng, typename Keyer>
			// auto make_sparse_lattice_store(const Rng   &prm_rng,  ///< TODOCUMENT
			//                                const Keyer &prm_keyer ///< TODOCUMENT
			//                                ) {
			// 	auto the_store = make_sparse_empty_lattice_store( prm_rng, prm_keyer );
			// 	for (const auto &value : prm_rng) {
			// 		the_store.push_back_entry_to_cell( prm_keyer.make_key( value ), value );
			// 	}
			// 	return the_store;
			// }

			// /// \brief TODOCUMENT
			// template <typename Rng, typename Keyer, typename Crit>
			// auto make_dense_lattice_store(const Rng   &prm_rng,   ///< TODOCUMENT
			//                               const Keyer &prm_keyer, ///< TODOCUMENT
			//                               const Crit  &prm_crit   ///< TODOCUMENT
			//                               ) {
			// 	auto the_store = make_dense_empty_lattice_store( prm_rng, prm_keyer, prm_crit );
			// 	for (const auto &value : prm_rng) {
			// 		const auto close_keys = prm_keyer.make_close_keys( value, prm_crit );
			// 		for (const auto &the_key : common::cross( close_keys ) ) {
			// 			the_store.push_back_entry_to_cell( the_key, value );
			// 		}
			// 	}
			// 	return the_store;
			// }

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_LATTICE_STORE_HPP
