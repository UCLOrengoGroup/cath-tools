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

#ifndef _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_LATTICE_STORE_H
#define _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_LATTICE_STORE_H

#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/join.hpp>
#include <boost/range/numeric.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/cpp17/apply.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/size_t_literal.hpp"
#include "common/tuple/mins_maxs_tuple_pair_mins_maxs_element.hpp"
#include "common/tuple/tuple_increment.hpp"
#include "common/tuple/tuple_lattice_index.hpp"
#include "common/tuple/tuple_mins_maxs_element.hpp"
#include "common/tuple/tuple_multiply_args.hpp"
#include "common/tuple/tuple_subtract.hpp"
#include "common/tuple/tuple_within_range.hpp"
#include "scan/detail/scan_type_aliases.hpp"

#include <utility>

using namespace cath::common::literals;

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

				info_quantity get_info_size() const;
			};

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_lattice_store<Key, Cell>::find_cell(const Key &arg_key ///< TODOCUMENT
			                                                           ) -> Cell & {
				return const_cast<Cell &>( static_cast<const scan_index_lattice_store &>( *this ).find_cell( arg_key ) );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_lattice_store<Key, Cell>::find_cell(const Key &arg_key ///< TODOCUMENT
			                                                           ) const -> const Cell & {
				return the_store[
					debug_numeric_cast<size_t>(
						common::tuple_lattice_index(
							common::tuple_subtract( arg_key, mins_key ),
							nums_of_cells_key
						)
					)
				].second;
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline scan_index_lattice_store<Key, Cell>::scan_index_lattice_store(const Key &arg_mins_key, ///< TODOCUMENT
			                                                                     const Key &arg_maxs_key  ///< TODOCUMENT
			                                                                     ) : mins_key         { arg_mins_key                                        },
			                                                                         nums_of_cells_key{
			                                                                         	common::tuple_increment(
			                                                                         		common::tuple_subtract( arg_maxs_key, arg_mins_key )
			                                                                         	)
			                                                                         },
			                                                                         the_store        {
			                                                                         	boost::numeric_cast<size_t>( common::tuple_multiply_args( nums_of_cells_key ) )
			                                                                         } {
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline void scan_index_lattice_store<Key, Cell>::push_back_entry_to_cell(const Key     &arg_key, ///< TODOCUMENT
			                                                                         const value_t &arg_data ///< TODOCUMENT
			                                                                         ) {
				find_cell( arg_key ).push_back( arg_data );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			template <typename... Ts>
			inline void scan_index_lattice_store<Key, Cell>::emplace_back_entry_to_cell(const Key  &    arg_key, ///< TODOCUMENT
			                                                                            Ts        &&... arg_data ///< TODOCUMENT
			                                                                            ) {
				find_cell( arg_key ).emplace_back( std::forward<Ts>( arg_data )... );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline bool scan_index_lattice_store<Key, Cell>::has_matches(const Key &arg_key ///< TODOCUMENT
			                                                             ) const {
				return common::tuple_within_range( common::tuple_subtract( arg_key, mins_key ), nums_of_cells_key );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_lattice_store<Key, Cell>::find_matches(const Key &arg_key ///< TODOCUMENT
			                                                              ) const -> const Cell & {
				return find_cell( arg_key );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_lattice_store<Key, Cell>::begin() const -> const_iterator {
				return common::cbegin( the_store );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_lattice_store<Key, Cell>::end() const -> const_iterator {
				return common::cend( the_store );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			info_quantity scan_index_lattice_store<Key, Cell>::get_info_size() const {
				const auto num_bytes =
					  sizeof( std::decay_t< decltype( *this ) > )
					+ sizeof( Cell ) * the_store.size()
					+ sizeof( value_t ) * boost::accumulate(
						the_store
							| boost::adaptors::map_values
							| boost::adaptors::transformed( [] (const Cell &x) { return x.size(); } ),
						0_z
					);
				return num_bytes * boost::units::information::bytes;
			}

			/// \brief TODOCUMENT
			template <typename Cell, typename Key>
			scan_index_lattice_store<Key, Cell> make_scan_index_lattice_store(const Key &arg_mins_key, ///< TODOCUMENT
			                                                                  const Key &arg_maxs_key  ///< TODOCUMENT
			                                                                  ) {
				return { arg_mins_key, arg_maxs_key };
			}

			/// \brief TODOCUMENT
			template <typename Rng, typename Keyer>
			auto make_sparse_empty_lattice_store(const Rng   &arg_rng,  ///< TODOCUMENT
			                                     const Keyer &arg_keyer ///< TODOCUMENT
			                                     ) {
				using value_t = common::range_value_t< Rng >;

				const auto mins_maxs = common::tuple_mins_maxs_element(
					arg_rng
						| boost::adaptors::transformed( [&] (const value_t &value) { return arg_keyer.make_key( value ); } )
				);

				return make_scan_index_lattice_store< std::vector<value_t> >( mins_maxs.first, mins_maxs.second );
			}

			/// \brief TODOCUMENT
			template <typename Rng, typename Keyer, typename Crit>
			auto make_dense_empty_lattice_store(const Rng   &arg_rng,   ///< TODOCUMENT
			                                    const Keyer &arg_keyer, ///< TODOCUMENT
			                                    const Crit  &arg_crit   ///< TODOCUMENT
			                                    ) {
				using value_t = common::range_value_t< Rng >;

				const auto mins_maxs = common::mins_maxs_tuple_pair_mins_maxs_element(
					arg_rng
						| boost::adaptors::transformed( [&] (const value_t &value) {
							return std::make_pair(
								arg_keyer.make_min_close_key( value, arg_crit ),
								arg_keyer.make_max_close_key( value, arg_crit )
							);
						} )
				);

				return make_scan_index_lattice_store< std::vector<value_t> >( mins_maxs.first, mins_maxs.second );
			}

			/// \brief TODOCUMENT
			template <typename Rng, typename Keyer>
			auto make_sparse_lattice_store(const Rng   &arg_rng,  ///< TODOCUMENT
			                               const Keyer &arg_keyer ///< TODOCUMENT
			                               ) {
				auto the_store = make_sparse_empty_lattice_store( arg_rng, arg_keyer );
				for (const auto &value : arg_rng) {
					the_store.push_back_entry_to_cell( arg_keyer.make_key( value ), value );
				}
				return the_store;
			}

			/// \brief TODOCUMENT
			template <typename Rng, typename Keyer, typename Crit>
			auto make_dense_lattice_store(const Rng   &arg_rng,   ///< TODOCUMENT
			                              const Keyer &arg_keyer, ///< TODOCUMENT
			                              const Crit  &arg_crit   ///< TODOCUMENT
			                              ) {
				auto the_store = make_dense_empty_lattice_store( arg_rng, arg_keyer, arg_crit );
				for (const auto &value : arg_rng) {
					const auto close_keys = arg_keyer.make_close_keys( value, arg_crit );
					for (const auto &the_key : common::cross( close_keys ) ) {
						the_store.push_back_entry_to_cell( the_key, value );
					}
				}
				return the_store;
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
