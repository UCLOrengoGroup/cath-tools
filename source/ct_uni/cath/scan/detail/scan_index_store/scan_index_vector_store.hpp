/// \file
/// \brief The scan_index_vector_store class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_VECTOR_STORE_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_VECTOR_STORE_HPP

#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/numeric.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"

using namespace ::cath::common::literals;

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			class scan_index_vector_store final {
			private:
				/// \brief TODOCUMENT
				using key_cell_pair = std::pair<Key, Cell>;

				/// \brief TODOCUMENT
				using key_cell_pair_vec = std::vector<key_cell_pair>;

				/// \brief TODOCUMENT
				using value_t = common::range_value_t<Cell>;

				/// \brief TODOCUMENT
				key_cell_pair_vec the_store;

				/// \brief TODOCUMENT
				Cell empty_cell;

				/// \brief TODOCUMENT
				Cell & find_or_create_cell(const Key &);

			public:
				/// \brief TODOCUMENT
				using const_iterator = typename key_cell_pair_vec::const_iterator;

				void push_back_entry_to_cell(const Key &,
				                             const value_t &);

				template <typename... Ts>
				void emplace_back_entry_to_cell(const Key &,
				                                Ts &&...);

				const Cell & find_matches(const Key &) const;

				const_iterator begin() const;
				const_iterator end() const;

				info_quantity get_info_size() const;
			};

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_vector_store<Key, Cell>::find_or_create_cell(const Key &prm_key ///< TODOCUMENT
			                                                                    ) -> Cell & {
				const auto cell_itr = boost::range::lower_bound(
					the_store,
					prm_key,
					[] (const std::pair<Key, Cell> &x, const Key &y) {/* TD< decltype( x ) > fred; */return x.first < y; }
				);
				if ( cell_itr != common::cend( the_store ) && cell_itr->first == prm_key ) {
					return cell_itr->second;
				}
				return the_store.insert( cell_itr, make_pair( prm_key, Cell{} ) )->second;
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline void scan_index_vector_store<Key, Cell>::push_back_entry_to_cell(const Key     &prm_key, ///< TODOCUMENT
			                                                                        const value_t &prm_data ///< TODOCUMENT
			                                                                        ) {
				find_or_create_cell( prm_key ).push_back( prm_data );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			template <typename... Ts>
			inline void scan_index_vector_store<Key, Cell>::emplace_back_entry_to_cell(const Key  &    prm_key, ///< TODOCUMENT
			                                                                           Ts        &&... prm_data ///< TODOCUMENT
			                                                                           ) {
				find_or_create_cell( prm_key ).emplace_back( std::forward<Ts>( prm_data )... );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline auto scan_index_vector_store<Key, Cell>::find_matches(const Key &prm_key ///< TODOCUMENT
			                                                             ) const -> const Cell & {
				const auto cell_itr = boost::range::lower_bound(
					the_store,
					prm_key,
					[] (const std::pair<Key, Cell> &x, const Key &y) { return x.first < y; }
				);
				return ( cell_itr == common::cend( the_store ) || cell_itr->first != prm_key ) ? empty_cell
				                                                                               : cell_itr->second;
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_vector_store<Key, Cell>::begin() const -> const_iterator {
				return common::cbegin( the_store );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_vector_store<Key, Cell>::end() const -> const_iterator {
				return common::cend( the_store );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			info_quantity scan_index_vector_store<Key, Cell>::get_info_size() const {
				const auto num_bytes = sizeof( value_t ) * boost::accumulate(
					the_store
						| boost::adaptors::map_values
						| boost::adaptors::transformed( [] (const Cell &x) { return x.size(); } ),
					0_z
				);
				return num_bytes * boost::units::information::bytes;
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
