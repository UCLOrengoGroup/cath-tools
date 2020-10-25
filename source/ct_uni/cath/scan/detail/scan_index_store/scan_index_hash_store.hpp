/// \file
/// \brief The scan_index_hash_store class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_HASH_STORE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_HASH_STORE_HPP

#include <boost/numeric/conversion/cast.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/type_traits.hpp"
#include "cath/scan/detail/scan_index_store/detail/hash_tuple.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"

#include <iostream> // ***** TEMPORARY *****
#include <limits>
#include <tuple>
#include <unordered_map>

namespace cath {
	namespace scan {
		namespace detail {
			namespace detail {

				/// \brief TODOCUMENT
				struct empty_key_maker final {
					/// \brief TODOCUMENT
					template <typename... Ts>
					constexpr auto operator()(const Ts &...
					                          ) const {
						return std::make_tuple( std::numeric_limits<Ts>::max()... );
					}
				};
			} // namespace detail

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			class scan_index_hash_store final {
			private:
				/// \brief TODOCUMENT
				using key_hash = hash_tuple::hash<Key>;

				/// \brief TODOCUMENT
				using value_t = common::range_value_t<Cell>;

				/// \brief TODOCUMENT
				std::unordered_map<Key, Cell, key_hash> the_store;

				/// \brief TODOCUMENT
				Cell empty_cell{};

				long long unsigned int num_adds = 0;

			public:
				using const_iterator = typename std::unordered_map<Key, Cell, key_hash>::const_iterator;

				scan_index_hash_store() {
//					const auto empty_key = common::apply( detail::empty_key_maker(), Key() );
//					the_store.set_empty_key( empty_key );

					// the_store.rehash( 131072 );
					// std::cerr << "scan_index_hash_store's ctor currently uses a hard-coded rehash to 131072 buckets!\n";
				}

				template <typename T>
				inline void push_back_entry_to_cell(const Key  &,
				                                    T  &&);

				template <typename... Ts>
				inline void emplace_back_entry_to_cell(const Key  &,
				                                       Ts &&...);

				const Cell & find_matches(const Key &) const;

				info_quantity get_info_size() const;

				const_iterator begin() const;
				const_iterator end() const;

				void summarize() const {
					// BOOST_LOG_TRIVIAL( warning )<< "scan_index_hash_store : size is            : " << the_store.size();
					// BOOST_LOG_TRIVIAL( warning )<< "scan_index_hash_store : max_size is        : " << the_store.max_size();
					// BOOST_LOG_TRIVIAL( warning )<< "scan_index_hash_store : bucket_count is    : " << the_store.bucket_count();
					// BOOST_LOG_TRIVIAL( warning )<< "scan_index_hash_store : load_factor is     : " << the_store.load_factor();
					// BOOST_LOG_TRIVIAL( warning )<< "scan_index_hash_store : max_load_factor is : " << the_store.max_load_factor();
					// BOOST_LOG_TRIVIAL( warning )<< "scan_index_hash_store : num_adds is        : " << num_adds;
				}
			};

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			template <typename T>
			inline void scan_index_hash_store<Key, Cell>::push_back_entry_to_cell(const Key  &prm_key, ///< TODOCUMENT
			                                                                      T         &&prm_data ///< TODOCUMENT
			                                                                      ) {
				the_store[ prm_key ].push_back( std::forward<T>( prm_data ) );
				++num_adds;
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			template <typename... Ts>
			inline void scan_index_hash_store<Key, Cell>::emplace_back_entry_to_cell(const Key  &    prm_key, ///< TODOCUMENT
			                                                                         Ts        &&... prm_data ///< TODOCUMENT
			                                                                         ) {
				the_store[ prm_key ].emplace_back( std::forward<Ts>( prm_data )... );
				++num_adds;
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			inline const Cell & scan_index_hash_store<Key, Cell>::find_matches(const Key &prm_key ///< TODOCUMENT
			                                                                   ) const {
				const auto &cell_itr = the_store.find( prm_key );
				return ( cell_itr == common::cend( the_store ) ) ? empty_cell : cell_itr->second;
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			info_quantity scan_index_hash_store<Key, Cell>::get_info_size() const {
				const auto num_bytes =
					  sizeof( common::remove_cvref_t< decltype( *this ) > )
					+ sizeof( Cell    ) * the_store.size()
					+ sizeof( value_t ) * num_adds;
				return num_bytes * boost::units::information::bytes;
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_hash_store<Key, Cell>::begin() const -> const_iterator {
				return common::cbegin( the_store );
			}

			/// \brief TODOCUMENT
			template <typename Key, typename Cell>
			auto scan_index_hash_store<Key, Cell>::end() const -> const_iterator {
				return common::cend  ( the_store );
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_HASH_STORE_HPP
