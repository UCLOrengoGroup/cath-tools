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

#ifndef _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_HASH_STORE_H
#define _CATH_TOOLS_SOURCE_SCAN_DETAIL_SCAN_INDEX_STORE_SCAN_INDEX_HASH_STORE_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "scan/detail/res_pair/multi_struc_res_rep_pair_list.hpp"
#include "scan/detail/scan_index_store/detail/hash_tuple.hpp"
#include "scan/detail/scan_type_aliases.hpp"

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
			template <typename Key>
			class scan_index_hash_store final {
			private:
				/// \brief TODOCUMENT
				using key_hash = hash_tuple::hash<Key>;

				/// \brief TODOCUMENT
				std::unordered_map<Key, multi_struc_res_rep_pair_list, key_hash> the_store;

				/// \brief TODOCUMENT
				const multi_struc_res_rep_pair_list empty_cell{};

				long long unsigned int num_adds = 0;

			public:
				scan_index_hash_store() {
//					const auto empty_key = common::apply( detail::empty_key_maker(), Key() );
//					the_store.set_empty_key( empty_key );
					the_store.rehash( 131072 );
					std::cerr << "scan_index_hash_store's ctor currently uses a hard-coded rehash to 131072 buckets!\n";
				}

				template <typename T>
				inline void push_back_entry_to_cell(const Key  &,
				                                    T  &&);

				template <typename... Ts>
				inline void emplace_back_entry_to_cell(const Key  &,
				                                       Ts &&...);

				const multi_struc_res_rep_pair_list & find_matches(const Key &) const;

				info_quantity get_info_size() const;

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
			template <typename Key>
			template <typename T>
			inline void scan_index_hash_store<Key>::push_back_entry_to_cell(const Key  &arg_key, ///< TODOCUMENT
			                                                                T         &&arg_data ///< TODOCUMENT
			                                                                ) {
				the_store[ arg_key ].push_back( std::forward<T>( arg_data ) );
				++num_adds;
			}

			/// \brief TODOCUMENT
			template <typename Key>
			template <typename... Ts>
			inline void scan_index_hash_store<Key>::emplace_back_entry_to_cell(const Key  &    arg_key, ///< TODOCUMENT
			                                                                   Ts        &&... arg_data ///< TODOCUMENT
			                                                                   ) {
				the_store[ arg_key ].emplace_back( std::forward<Ts>( arg_data )... );
				++num_adds;
			}

			/// \brief TODOCUMENT
			template <typename Key>
			inline const multi_struc_res_rep_pair_list & scan_index_hash_store<Key>::find_matches(const Key &arg_key ///< TODOCUMENT
			                                                                                      ) const {
				const auto &cell_itr = the_store.find( arg_key );
				return ( cell_itr == common::cend( the_store ) ) ? empty_cell : cell_itr->second;
			}

			/// \brief TODOCUMENT
			template <typename Key>
			info_quantity scan_index_hash_store<Key>::get_info_size() const {

				const info_value num_bytes = boost::numeric_cast<info_value>( num_adds * sizeof( multi_struc_res_rep_pair ) );
				return num_bytes * boost::units::information::bytes;
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
