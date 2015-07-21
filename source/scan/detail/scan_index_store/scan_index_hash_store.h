/// \file
/// \brief The scan_index_hash_store class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef SCAN_INDEX_HASH_STORE_H_INCLUDED
#define SCAN_INDEX_HASH_STORE_H_INCLUDED

//#include <boost/log/trivial.hpp> // ***** TEMPORARY *****
#include <boost/numeric/conversion/cast.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

//#include "common/algorithm/transform_tuple.h" // ***** TEMPORARY *****
#include "common/c++14/cbegin_cend.h"
#include "scan/detail/res_pair/multi_struc_res_rep_pair.h"
#include "scan/detail/res_pair/multi_struc_res_rep_pair_list.h"
#include "scan/detail/scan_index_store/detail/hash_tuple.h"
#include "scan/detail/scan_type_aliases.h"

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
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			class scan_index_hash_store final {
			private:
				/// \brief TODOCUMENT
				using key_hash = hash_tuple::hash<KEY>;

				/// \brief TODOCUMENT
				std::unordered_map<KEY, multi_struc_res_rep_pair_list, key_hash> the_store;

				/// \brief TODOCUMENT
				const multi_struc_res_rep_pair_list empty_cell{};

				long long unsigned int num_adds = 0;

			public:
				scan_index_hash_store() {
//					const auto empty_key = common::transform_tuple( KEY(), detail::empty_key_maker() );
//					the_store.set_empty_key( empty_key );
					the_store.rehash( 131072 );
				}
				void add_entry(const KEY &,
				               const multi_struc_res_rep_pair &);
				const multi_struc_res_rep_pair_list & find_matches(const KEY &) const;

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
			template <typename KEY>
			inline void scan_index_hash_store<KEY>::add_entry(const KEY                      &arg_key,     ///< TODOCUMENT
			                                                  const multi_struc_res_rep_pair &arg_res_pair ///< TODOCUMENT
			                                                  ) {
				the_store[ arg_key ].emplace_back( arg_res_pair );
				++num_adds;
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			inline const multi_struc_res_rep_pair_list & scan_index_hash_store<KEY>::find_matches(const KEY &arg_key ///< TODOCUMENT
			                                                                                      ) const {
				const auto &cell_itr = the_store.find( arg_key );
				return ( cell_itr == common::cend( the_store ) ) ? empty_cell : cell_itr->second;
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			info_quantity scan_index_hash_store<KEY>::get_info_size() const {

				const info_value num_bytes = boost::numeric_cast<info_value>( num_adds * sizeof( multi_struc_res_rep_pair ) );
				return num_bytes * boost::units::information::bytes;
			}

		}
	}
}

#endif
