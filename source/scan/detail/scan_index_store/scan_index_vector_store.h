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

#ifndef SCAN_INDEX_VECTOR_STORE_H_INCLUDED
#define SCAN_INDEX_VECTOR_STORE_H_INCLUDED

//#include <boost/log/expressions.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/range/numeric.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "common/c++14/cbegin_cend.h"
#include "common/size_t_literal.h"
//#include "exception/not_implemented_exception.h"
#include "scan/detail/res_pair/multi_struc_res_rep_pair.h"
#include "scan/detail/res_pair/multi_struc_res_rep_pair_list.h"
#include "scan/detail/scan_type_aliases.h"

using namespace cath::common::literals;

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename KEY>
			class scan_index_vector_store final {
			private:
				/// \brief TODOCUMENT
				key_multi_res_pair_list_pair_vec<KEY> the_store;

				/// \brief TODOCUMENT
				multi_struc_res_rep_pair_list empty_cell;

				multi_struc_res_rep_pair_list & find_or_create_cell(const KEY &);

			public:
				/// \brief TODOCUMENT
				using const_iterator = key_multi_res_pair_list_pair_vec_citr<KEY>;

				void add_entry(const KEY &,
				               const multi_struc_res_rep_pair &);
				const multi_struc_res_rep_pair_list & find_matches(const KEY &) const;

				const_iterator begin() const;
				const_iterator end() const;

				info_quantity get_info_size() const;
			};

			/// \brief TODOCUMENT
			template <typename KEY>
			inline multi_struc_res_rep_pair_list & scan_index_vector_store<KEY>::find_or_create_cell(const KEY &arg_key ///< TODOCUMENT
			                                                                                         ) {
				const auto cell_itr = boost::range::lower_bound(
					the_store,
					arg_key,
					[] (const key_multi_res_pair_list_pair<KEY> &x, const KEY &y) { return x.first < y; }
				);
				if ( cell_itr != common::cend( the_store ) && cell_itr->first == arg_key ) {
					return cell_itr->second;
				}
				else {
					return the_store.insert( cell_itr, make_pair( arg_key, multi_struc_res_rep_pair_list{} ) )->second;
				}
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			inline void scan_index_vector_store<KEY>::add_entry(const KEY                      &arg_key,     ///< TODOCUMENT
			                                                    const multi_struc_res_rep_pair &arg_res_pair ///< TODOCUMENT
			                                                    ) {
				find_or_create_cell( arg_key ).push_back( arg_res_pair );
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			inline const multi_struc_res_rep_pair_list & scan_index_vector_store<KEY>::find_matches(const KEY &arg_key ///< TODOCUMENT
			                                                                                        ) const {
				const auto cell_itr = boost::range::lower_bound(
					the_store,
					arg_key,
					[] (const key_multi_res_pair_list_pair<KEY> &x, const KEY &y) { return x.first < y; }
				);
				return ( cell_itr == common::cend( the_store ) || cell_itr->first != arg_key ) ? empty_cell
				                                                                               : cell_itr->second;
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			auto scan_index_vector_store<KEY>::begin() const -> const_iterator {
				return common::cbegin( the_store );
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			auto scan_index_vector_store<KEY>::end() const -> const_iterator {
				return common::cend( the_store );
			}

			/// \brief TODOCUMENT
			template <typename KEY>
			info_quantity scan_index_vector_store<KEY>::get_info_size() const {
				const auto num_bytes = sizeof( multi_struc_res_rep_pair ) * boost::accumulate(
					the_store
						| boost::adaptors::map_values
						| boost::adaptors::transformed( [] (const multi_struc_res_rep_pair_list &x) { return x.size(); } ),
					0_z
				);
				return num_bytes * boost::units::information::bytes;
			}

		}
	}
}

#endif
