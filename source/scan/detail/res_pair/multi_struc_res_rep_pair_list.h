/// \file
/// \brief The multi_struc_res_rep_pair_list class header

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

#ifndef MULTI_STRUC_RES_REP_PAIR_LIST_H_INCLUDED
#define MULTI_STRUC_RES_REP_PAIR_LIST_H_INCLUDED

#include "common/c++14/cbegin_cend.h"
#include "scan/detail/scan_type_aliases.h"

#include <cstddef>
#include <type_traits>

namespace cath { namespace scan { namespace detail { class multi_struc_res_rep_pair; } } }

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Store an ordered list of multi_struc_res_rep_pair objects
			///
			/// This is useful for implementing a cell of reps for all-vs-all scanning
			class multi_struc_res_rep_pair_list final {
			private:
				/// \brief The multi_struc_res_rep_pairs, stored in a vector
				multi_struc_res_rep_pair_vec multi_struc_res_rep_pairs;

			public:
				/// \brief const_iterator type alias at part of making this a range over the multi_struc_res_rep_pairs
				using const_iterator = multi_struc_res_rep_pair_vec_citr;

				/// \brief Default ctor
				multi_struc_res_rep_pair_list() = default;
				multi_struc_res_rep_pair_list(const multi_struc_res_rep_pair_vec &);

				bool empty() const;
				size_t size() const;
				const multi_struc_res_rep_pair & operator[](const size_t &) const;

				template <class... Ts>
				void emplace_back(Ts&& ...);
				void push_back(const multi_struc_res_rep_pair &);

				const_iterator begin() const;
				const_iterator end() const;
			};

			/// \brief Ctor from a vector of multi_struc_res_rep_pair objects
			inline multi_struc_res_rep_pair_list::multi_struc_res_rep_pair_list(const multi_struc_res_rep_pair_vec &arg_multi_struc_res_rep_pairs ///< The vector of multi_struc_res_rep_pairs from which to construct the multi_struc_res_rep_pair_list
			                                                                    ) : multi_struc_res_rep_pairs ( arg_multi_struc_res_rep_pairs ) {
			}

			/// \brief Return whether this multi_struc_res_rep_pair_list is empty
			inline bool multi_struc_res_rep_pair_list::empty() const {
				return multi_struc_res_rep_pairs.empty();
			}

			/// \brief Return the number of multi_struc_res_rep_pair entries
			inline size_t multi_struc_res_rep_pair_list::size() const {
				return multi_struc_res_rep_pairs.size();
			}

			/// \brief Standard subscript operator
			inline const multi_struc_res_rep_pair & multi_struc_res_rep_pair_list::operator[](const size_t &arg_index ///< The index of the entry to be queried
			                                                                                  ) const {
				return multi_struc_res_rep_pairs[ arg_index ];
			}

			/// \brief TODOCUMENT
			template <class... Ts>
			void multi_struc_res_rep_pair_list::emplace_back(Ts&& ... arg_values ///< TODOCUMENT
			                                                 ) {
				multi_struc_res_rep_pairs.emplace_back( std::forward<Ts>( arg_values ) ... );
			}

			/// \brief TODOCUMENT
			inline void multi_struc_res_rep_pair_list::push_back(const multi_struc_res_rep_pair &arg_res_pair ///< TODOCUMENT
			                                                     ) {
				multi_struc_res_rep_pairs.push_back( arg_res_pair );
			}

			/// \brief Standard const begin() method
			inline auto multi_struc_res_rep_pair_list::begin() const -> const_iterator {
				return common::cbegin( multi_struc_res_rep_pairs );
			}

			/// \brief Standard const end() method
			inline auto multi_struc_res_rep_pair_list::end() const -> const_iterator {
				return common::cend  ( multi_struc_res_rep_pairs );
			}

		}
	}
}

#endif
