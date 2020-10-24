/// \file
/// \brief The discont_hits_index_by_start class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_DISCONT_HITS_INDEX_BY_START_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_DISCONT_HITS_INDEX_BY_START_HPP

#include "cath/resolve_hits/calc_hit.hpp"
#include "cath/resolve_hits/calc_hit_list.hpp"

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Index the discontiguous hits in a calc_hit_list, ordering them by their starts
			///
			/// Calculating this once at the start a resolve() enables much faster calculation
			/// of the discontiguous domains that right-intersperse a set of domains in a mask.
			class discont_hits_index_by_start final {
			private:
				/// \brief A reference to the calc_hit_list this is indexing
				std::reference_wrapper<const calc_hit_list> the_hits;

				/// \brief A list of the start's of the discontiguous hits and their indices in the calc_hit_list,
				///        sorted by their starts
				res_arr_idx_pair_vec disconts;

				static res_arr_idx_pair_vec calc_disconts(const calc_hit_list &);

			public:
				explicit discont_hits_index_by_start(const calc_hit_list &);

				boost::integer_range<size_t> get_index_indices_of_disconts_in_range(const seq::seq_arrow &,
				                                                                    const seq::seq_arrow &) const;

				size_t size() const;
				const calc_hit & get_discont_hit_of_index_index(const size_t &) const;
			};

			std::string to_string(const discont_hits_index_by_start &);
		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_DISCONT_HITS_INDEX_BY_START_HPP
