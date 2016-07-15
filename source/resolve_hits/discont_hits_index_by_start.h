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

#ifndef DISCONT_HITS_INDEX_BY_START_H_INCLUDED
#define DISCONT_HITS_INDEX_BY_START_H_INCLUDED

#include "resolve_hits/hit.h"
#include "resolve_hits/hit_list.h"

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Index the discontiguous hits in a hit_list, ordering them by their starts
			///
			/// Calculating this once at the start a resolve() enables much faster calculation
			/// of the discontiguous domains that right-intersperse a set of domains in a mask.
			class discont_hits_index_by_start final {
			private:
				/// \brief A reference to the hit_list this is indexing
				std::reference_wrapper<const hit_list> the_hits;

				/// \brief A list of the start's of the discontiguous hits and their indices in the hit_list,
				///        sorted by their starts
				res_arr_idx_pair_vec disconts;

				static res_arr_idx_pair_vec calc_disconts(const hit_list &);

			public:
				explicit discont_hits_index_by_start(const hit_list &);

				boost::integer_range<size_t> get_index_indices_of_disconts_in_range(const res_arrow &,
				                                                                    const res_arrow &) const;

				const hit & get_discont_hit_of_index_index(const size_t &) const;

			};


		}
	}
}

#endif
