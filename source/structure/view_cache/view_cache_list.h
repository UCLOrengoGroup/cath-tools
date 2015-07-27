/// \file
/// \brief The view_cache_list class header

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

#ifndef VIEW_CACHE_LIST_H_INCLUDED
#define VIEW_CACHE_LIST_H_INCLUDED

#include "common/type_aliases.h"
#include "ssap/context_res.h"
#include "structure/structure_type_aliases.h"
#include "structure/view_cache/view_cache.h"

#include <vector>

namespace cath { namespace geom { class coord; } }
namespace cath { class protein_list; }

namespace cath {
	namespace index {

		/// \brief Caches of views (ie vectors implemented as coords) between pairs of residues for a list of residue lists
		///        (most likely: a list of proteins)
		///
		/// This is essentially a wrapper for vector<view_cache>
		class view_cache_list final {
		private:
			/// \brief The view_caches for the individual residues lists (probably proteins)
			view_cache_vec view_caches;

			static view_cache_vec build_caches(const protein_list &);

		public:
			view_cache_list(const protein_list &);
			
			const view_cache & get_view_cache(const size_t &) const;
		};

		/// \brief Getter to retrieve the view_cache for the specified index
		inline const view_cache & view_cache_list::get_view_cache(const size_t &arg_index ///< The index of the residue list (eg protein) of interest
		                                                          ) const {
			return view_caches[ arg_index ];
		}

		/// \brief Get the view in the specified view_cache_list for the specified
		///        residue list (eg protein) index, from_residue index and to_residue index
		///
		/// \relates view_cache_list
		inline const geom::coord & get_view(const view_cache_list &arg_cache_list, ///< The view_cache_list to query
		                                    const size_t          &arg_index,      ///< The index of the residue list (eg protein) of interest
		                                    const size_t          &arg_from_index, ///< The index of the from_residue in the residue list (eg protein)
		                                    const size_t          &arg_to_index    ///< The index of the to_residue   in the residue list (eg protein)
		                                    ) {
			return arg_cache_list.get_view_cache( arg_index ).get_view( arg_from_index, arg_to_index );
		}

		/// \brief Retrieve two from/to residue pairs from the specified view_cache_list and return the residue_context between them
		///
		/// \relates view_cache_list
		inline float_score_type get_residue_context(const view_cache_list &arg_cache_list,   ///< The view_cache_list to query
		                                            const size_t          &arg_index_a,      ///< The index of the first  residue list (eg protein) of interest
		                                            const size_t          &arg_index_b,      ///< The index of the second residue list (eg protein) of interest
		                                            const size_t          &arg_from_index_a, ///< The index of the first  from_residue in the first  residue list (eg protein)
		                                            const size_t          &arg_from_index_b, ///< The index of the second from_residue in the second residue list (eg protein)
		                                            const size_t          &arg_to_index_a,   ///< The index of the first  to_residue   in the first  residue list (eg protein)
		                                            const size_t          &arg_to_index_b    ///< The index of the second to_residue   in the second residue list (eg protein)
		                                            ) {
			/// \todo Consider replacing simplified_context_res_vec() with context_res_vec<false, distance_score_formula::SIMPLIFIED>()
			return simplified_context_res_vec(
				get_view( arg_cache_list, arg_index_a, arg_from_index_a, arg_to_index_a ),
				get_view( arg_cache_list, arg_index_b, arg_from_index_b, arg_to_index_b )
			);
		}

	}

}

#endif
