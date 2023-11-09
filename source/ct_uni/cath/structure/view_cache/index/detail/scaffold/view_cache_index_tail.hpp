/// \file
/// \brief The view_cache_index_tail class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_SCAFFOLD_VIEW_CACHE_INDEX_TAIL_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_SCAFFOLD_VIEW_CACHE_INDEX_TAIL_HPP

#include <boost/range/algorithm_ext/for_each.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>
#include <boost/tuple/tuple.hpp>

#include "cath/structure/view_cache/index/detail/vcie_match_criteria.hpp"
#include "cath/structure/view_cache/index/view_cache_index_entry.hpp"

#include <iostream> // ***** TEMPORARY *****
#include <vector>

namespace cath::index::detail {

	/// \brief TODOCUMENT
	class view_cache_index_tail final {
	private:
		/// \brief TODOCUMENT
		view_cache_index_entry_vec entries;

	public:
		void store(const view_cache_index_entry &);

		void store(const view_cache_index_entry &,
		           const boost::tuples::null_type &);

		[[nodiscard]] bool   empty() const;
		[[nodiscard]] size_t get_num_cells() const;
		[[nodiscard]] size_t dims_remaining() const;

		template <typename ACTN>
		void perform_action_on_all_match_at_leaves(ACTN &) const;

		template <typename ACTN>
		void perform_action_on_matches(const view_cache_index_entry &,
		                               const vcie_match_criteria &,
		                               ACTN &) const;

		/// \brief TODOCUMENT
		template <typename ACTN>
		void perform_action_on_all_match_at_nodes(const view_cache_index_tail &,
		                                          const vcie_match_criteria &,
		                                          ACTN &) const;

		/// \brief TODOCUMENT
		constexpr static size_t num_dims_remaining = 0;
	};

	/// \brief TODOCUMENT
	inline void view_cache_index_tail::store(const view_cache_index_entry &prm_entry ///< TODOCUMENT
	                                         ) {
		entries.push_back( prm_entry );
	}

	/// \brief TODOCUMENT
	inline void view_cache_index_tail::store(const view_cache_index_entry   &prm_entry,  ///< TODOCUMENT
	                                         const boost::tuples::null_type &/*prm_dfs*/ ///< TODOCUMENT
	                                         ) {
		store( prm_entry );
	}

	/// \brief TODOCUMENT
	inline bool view_cache_index_tail::empty() const {
		return entries.empty();
	}

	/// \brief TODOCUMENT
	inline size_t view_cache_index_tail::get_num_cells() const {
		return entries.size();
	}

	/// \brief TODOCUMENT
	template <typename ACTN>
	void view_cache_index_tail::perform_action_on_all_match_at_leaves(ACTN &prm_action ///< TODOCUMENT
	                                                                  ) const {
		for (const view_cache_index_entry &entry : entries) {
			prm_action( entry );
		}
	}

	template <typename ACTN>
	void view_cache_index_tail::perform_action_on_matches(const view_cache_index_entry &prm_entry,    ///< TODOCUMENT
	                                                      const vcie_match_criteria    &prm_criteria, ///< TODOCUMENT
	                                                      ACTN                         &prm_action    ///< TODOCUMENT
	                                                      ) const {
		for (const view_cache_index_entry &entry : entries) {
			if ( prm_criteria( prm_entry, entry ) ) {
				// const double sq_dist = squared_distance( prm_entry, entry );
				// const double score   = ( 10.0 / ( sq_dist + 10.0 ) );
					// std::cerr << "index_match ";
				// std::cerr << std::right << std::setw( 4 ) << prm_entry.get_from_index();
				// std::cerr << "";
				// std::cerr << std::right << std::setw( 4 ) << prm_entry.get_to_index();
				// std::cerr << " ";
				// std::cerr << std::right << std::setw( 4 ) << entry.get_from_index();
				// std::cerr << " ";
				// std::cerr << std::right << std::setw( 4 ) << entry.get_to_index();
				// std::cerr << " " << score;
				// std::cerr << std::endl;
				prm_action( prm_entry, entry );
			}
		}
	}

	/// \brief TODOCUMENT
	template <typename ACTN>
	void view_cache_index_tail::perform_action_on_all_match_at_nodes(const view_cache_index_tail &prm_match_tail, ///< TODOCUMENT
	                                                                 const vcie_match_criteria   &prm_criteria,   ///< TODOCUMENT
	                                                                 ACTN                        &prm_action      ///< TODOCUMENT
	                                                                 ) const {
		// std::cerr << "Comparing "           << entries.size()
		//           << " query entries with " << prm_match_tail.entries.size()
		//           << " match entries"       << "\n";

		for (const view_cache_index_entry &x : entries) {
			for (const view_cache_index_entry &y : prm_match_tail.entries) {
				if ( prm_criteria( x, y ) ) {
					prm_action( x, y );
				}
			}
		}
	}

} // namespace cath::index::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_SCAFFOLD_VIEW_CACHE_INDEX_TAIL_HPP
