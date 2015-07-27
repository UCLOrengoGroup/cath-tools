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

#ifndef VIEW_CACHE_INDEX_TAIL_H_INCLUDED
#define VIEW_CACHE_INDEX_TAIL_H_INCLUDED

#include <boost/range/algorithm_ext/for_each.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>
#include <boost/tuple/tuple.hpp>

#include "common/c++14/cbegin_cend.h"
#include "structure/view_cache/index/detail/vcie_match_criteria.h"
#include "structure/view_cache/index/view_cache_index_entry.h"

#include <iostream> // ***** TEMPORARY *****
#include <vector>

namespace cath {
	namespace index {
		namespace detail {

			/// \brief TODOCUMENT
			class view_cache_index_tail final {
			private:
				/// \brief TODOCUMENT
				view_cache_index_entry_vec entries;

			public:
				void store(const view_cache_index_entry &);

				void store(const view_cache_index_entry &,
				           const boost::tuples::null_type &);

				bool empty() const;
				size_t get_num_cells() const;
				size_t dims_remaining() const;

				template <typename ACTN>
				void perform_action_on_all_match_at_leaves(ACTN &) const;

				template <typename ACTN>
				void perform_action_on_matches(const view_cache_index_entry &,
				                               const detail::vcie_match_criteria &,
				                               ACTN &) const;

				/// \brief TODOCUMENT
				template <typename ACTN>
				void perform_action_on_all_match_at_nodes(const view_cache_index_tail &,
				                                          const detail::vcie_match_criteria &,
				                                          ACTN &) const;

				/// \brief TODOCUMENT
				constexpr static size_t num_dims_remaining = 0;
			};

			/// \brief TODOCUMENT
			inline void view_cache_index_tail::store(const view_cache_index_entry &arg_entry ///< TODOCUMENT
			                                         ) {
				entries.push_back( arg_entry );
			}

			/// \brief TODOCUMENT
			inline void view_cache_index_tail::store(const view_cache_index_entry   &arg_entry,  ///< TODOCUMENT
			                                         const boost::tuples::null_type &/*arg_dfs*/ ///< TODOCUMENT
			                                         ) {
				store( arg_entry );
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
			void view_cache_index_tail::perform_action_on_all_match_at_leaves(ACTN &arg_action ///< TODOCUMENT
			                                                                  ) const {
				for (const view_cache_index_entry &entry : entries) {
					arg_action( entry );
				}
			}

			template <typename ACTN>
			void view_cache_index_tail::perform_action_on_matches(const view_cache_index_entry      &arg_entry,    ///< TODOCUMENT
			                                                      const detail::vcie_match_criteria &arg_criteria, ///< TODOCUMENT
			                                                      ACTN                              &arg_action    ///< TODOCUMENT
			                                                      ) const {
				for (const view_cache_index_entry &entry : entries) {
					if ( arg_criteria( arg_entry, entry ) ) {
						// const double sq_dist = squared_distance( arg_entry, entry );
						// const double score   = ( 10.0 / ( sq_dist + 10.0 ) );
	
						// std::cerr << "index_match ";
						// std::cerr << std::right << std::setw( 4 ) << arg_entry.get_from_index();
						// std::cerr << "";
						// std::cerr << std::right << std::setw( 4 ) << arg_entry.get_to_index();
						// std::cerr << " ";
						// std::cerr << std::right << std::setw( 4 ) << entry.get_from_index();
						// std::cerr << " ";
						// std::cerr << std::right << std::setw( 4 ) << entry.get_to_index();
						// std::cerr << " " << score;
						// std::cerr << std::endl;
						arg_action( arg_entry, entry );
					}
				}
			}

			/// \brief TODOCUMENT
			template <typename ACTN>
			void view_cache_index_tail::perform_action_on_all_match_at_nodes(const view_cache_index_tail       &arg_match_tail, ///< TODOCUMENT
			                                                                 const detail::vcie_match_criteria &arg_criteria,   ///< TODOCUMENT
			                                                                 ACTN                              &arg_action      ///< TODOCUMENT
			                                                                 ) const {
				// std::cerr << "Comparing "           << entries.size()
				//           << " query entries with " << arg_match_tail.entries.size()
				//           << " match entries"       << "\n";

				for (const view_cache_index_entry &x : entries) {
					for (const view_cache_index_entry &y : arg_match_tail.entries) {
						if ( arg_criteria( x, y ) ) {
							arg_action( x, y );
						}
					}
				}
			}

		}
	}
}

#endif
