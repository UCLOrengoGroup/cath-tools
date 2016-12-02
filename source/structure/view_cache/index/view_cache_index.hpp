/// \file
/// \brief The view_cache_index class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_INDEX_VIEW_CACHE_INDEX_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_INDEX_VIEW_CACHE_INDEX_H

#include <boost/optional/optional_io.hpp>
#include <boost/optional.hpp>

#include "structure/view_cache/index/detail/pair_scan_action.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_dirn.hpp"
#include "structure/view_cache/index/detail/scaffold/view_cache_index_layer.hpp"
#include "structure/view_cache/index/detail/scaffold/view_cache_index_tail.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_linear_from_phi.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_linear_from_psi.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_linear_to_phi.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_linear_to_psi.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_linear_x.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_linear_y.hpp"
#include "structure/view_cache/index/detail/dims/view_cache_index_dim_linear_z.hpp"
#include "structure/view_cache/index/view_cache_index_entry.hpp"

#include <chrono>

namespace cath { namespace align { class alignment; } }
namespace cath { class protein; }
namespace cath { namespace index { class quad_find_action; } }
namespace cath { namespace index { class quad_find_action_check; } }

namespace cath {
	namespace index {
		namespace detail {

			/// \brief TODOCUMENT
			using standard_vci_nested_layers = 
			                                   view_cache_index_layer< view_cache_index_dim_dirn,
			                                   view_cache_index_layer< view_cache_index_dim_linear_from_phi,
			                                   // view_cache_index_layer< view_cache_index_dim_linear_from_psi,
			                                   // view_cache_index_layer< view_cache_index_dim_linear_to_phi,
			                                   // view_cache_index_layer< view_cache_index_dim_linear_to_psi,
			                                   view_cache_index_layer< view_cache_index_dim_linear_x,
			                                   view_cache_index_layer< view_cache_index_dim_linear_y,
			                                   view_cache_index_layer< view_cache_index_dim_linear_z,
			                                   
			                                   view_cache_index_tail
			                                   > > > > >;
		} // namespace detail


		/// \brief TODOCUMENT
		class view_cache_index final {
		public:
			/// \brief TODOCUMENT
			using dim_tuple = boost::tuple<
			                                detail::view_cache_index_dim_dirn,
			                                detail::view_cache_index_dim_linear_from_phi,
			                                // detail::view_cache_index_dim_linear_from_psi,
			                                // detail::view_cache_index_dim_linear_to_phi,
			                                // detail::view_cache_index_dim_linear_to_psi,
			                                detail::view_cache_index_dim_linear_x,
			                                detail::view_cache_index_dim_linear_y,
			                                detail::view_cache_index_dim_linear_z,
			                                detail::view_cache_index_tail> ;
		private:

			/// \brief TODOCUMENT
			dim_tuple dim_defaults;

			/// \brief TODOCUMENT
			detail::standard_vci_nested_layers the_index;

		public:
			explicit view_cache_index(const dim_tuple &);

			void store(const view_cache_index_entry &);

			template <typename ACTN>
			void perform_action_on_matches(const view_cache_index_entry &,
			                               const detail::vcie_match_criteria &,
			                               ACTN &) const;

			template <typename ACTN>
			void perform_action_on_all_match_at_leaves(const view_cache_index &,
			                                           const detail::vcie_match_criteria &,
			                                           ACTN &) const;

			template <typename ACTN>
			void perform_action_on_all_match_at_nodes(const view_cache_index &,
			                                          const detail::vcie_match_criteria &,
			                                          ACTN &) const;
		};

		/// \brief TODOCUMENT
		template <typename ACTN>
		void view_cache_index::perform_action_on_matches(const view_cache_index_entry      &arg_entry,    ///< TODOCUMENT
		                                                 const detail::vcie_match_criteria &arg_criteria, ///< TODOCUMENT
		                                                 ACTN                              &arg_action    ///< TODOCUMENT
		                                                 ) const {
			the_index.perform_action_on_matches( arg_entry, arg_criteria, arg_action );
		}

		/// \brief TODOCUMENT
		template <typename ACTN>
		void view_cache_index::perform_action_on_all_match_at_leaves(const view_cache_index            &arg_search_index, ///< TODOCUMENT
		                                                             const detail::vcie_match_criteria &arg_criteria,     ///< TODOCUMENT
		                                                             ACTN                              &arg_action        ///< TODOCUMENT
		                                                             ) const {
			const detail::pair_scan_action<view_cache_index> the_pair_scan_action( arg_search_index, arg_criteria, arg_action );
			the_index.perform_action_on_all_match_at_leaves( the_pair_scan_action );
		}

		/// \brief TODOCUMENT
		template <typename ACTN>
		void view_cache_index::perform_action_on_all_match_at_nodes(const view_cache_index            &arg_search_index, ///< TODOCUMENT
		                                                            const detail::vcie_match_criteria &arg_criteria,     ///< TODOCUMENT
		                                                            ACTN                              &arg_action        ///< TODOCUMENT
		                                                            ) const {
			the_index.perform_action_on_all_match_at_nodes( arg_search_index.the_index, arg_criteria, arg_action );
		}

		view_cache_index build_view_cache_index(const double &,
		                                        const detail::angle_type &,
		                                        const detail::angle_type &,
		                                        const protein &,
		                                        const detail::vcie_match_criteria &);

		std::chrono::high_resolution_clock::duration process_quads_indexed(const protein &,
		                                                                   const protein &,
		                                                                   const double &,
		                                                                   const detail::angle_type &,
		                                                                   const detail::angle_type &,
		                                                                   const detail::vcie_match_criteria &,
		                                                                   quad_find_action_check &);

		std::chrono::high_resolution_clock::duration process_quads_complete(const protein &,
		                                                                    const protein &,
		                                                                    const double &,
		                                                                    const detail::vcie_match_criteria &,
		                                                                    quad_find_action &);
	} // namespace index
} // namespace cath

#endif
