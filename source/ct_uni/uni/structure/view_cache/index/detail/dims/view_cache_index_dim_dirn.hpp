/// \file
/// \brief The view_cache_index_dim_dirn class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_DIRN_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_DIRN_H

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/begin.hpp>

#include "common/debug_numeric_cast.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/not_implemented_exception.hpp"
#include "structure/view_cache/index/detail/vcie_match_criteria.hpp"

#include <cstddef>
#include <type_traits>
#include <vector>

namespace cath {
	namespace index {
		namespace detail {

			/// \brief TODOCUMENT
			class view_cache_index_dim_dirn final {
			private:
				static bool get_increases(const view_cache_index_entry &);

				/// \brief const-agnostic implementation of cell_at_value()
				///
				/// See GSL rule: Pro.Type.3: Don't use const_cast to cast away const (i.e., at all)
				/// (https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Pro-type-constcast)
				template <typename CELLS, typename DimDirn>
				static auto cell_at_value_impl(DimDirn      &/*prm_dim_dirn*/, ///< TODOCUMENT
				                               CELLS        &prm_cells,        ///< TODOCUMENT
				                               const bool   &prm_increases     ///< TODOCUMENT
				                               ) -> std::conditional_t<
				                                    	std::is_const< DimDirn >::value,
				                                    	const typename CELLS::value_type &,
				                                    	typename CELLS::value_type &
				                                    > {
#ifndef NDEBUG
					if ( prm_cells.empty() ) {
						BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot get entry at_value() with no populated cells"));
					}
#endif
					return prm_increases ? prm_cells.back()
					                     : prm_cells.front();
				}

				template <typename CELLS>
				typename CELLS::value_type & cell_at_value(CELLS &,
				                                           const typename CELLS::value_type &,
				                                           const bool &);

				template <typename CELLS>
				const typename CELLS::value_type & cell_at_value(const CELLS &,
				                                                 const bool &) const;

			public:
				template <typename CELLS, typename DEFAULTS>
				void store(const view_cache_index_entry &,
				           CELLS &,
				           const DEFAULTS &);

				template <typename CELLS, typename ACTN>
				void perform_action_on_matches(const view_cache_index_entry &,
				                               const CELLS &,
				                               const detail::vcie_match_criteria &,
				                               ACTN &) const;

				template <typename CELLS, typename ACTN>
				inline void perform_action_on_all_match_at_nodes(const CELLS &,
				                                                 const view_cache_index_dim_dirn &,
				                                                 const CELLS &,
				                                                 const detail::vcie_match_criteria &,
				                                                 ACTN &) const;
			};

			/// \brief Static function to return
			///
			/// \pre The view_cache_index_entry's from_index must not equal its to_index
			///      else an invalid_argument_exception will be thrown (unless NDEBUG is set)
			inline bool view_cache_index_dim_dirn::get_increases(const view_cache_index_entry &prm_entry ///< TODOCUMENT
			                                                     ) {
				// Grab the from_index and to_index from the view_cache_index_entry
				const size_t &from_index = prm_entry.get_from_index();
				const size_t &to_index   = prm_entry.get_to_index();
#ifndef NDEBUG
				// If in debug mode, check that the from index and to index aren't equal
				if ( from_index == to_index ) {
					BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot handle view_cache_index_entry for which the from index equals the to index"));
				}
#endif
				// Return whether the to_index is greater than the from_index
				return ( to_index > from_index );
			}

			

			/// \brief TODOCUMENT
			template <typename CELLS>
			inline typename CELLS::value_type & view_cache_index_dim_dirn::cell_at_value(CELLS                            &prm_cells,        ///< TODOCUMENT
			                                                                             const typename CELLS::value_type &prm_default_cell, ///< TODOCUMENT
			                                                                             const bool                       &prm_increases     ///< TODOCUMENT
			                                                                             ) {
				// If the cells haven't already been initialised, then create two cells
				// (one for entries with increasing indices; one for entries with decreasing indices)
				if ( prm_cells.empty() ) {
					prm_cells.assign( 2, prm_default_cell );
				}

				return cell_at_value_impl( *this, prm_cells, prm_increases );
			}

			/// \brief TODOCUMENT
			template <typename CELLS>
			inline const typename CELLS::value_type & view_cache_index_dim_dirn::cell_at_value(const CELLS  &prm_cells,    ///< TODOCUMENT
			                                                                                   const bool   &prm_increases ///< TODOCUMENT
			                                                                                   ) const {
				return cell_at_value_impl( *this, prm_cells, prm_increases );
			}

			/// \brief TODOCUMENT
			template <typename CELLS, typename DEFAULTS>
			inline void view_cache_index_dim_dirn::store(const view_cache_index_entry &prm_entry,   ///< TODOCUMENT
			                                             CELLS                        &prm_cells,   ///< TODOCUMENT
			                                             const DEFAULTS               &prm_defaults ///< TODOCUMENT
			                                             ) {
				const bool increases = get_increases( prm_entry );
				cell_at_value( prm_cells, typename CELLS::value_type{ prm_defaults.get_head() }, increases ).store(
					prm_entry,
					prm_defaults.get_tail()
				);
			}

			/// \brief TODOCUMENT
			template <typename CELLS, typename ACTN>
			inline void view_cache_index_dim_dirn::perform_action_on_matches(const view_cache_index_entry      &prm_entry,    ///< TODOCUMENT
			                                                                 const CELLS                       &prm_cells,    ///< TODOCUMENT
			                                                                 const detail::vcie_match_criteria &prm_criteria, ///< TODOCUMENT
			                                                                 ACTN                              &prm_action    ///< TODOCUMENT
			                                                                 ) const {
				const bool increases = get_increases( prm_entry );
				cell_at_value( prm_cells, increases ).perform_action_on_matches( prm_entry, prm_criteria, prm_action );
				if ( ! prm_criteria.get_require_matching_directions() ) {
					const bool decreases = ! increases;
					cell_at_value( prm_cells, decreases ).perform_action_on_matches( prm_entry, prm_criteria, prm_action );
				}
			}

			/// \brief TODOCUMENT
			template <typename CELLS, typename ACTN>
			inline void view_cache_index_dim_dirn::perform_action_on_all_match_at_nodes(const CELLS                       &prm_query_cells, ///< TODOCUMENT
			                                                                            const view_cache_index_dim_dirn   &prm_query_dim,   ///< TODOCUMENT
			                                                                            const CELLS                       &prm_match_cells, ///< TODOCUMENT
			                                                                            const detail::vcie_match_criteria &prm_criteria,    ///< TODOCUMENT
			                                                                            ACTN                              &prm_action       ///< TODOCUMENT
			                                                                            ) const {
				if ( ! prm_query_cells.empty() && ! prm_match_cells.empty() ) {
					for (const bool &increases : { false, true } ) {
						const typename CELLS::value_type &query_cell =               cell_at_value( prm_query_cells, increases );
						const typename CELLS::value_type &match_cell = prm_query_dim.cell_at_value( prm_match_cells, increases );
						if ( ! query_cell.empty() && !match_cell.empty() ) {
							std::cerr << "In dimension dirn, scanning for "
							          << ( increases ? "increases" : "decreases" )
							          << ". query_cell.size_" << query_cell.get_num_cells()
							          << ", match_cell.size_" << match_cell.get_num_cells()
							          << "\n";
							query_cell.perform_action_on_all_match_at_nodes( match_cell, prm_criteria, prm_action );
						}
					}
				}
			}

		} // namespace detail
	} // namespace index
} // namespace cath

#endif

